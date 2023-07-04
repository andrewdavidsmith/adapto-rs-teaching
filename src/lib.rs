/* MIT License
 *
 * Copyright (c) 2023 Andrew Smith
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

use std::io::{Write,Read};
use std::process;
use std::cmp::{max, min};


// ADS: this works well, but imports way more than is
// needed. Considering `gzp`
use rust_htslib::bgzf;
use rust_htslib::tpool::ThreadPool;


/// The prefix function for the KMP algorithm
fn kmp_prefix_function(p: &[u8]) -> Vec<usize> {
    let n = p.len();
    let mut sp = vec![0 as usize; n];
    let mut k = 0usize;
    for i in 1..n {
        while k > 0 && p[k] != p[i] {
            k = sp[k - 1];
        }
        if p[k] == p[i] {
            k += 1;
        }
        sp[i] = k;
    }
    sp
}


/// The KMP algorithm that returns the first full match or the start
/// of any suffix match to the pattern (i.e. adaptor).
fn kmp(adaptor: &[u8], sp: &[usize], read: &[u8], m: usize) -> usize {
    let n = adaptor.len();
    let mut j: usize = 0;
    let mut i: usize = 0;
    while i < m {
        // look for the longest prefix of P that is the same as a
        // suffix of P[1..j - 1] AND has a different next character
        while j > 0 && adaptor[j] != read[i] {
            j = sp[j - 1];
        }
        // check if the character matches
        if adaptor[j] == read[i] {
            j += 1;
        }
        // if we have already successfully compared all positions in
        // P, then we have found a match
        if j == n {
            return i - n + 1;
        }
        i += 1;
    }
    // if we have not found a full match, then return the maximum
    // prefix match of the pattern
    i - j
}


/// Find the positions in the read of the first non-N and last non-N.
fn trim_n_ends(read: &[u8]) -> (usize, usize) {
    let start = match read.iter().position(|&x| x != b'N') {
        Some(x) => x,
        _ => 0,
    };
    let stop = match read.iter().rposition(|&x| x != b'N') {
        Some(x) => x + 1,
        _ => 0,
    };
    (start, stop)
}


/// Find the positions in the read where quality scores indicate the
/// read should be trimmed. This is copied from cutadapt source.
fn qual_trim(qual: &[u8], cut_front: i32, cut_back: i32)
             -> (usize, usize) {

    const QUAL_BASE: i32 = 33; // assumes base quality starts at 33

    /* ADS: COPIED FROM cutadapt SOURCE */
    let n = qual.len();

    //  find trim position for 5' end
    let mut start: usize = 0;
    let mut s: i32 = 0;
    let mut max_qual: i32 = 0;

    if cut_front > 0 {
        let cut_front = cut_front + QUAL_BASE;
        for i in 0..n {
            s += (cut_front + QUAL_BASE) - qual[i] as i32;
            if s < 0 {
                break;
            }
            if s > max_qual {
                max_qual = s;
                start = i + 1;
            }
        }
    }
    // same for 3' end
    let mut stop: usize = n;
    max_qual = 0;
    s = 0;
    let cut_back = cut_back + QUAL_BASE;
    for i in (0..n).rev() {
        s += cut_back - qual[i] as i32;
        if s < 0 {
            break;
        }
        if s > max_qual {
            max_qual = s;
            stop = i;
        }
    }
    if start >= stop {
        (start, stop) = (0, 0)
    }
    (start as usize, stop as usize)
}


/// FQBuf is just a buffer that keeps a cursor (position) and how much
/// of the buffer is filled.
#[derive(Default)]
struct FQBuf {
    buf: Vec<u8>,
    pos: usize, // current offset into buf, must always be <= `filled`.
    filled: usize,
}


impl FQBuf {
    fn shift(&mut self) {
        let mut j = 0;
        for i in self.pos..self.filled {
            self.buf[j] = self.buf[i];
            j += 1;
        }
        self.filled = j;
        self.pos = 0;
    }
    fn next_newline(&self, offset: usize) -> usize {
        for i in offset..self.filled {
            if self.buf[i] == b'\n' {
                return i;
            }
        }
        0
    }
}


/// FQRec is a FASTQ record that represents the position of the start
/// of the name (n), the start of the read sequence (r), the start of
/// the other name, the one with the "+" (o), and the start of the
/// quality scores (q). The `start` and `stop` variables are used to
/// store the offsets of trimmed ends for the read and quality scores
/// strings.
#[derive(Default)]
struct FQRec {
    n: usize, // start of "name"
    r: usize, // start of "read"
    o: usize, // start of "other"
    q: usize, // start of "quality" scores
    start: usize, // *start* of good part of seq
    stop: usize, // *stop* of good part of seq
}


impl FQRec {
    fn set_start_stop(&mut self,
                      adaptor: &[u8],
                      sp: &Vec<usize>, cutoff: u8, buf: &FQBuf) {
        let x = buf.pos - 1;
        let (qstart, qstop) = qual_trim(&buf.buf[self.q..x], 0, cutoff as i32);
        // consecutive N values at both ends
        let x = self.o - 1;
        let (nstart, nstop) = trim_n_ends(&buf.buf[self.r..x]);
        // do not allow any N or low qual bases to interfere with adaptor
        self.stop = min(qstop, nstop);
        // find the adaptor at the 3' end
        let adaptor_start = kmp(adaptor, &sp, &buf.buf[self.r..x], self.stop);
        self.stop = min(self.stop, adaptor_start);
        let x = self.r + self.stop;
        let (_, nstop) = trim_n_ends(&buf.buf[self.r..x]);
        self.stop = min(self.stop, nstop);
        self.start = min(max(qstart, nstart), self.stop);
    }
    fn find_record_from_buf(&mut self, buf: &mut FQBuf) -> bool {
        // ADS: does not fail if a record is broken, but will ignore it
        self.r = buf.next_newline(buf.pos) + 1;
        self.o = buf.next_newline(self.r) + 1;
        self.q = buf.next_newline(self.o) + 1;
        self.n = buf.next_newline(self.q) + 1; // to swap later!!!
        if self.r > 1 && self.o > 1 && self.q > 1 && self.n > 1 {
            (self.n, buf.pos) = (buf.pos, self.n);
            return true;
        }
        else {
            return false;
        }
    }
    fn write<W: Write>(&mut self, buf: &mut FQBuf, writer: &mut W) {
        buf.buf[self.o + 1] = b'\n';
        buf.buf[self.r + self.stop] = b'\n';
        buf.buf[self.q + self.stop] = b'\n';
        self.stop += 1;
        /* ADS: below, this will compile with bgzf, but it segfaults */
        // use std::io::IoSlice;
        // writer.write_vectored(
        //     &[IoSlice::new(&buf.buf[self.n..self.r]),
        //       IoSlice::new(&buf.buf[(self.r + self.start)..
        //                             (self.r + self.stop)]),
        //       IoSlice::new(&buf.buf[self.o..(self.o + 2)]),
        //       IoSlice::new(&buf.buf[(self.q + self.start)..
        //                             (self.q + self.stop)]),
        //     ]).unwrap();
        writer.write(&buf.buf[self.n..self.r]).unwrap();
        writer.write(&buf.buf[(self.r + self.start)..(self.r + self.stop)]).unwrap();
        writer.write(&buf.buf[self.o..(self.o + 2)]).unwrap();
        writer.write(&buf.buf[(self.q + self.start)..(self.q + self.stop)]).unwrap();
    }
}


pub fn process_reads(
    _zip: bool,
    n_threads: u32,
    buffer_size: usize,
    adaptor: &[u8],
    input: &String,
    output: &String,
    cutoff: u8
) -> Result<(), std::io::Error> {

    let sp = kmp_prefix_function(adaptor);

    // open the bgzf files for reading and writing
    let mut reader = bgzf::Reader::from_path(input).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });
    let mut writer = bgzf::Writer::from_path(output).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });

    // make a thread pool and give it to the input and output files
    let tpool = match ThreadPool::new(n_threads) {
        Ok(p) => p,
        Err(error) => {
            eprintln!("failed to acquire threads: {error}");
            process::exit(1);
        }
    };
    reader.set_thread_pool(&tpool).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });
    writer.set_thread_pool(&tpool).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });

    let mut buf: FQBuf = Default::default();
    buf.buf.resize(buffer_size, b'\0');

    let mut fq: FQRec = Default::default();

    buf.filled += reader.read(&mut buf.buf[buf.filled..]).unwrap();

    let mut found_rec = false;
    let mut data_exhausted = false;

    loop {
        loop {
            found_rec = fq.find_record_from_buf(&mut buf);
            if found_rec || data_exhausted {
                break;
            }
            buf.shift(); // reload the buffer; not circular...
            buf.filled += reader.read(&mut buf.buf[buf.filled..]).unwrap();
            data_exhausted = buf.filled < buf.buf.len();
        }
        if !found_rec && data_exhausted {
            break;
        }
        fq.set_start_stop(&adaptor, &sp, cutoff, &buf);
        fq.write(&mut buf, &mut writer);
    }
    writer.flush().unwrap();

    Ok(())
}
