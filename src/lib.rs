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

use std::fs::File;
use std::io::{prelude::*,BufReader,BufWriter,Write};
use std::process;
use std::cmp::{max, min};


const ADAPTOR: &[u8] = b"AGATCGGAAGAGC";
const QUAL_BASE: u8 = 33;


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


fn qual_trim(qual: &[u8], cut_front: i32, cut_back: i32, base: i32)
             -> (usize, usize) {
    /* ADS: COPIED FROM cutadapt SOURCE */
    let n = qual.len();
    //  find trim position for 5' end
    let mut start: usize = 0;
    let mut s: i32 = 0;
    let mut max_qual: i32 = 0;
    if cut_front > 0 {
        for i in 0..n {
            s += cut_front - (qual[i] as i32 - base);
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
    for i in (0..n).rev() {
        s += cut_back - (qual[i] as i32 - base);
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


#[derive(Default)]
struct FQRec {
    data: Vec<u8>,
    n: usize, // end of "name"
    r: usize, // end of "read"
    o: usize, // end of "other"
    q: usize, // end of "quality" scores
    start: usize, // start of good part of seq
    stop: usize, // stop of good part of seq
}


impl FQRec {
    fn set_start_stop(&mut self, sp: &Vec<usize>, cutoff: u8) {
        // quality score cutoff positions at both ends
        let (qstart, qstop) = qual_trim(&self.data[self.o..(self.q - 1)], 0,
                                        cutoff as i32, QUAL_BASE as i32);
        // consecutive N values at both ends
        let (nstart, nstop) = trim_n_ends(&self.data[self.n..(self.r - 1)]);
        // do not allow any N or low qual bases to interfere with adaptor
        self.stop = min(qstop, nstop);
        // find the adaptor at the 3' end
        let adaptor_start =
            kmp(ADAPTOR, &sp, &self.data[self.n..(self.r - 1)], self.stop);
        self.stop = min(self.stop, adaptor_start);
        self.start = min(max(qstart, nstart), self.stop);
    }
    fn read_from<R: BufRead>(&mut self, rdr: &mut R) -> bool {
        // ADS: does not fail if a record is broken, but will ignore it
        self.data.clear();
        self.n = rdr.read_until(b'\n', &mut self.data).expect("read fail");
        self.r = self.n;
        self.r += rdr.read_until(b'\n', &mut self.data).expect("read fail");
        self.o = self.r;
        self.o += rdr.read_until(b'\n', &mut self.data).expect("read fail");
        self.q = self.o;
        self.q += rdr.read_until(b'\n', &mut self.data).expect("read fail");
        self.n != 0 && self.r != 0 && self.o != 0 && self.q != 0
    }
    fn write<W: Write>(&mut self, writer: &mut W) {
        self.data[self.r + 1] = b'\n';
        self.data[self.n + self.stop] = b'\n';
        self.data[self.o + self.stop] = b'\n';
        self.stop += 1;
        use std::io::IoSlice;
        writer.write_vectored(
            &[IoSlice::new(&self.data[..self.n]),
              IoSlice::new(&self.data[(self.n + self.start)..
                                      (self.n + self.stop)]),
              IoSlice::new(&self.data[self.r..self.r + 2]),
              IoSlice::new(&self.data[(self.o + self.start)..
                                      (self.o + self.stop)])]).unwrap();
    }
}


pub fn process_reads(input: &String, output: &String, cutoff: u8)
                     -> Result<(), std::io::Error> {

    const BUFFER_SIZE: usize = 128*1024;
    const FQR_BUFFER_SIZE: usize = 4096;

    let sp = kmp_prefix_function(ADAPTOR);

    // ADS: buffered reader and buffer within FQRec is redundant
    let mut reader =
        BufReader::with_capacity(BUFFER_SIZE,
                                 File::open(input).unwrap_or_else(|err| {
                                     eprintln!("{err}");
                                     process::exit(1);
                                 }));

    let mut writer =
        BufWriter::with_capacity(BUFFER_SIZE,
                                 File::create(output).unwrap_or_else(|err| {
                                     eprintln!("{err}");
                                     process::exit(1);
                                 }));

    let mut fq: FQRec = Default::default();
    fq.data.reserve(FQR_BUFFER_SIZE);
    // iterate over lines in the fastq file
    loop {
        if !fq.read_from(&mut reader) {
            break;
        }
        fq.set_start_stop(&sp, cutoff);
        fq.write(&mut writer);
    }
    writer.flush().unwrap();

    Ok(())
}
