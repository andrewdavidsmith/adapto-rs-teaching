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


fn kmp(adaptor: &[u8], read: &Vec<u8>, sp: &[usize]) -> usize {

    let n = adaptor.len();
    let m = read.len();

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


fn trim_n_ends(read: &Vec<u8>) -> (usize, usize) {
    (read.iter().position(|&x| x != b'N').unwrap(),
     read.iter().rposition(|&x| x != b'N').unwrap() + 1)
}


fn trim_qual_ends(qual: &Vec<u8>, cutoff: u8) -> (usize, usize) {
    let start = match qual.iter().position(|&x| x >= cutoff) {
        Some(x) => x,
        _ => 0,
    };
    let stop = match qual.iter().rposition(|&x| x >= cutoff) {
        Some(x) => x + 1,
        _ => 0,
    };
    (start, stop)
}


fn start_stop(sp: &Vec<usize>, read: &Vec<u8>, qual: &Vec<u8>, cutoff: u8)
              -> (usize, usize) {
    // quality score cutoff positions at both ends
    let (qstart, qstop) = trim_qual_ends(&qual, cutoff + QUAL_BASE);
    // consecutive N values at both ends
    let (nstart, nstop) = trim_n_ends(&read);
    // find the adaptor at the 3' end
    let a = kmp(ADAPTOR, &read, &sp);
    (max(qstart, nstart), min(min(qstop, nstop), a))
}


fn report<W: Write>(
    out: &mut BufWriter<W>,
    the_name: &Vec<u8>,
    the_read: &Vec<u8>,
    the_qual: &Vec<u8>,
    start: usize,
    stop: usize,
) {
    out.write(&the_name).unwrap();
    out.write(&[b'\n']).unwrap();
    out.write(&the_read[start..stop]).unwrap();
    out.write(&[b'\n']).unwrap();
    out.write(&[b'+', b'\n']).unwrap();
    out.write(&the_qual[start..stop]).unwrap();
    out.write(&[b'\n']).unwrap();
}


pub fn process_reads(
    _verbose: bool,
    input: &String,
    output: &String,
    cutoff: u8,
) -> Result<(), std::io::Error> {

    let sp = kmp_prefix_function(ADAPTOR);

    // setup the input file
    let mut in_file = BufReader::new(File::open(input).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    }));

    // setup the output stream
    let mut out = BufWriter::new(File::create(output).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    }));

    let mut line_idx: usize = 0;

    // iterate over lines in the fastq file
    let mut read = vec![];
    let mut name = vec![];
    let mut buf = vec![];
    loop {

        let n_bytes = in_file.read_until(b'\n', &mut buf)
            .expect("failed reading fastq");
        if n_bytes == 0 || buf[n_bytes - 1] != b'\n' {
            break;
        }
        buf.truncate(n_bytes - 1);

        if line_idx % 4 == 0 {
            let t = buf.iter().position(|&x| x == b' ' || x == b'\t').unwrap();
            buf.truncate(t);
            name = buf.clone();
        }
        else if line_idx % 4 == 1 {
            read = buf.clone();
        }
        // do nothing for '+' line
        else if line_idx % 4 == 3 {
            let (a, b) = start_stop(&sp, &read, &buf, cutoff);
            if a < b {
                report(&mut out, &name, &read, &buf, a, b);
            }
        }
        line_idx += 1;
        buf.clear();
    }
    out.flush().unwrap();

    Ok(())
}


pub fn process_reads_pe(
    _verbose: bool,
    input1: &String,
    input2: &String,
    output1: &String,
    output2: &String,
    cutoff: u8,
) -> Result<(), std::io::Error> {

    let sp = kmp_prefix_function(ADAPTOR);

    // setup the first input stream
    let mut in_file1 = BufReader::new(File::open(input1).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    }));

    // setup the first output stream
    let mut out1 = BufWriter::new(File::create(output1).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    }));

    // setup the second input stream
    let mut in_file2 = BufReader::new(File::open(input2).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    }));

    // setup the second output stream
    let mut out2 = BufWriter::new(File::create(output2).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    }));

    let mut line_idx: usize = 0;

    // iterate over lines in the fastq file
    let mut read1 = vec![];
    let mut read2 = vec![];
    let mut name1 = vec![];
    let mut name2 = vec![];
    let mut buf1 = vec![];
    let mut buf2 = vec![];
    loop {

        let n_bytes1 = in_file1.read_until(b'\n', &mut buf1)
            .expect("failed reading end1 fastq");
        if n_bytes1 == 0 || buf1[n_bytes1 - 1] != b'\n' {
            break;
        }
        buf1.truncate(n_bytes1 - 1);

        let n_bytes2 = in_file2.read_until(b'\n', &mut buf2)
            .expect("failed reading end2 fastq");
        if n_bytes2 == 0 || buf2[n_bytes2 - 1] != b'\n' {
            break;
        }
        buf2.truncate(n_bytes2 - 1);

        if line_idx % 4 == 0 {
            let y = buf1.iter().position(|&x| x == b' ' || x == b'\t').unwrap();
            buf1.truncate(y);
            name1 = buf1.clone();

            let y = buf2.iter().position(|&x| x == b' ' || x == b'\t').unwrap();
            buf2.truncate(y);
            name2 = buf2.clone();
        }
        else if line_idx % 4 == 1 {
            read1 = buf1.clone();
            read2 = buf2.clone();
        }
        else if line_idx % 4 == 3 {

            let (mut a1, mut b1) = start_stop(&sp, &read1, &buf1, cutoff);
            let (mut a2, mut b2) = start_stop(&sp, &read2, &buf2, cutoff);

            let first_is_good = a1 < b1;
            let second_is_good = a2 < b2;

            if first_is_good || second_is_good {
                if !first_is_good {
                    (a1, b1) = (0, 1);
                    read1[0] = b'N';
                    buf1[0] = b'B';
                }
                report(&mut out1, &name1, &read1, &buf1, a1, b1);

                if !second_is_good {
                    (a2, b2) = (0, 1);
                    read2[0] = b'N';
                    buf2[0] = b'B';
                }
                report(&mut out2, &name2, &read2, &buf2, a2, b2);
            }
        }
        line_idx += 1;
        buf1.clear();
        buf2.clear();
    }
    out1.flush().unwrap();
    out2.flush().unwrap();

    Ok(())
}
