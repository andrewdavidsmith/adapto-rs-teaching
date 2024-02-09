/* MIT License
 *
 * Copyright (c) 2023-2024 Andrew Smith
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

/// Program to cut adaptors from sequenced reads. Accepts one adaptor
/// and will apply it to both ends in paired-end data. Removes Ns at
/// the end of reads. Removes low quality bases at ends of reads.
/// Output is compressed as bgzf. Input may be compressed as gz/bgzf
/// or not. Extra threads help with compressing output and
/// decompressing input.
use clap::Parser;
use file_format::FileFormat;
use num_cpus;
use std::error::Error;
use std::str::from_utf8;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Fastq input file
    #[structopt(required = true)]
    fastq: String,

    /// Paired-end input second fastq file
    #[structopt(required = false)]
    pfastq: Option<String>,

    /// Output file
    #[arg(short, long)]
    out: String,

    /// Second output file for paired-end reads
    #[structopt(required = false)]
    #[arg(short, long)]
    pout: Option<String>,

    /// Quality score cutoff
    #[arg(short, long, default_value_t = 20)]
    qual_cutoff: u8,

    /// Adaptor sequence
    #[arg(short, long, default_value = "AGATCGGAAGAGC")]
    adaptor: Option<String>,

    /// Keep all read prefixes (not implemented)
    #[arg(short, long, default_value_t = true)]
    keep_prefix: bool,

    /// Zip output files as BGZF format
    #[arg(short, long)]
    zip: bool,

    /// Threads to use
    #[arg(short, long, default_value_t = 1)]
    threads: u32,

    /// Buffer size for reading input
    #[arg(short, long, default_value_t = 256*1024)]
    buffer_size: usize,

    /// Be verbose
    #[arg(short, long)]
    verbose: bool,
}

fn is_readable(filename: &String) -> bool {
    use std::fs::File;
    let mut f = match File::open(&filename) {
        Ok(file) => file,
        _ => return false,
    };
    let mut byte = [0_u8];
    use std::io::Read;
    match f.read_exact(&mut byte) {
        Ok(_) => true,
        Err(_) => false,
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    if args.threads <= 0 {
        return Err("number of threads must be positive")?;
    }

    if args.buffer_size <= 0 {
        return Err("buffer size must be positive")?;
    }

    let adaptor = args.adaptor.unwrap().into_bytes();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads as usize)
        .build_global()
        .unwrap();

    if args.verbose {
        eprintln!("input file: {}", args.fastq);
        eprintln!("input file format: {}", FileFormat::from_file(&args.fastq)?);
        eprintln!("output file: {}", args.out);
        eprintln!("quality score cutoff: {}", args.qual_cutoff);
        eprintln!("adaptor sequence: {}", from_utf8(&adaptor)?);
        eprintln!("keep prefix: {}", args.keep_prefix);
        eprintln!("compress output: {}", args.zip);
        eprintln!("threads requested: {}", args.threads);
        eprintln!("detected cpu cores: {}", num_cpus::get());
        eprintln!("buffer size: {}", args.buffer_size);
        match (&args.pfastq, &args.pout) {
            (Some(x), Some(y)) => {
                eprintln!("input2 file: {}", x);
                eprintln!("input2 file format: {}", FileFormat::from_file(&x)?);
                eprintln!("output2 file: {}", y);
            }
            (Some(_), None) | (None, Some(_)) => {
                Err("paired end requires two input and output files")?;
            }
            (None, None) => (),
        }
    }

    // ADS: do this 1st so we don't waste time on end2 if end1 is bad
    if !is_readable(&args.fastq) {
        return Err(format!("input file not readable: {}", args.fastq))?;
    }

    use adapto_rs::remove_adaptors;

    if let (Some(pfastq), Some(pout)) = (args.pfastq, args.pout) {
        if !is_readable(&pfastq) {
            return Err(format!("input file not readable: {}", pfastq))?;
        }
        remove_adaptors(
            args.zip,
            args.threads,
            args.buffer_size,
            &adaptor,
            &pfastq,
            &pout,
            args.qual_cutoff,
        )?;
    }

    remove_adaptors(
        args.zip,
        args.threads,
        args.buffer_size,
        &adaptor,
        &args.fastq,
        &args.out,
        args.qual_cutoff,
    )
}
