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


/// Program to cut adaptors from sequenced reads. Accepts one adaptor
/// and will apply it to both ends in paired-end data. Removes Ns at
/// the end of reads. Removes low quality bases at ends of reads.
/// Output is compressed as bgzf. Input may be compressed as gz/bgzf
/// or not. Extra threads help with compressing output and
/// decompressing input.


use clap::Parser;
use std::process::ExitCode;
extern crate num_cpus;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Fastq input file
    #[structopt(required=true)]
    fastq: String,

    /// Paired-end input second fastq file
    #[structopt(required=false)]
    pfastq: Option<String>,

    /// Output file
    #[arg(short, long)]
    out: String,

    /// Second output file for paired-end reads
    #[structopt(required=false)]
    #[arg(short, long)]
    pout: Option<String>,

    /// Quality score cutoff
    #[arg(short, long, default_value_t = 20)]
    qual_cutoff: u8,

    /// Adaptor sequence
    #[arg(short, long, default_value = "AGATCGGAAGAGC")]
    adaptor: Option<String>,

    /// Adaptor sequence
    #[arg(short, long)]
    keep_header: bool,

    /// Zip the output (does nothing...)
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


fn main() -> ExitCode {

    let args = Args::parse();

    if args.threads <= 0 {
        eprintln!("number of threads must be positive");
        return ExitCode::FAILURE;
    }

    if args.buffer_size <= 0 {
        eprintln!("buffer size must be positive");
        return ExitCode::FAILURE;
    }

    let adaptor = args.adaptor.unwrap().into_bytes();

    if args.verbose {
        eprintln!("input file: {}", args.fastq);
        eprintln!("output file: {}", args.out);
        eprintln!("quality score cutoff: {}", args.qual_cutoff);
        use std::str::from_utf8;
        eprintln!("adaptor sequence: {}", from_utf8(&adaptor).unwrap());
        eprintln!("keep header: {}", if args.keep_header {"no"} else {"no"});
        eprintln!("compress output: {}", if args.zip {"yes"} else {"yes"});
        eprintln!("threads requested: {}", args.threads);
        eprintln!("detected cpu cores: {}", num_cpus::get());
        eprintln!("buffer size: {}", args.buffer_size);
        match (&args.pfastq, &args.pout) {
            (Some(x), Some(y)) => {
                eprintln!("input file2: {}", x);
                eprintln!("output file2: {}", y);
            },
            (Some(_), None) | (None, Some(_)) => {
                eprintln!("paired end requires two input and output files");
                return ExitCode::FAILURE;
            },
            (None, None) => {}
        }
    }

    use adaptrs::process_reads;
    if let (Some(pfastq), Some(pout)) = (args.pfastq, args.pout) {
        process_reads(args.zip, args.threads, args.buffer_size, &adaptor,
                      &pfastq, &pout, args.qual_cutoff).unwrap();
    }
    process_reads(args.zip, args.threads, args.buffer_size, &adaptor,
                  &args.fastq, &args.out, args.qual_cutoff).unwrap();

    ExitCode::SUCCESS
}
