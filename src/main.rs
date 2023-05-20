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


use clap::Parser;
use std::process::ExitCode;

use radapt::process_reads;
use radapt::process_reads_pe;


#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// fastq file
    #[structopt(required=true)]
    fastq: String,

    /// second fastq file for paired
    #[structopt(required=false)]
    pfastq: Option<String>,

    /// output file
    #[arg(short, long)]
    out: String,

    /// second output file for paired
    #[structopt(required=false)]
    #[arg(short, long)]
    pout: Option<String>,

    /// quality score cutoff
    #[arg(short, long, default_value_t = 20)]
    qual_cutoff: u8,

    /// be verbose
    #[arg(short, long)]
    verbose: bool,
}


fn main() -> ExitCode {

    let args = Args::parse();

    if let Some(pfastq) = args.pfastq {
        if let Some(pout) = args.pout {
            process_reads_pe(args.verbose, &args.fastq, &pfastq,
                             &args.out, &pout, args.qual_cutoff).unwrap();
        }
        else {
            eprintln!("specifying two inputs requires two outputs");
            return ExitCode::FAILURE;
        }
    }
    else {
        process_reads(args.verbose, &args.fastq,
                      &args.out, args.qual_cutoff).unwrap();
    }
    ExitCode::SUCCESS
}
