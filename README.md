# adapto-rs

Remove adaptors from short-read sequencing data. I needed a faster
tool to remove adaptor sequences from reads for my own in my own data
analysis. I wrote this in Rust because I wanted to experiment with the
libraries for parallelism. It turns out this application does not need
much in the way of parallelism. I think this code is reasonably fast.

I use the Knuth-Morris-Pratt algorithm to do the pattern matching
required for finding the exact match of the adaptor sequence with each
read. I think the approach is good, and it also works nicely to
automatically find the longest suffix-prefix match of read-adaptor.
The KMP algorithm is basically the same thing as a compiled regular
expression, but tailored for a regular expression that is simply a
string.

If you have `cargo` installed, you can build this code by doing:
```
cargo build --release
```
in the root of the source directory, and then run
```
./target/release/adapto-rs
```
to see the command line arguments.
