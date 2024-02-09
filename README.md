# adapto-rs

The code in this repo is for teaching. I decided to use this repo for
some code I wanted to show students. I will add features to a
different adaptor trimming tool, but in this repo I will try to keep
the code as small as possible so I can point quickly to some important
parts.

`adapto-rs` removes adaptors from short-read sequencing data. I wrote
this in Rust initially because I wanted to experiment with the
libraries for parallelism. It turned out this application does not
need much in the way of parallelism. I think this code is reasonably
fast.

For students: I use the Knuth-Morris-Pratt algorithm to do the pattern
matching required for finding the exact match of the adaptor sequence
with each read. I think the approach is good, and it also works nicely
to automatically find the longest suffix-prefix match of read-adaptor.
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
