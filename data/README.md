This folder is used to store hash files. The format of a hash file sretty
simple. Each file is the concatenation of 64-bit integers, each integer
represent a hash value.

CMake will download the files specified in `CMakeLists.txt`. The default file
is `test.hash` which contains 12.5 million 64-bit hash values.
You can also get the files of the SIGIR 2016 paper by removing the hash symbol
at the beginning of line 6 in `CMakeLists.txt`.
