The smallWig tool provides compression/decompression for WIG files. 
Copyright by Zhiying Wang, 2015.
Please contact zhiying@uci.edu for any questions and source files.

==============================================
											usage
==============================================											

To compress a file, please use:
 wig2smallwig [InputFile] [OutputFile] 
options: 
-m [N=0..9, uses 3 + 3*2^N MB memory, decompress should use same N] 
	context mixing 
-r [B, encode block size from 8 to 32]
	random access and encode by blocks of size 2^B 
-p [total number of processes] 
	parallel realization, only available if -r is enabled and -m is disabled


To decompress a file, please use:
To use:
 smallwig2wig [InputFile] [OutputFile] 
options: 
-s [ChrmName (e.g. chr1)] [Query Start (integer)] [Query End (integer)] 
	subsequence query 
-p [total number of processes]
	parallel realization, only available if -s is disabled and -r is enabled in encoding


==============================================
									  example
==============================================	
1. Compress a file using standard setup:
$ wig2smallwig in.wig out.swig

2. Compress a file allowing random query in the future:
$ wig2smallwig in.wig out.wig -r 9 13

3. Compress a file allowing random query in the future, and use 4 parallel processors:
$ wig2smallwig in.wig out.wig -r 9 13 -p 4

4. Compress a file as much as possible for archive purposes using the maximum amount of memory:
$ wig2smallwig in.wig out.wig -m 9

5. Decompress a whole file: 
$ smallwig2wig in.swig out.wig

6. Decompress a file only in chr1, start from location 300, end at location 500:
$ smallwig2wig in.swig out.wig -s chr1 300 500

7. Decompress a file with 4 parallel processors:
$ smallwig2wig in.swig out.wig -p 4
