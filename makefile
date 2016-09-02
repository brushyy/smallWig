all: parameters.h arithmetic.h encoder_by_blocks.cpp encoder.cpp getseqblock.c getmatlabseq.c getseq.c fmergeadd.c fmergesub.c fmerge.c prll.c fmergecount.c lpaq1.cpp wig2smallwig.cpp parameters.h arithmetic.h decoder_by_blocks.cpp decoder.cpp partial_decoder_by_blocks.cpp to_wig.c to_wig_by_blocks.c to_wig_by_blocks_partial.c to_wig_sub.c fsplit.c fmergesub.c fremove.c lpaq1.cpp prll.c smallwig2wig.cpp
	g++ -o  wig2smallwig wig2smallwig.cpp -lm -g -w
	g++ -o  smallwig2wig smallwig2wig.cpp -lm -g -w
	

