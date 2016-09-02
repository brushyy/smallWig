//encode data by blocks of size 2^BlockSize using arithmetic code
//also supports parallel implementation
//if nproc<=0, then run whole file
//if nproc>=, find corresponding sub files and encode




int encoder_by_blocks(char * argv[], int nproc, int iproc) {

	char FileName[400];//name of processed Wig file
	char Type[40];
	char T[40];
	
	strcpy(FileName, argv[1]);

	SPAM(("=================\n encoder_by_blocks.cpp start!\n"));
	
	char tmpFileName[400];
	//parameters that can be changed for testing
	//the RANGE_SIZE_IN_BITS depends on your alphabet size
	int alphabet_size;//256; //Please notice the size of the alphabet. 
	int data_size;
	//printf("data_size %d\n",data_size);
	int RANGE_SIZE_IN_BITS=30;
	//ranges_in_bits<32, it cannot be more than sizeof(int), which is 4 bytes=32 bits, otherwise makeRanges does not work
	//range_in_bits is the number of bits used to represent the probabilities
	// L=64-range_in_bits is the number of bits used to store the current range of low and high
	//Need for any probability p of a symbol, p>= 1/2^(L-3). So L-3>=range_in_bits, or range_in_bits<=30
	int *data;
	double entropy;	
	int expected_size;
	//int cm_size;
	int *cm;
	//int actual_size;
	int i=0;
	int j=0;
	int k=0;
	int tmp,tmpint;
	char str[2];
	int pcount=0; //count for tmp print output
	FILE *fpDiffSeqMatlab;
	int SetCodeSize = (int)(1)<<(BlockSize+4); //code size in bytes allocated for a block
	int BMASK = ((int)(1)<<BlockSize) - 1; //mask used to do modulo 2^BlockSize
	//int SMASK = ((int)(1)<<	SearchSize) - 1; //mask used to do modulo 2^SearchSize
	int BlockAdd = 0; //address of a encode block
	RangeMapper* rm_encode;

	//long int underflow; //count number of times low and high are too close
	
	//input blockaux to get chrmadd
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"BlockAux");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpBlockAux = fopen(tmpFileName, "r");
	if (fpBlockAux  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	int BlockAux[32];
	i=0;
	while (i<32 && !feof(fpBlockAux)){
		fscanf(fpBlockAux,"%d ",&BlockAux[i]);
		i++;
	}
	if (i==32){
		printf("too many chrms!\n");
		exit(1);
	}
	fclose(fpBlockAux);
	//remove(tmpFileName);
	

	
	
	
	//write to summary result file
	FILE *fpResult=fopen("result","a");
	if (fpResult == NULL) {
		fprintf(stderr,"Can't open output result file!\n");
		exit(1);
	}

	//start encoding
	for (int type=0; type<2; type++){
		clock_t start_encoding = clock();
	
		if (type==0){
			strcpy(Type,"DiffSeqMatlab");
			strcpy(T,"Diff");
			//alphabet_size = 6134;
		}
		else{
			strcpy(Type,"LenDiffSeqMatlab");
			strcpy(T,"LenDiff");
			//alphabet_size = 67955;
		}		
		
		//load data from file
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,Type);
		if (nproc>0){
			tmpint = sprintf(str,"%02d",iproc);
			strcat(tmpFileName,str);
		}
		fpDiffSeqMatlab = fopen(tmpFileName, "r");
		if (fpDiffSeqMatlab  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}
		
		///////////make ranges////////////
		alphabet_size = 0;
		data_size=100000000;
		data = (int*)malloc(data_size * sizeof(int));
		i=0;
		while (!feof(fpDiffSeqMatlab) && i<data_size){
			fscanf(fpDiffSeqMatlab,"%d ", &data[i]);
			/*if (data[i]<0){
				printf("Wrong range! data[%d]=%d\n",i,data[i]);
				exit(1);
			}
			if (data[i]+1>alphabet_size)
				alphabet_size = data[i]+1;*/
			i++;
		}
		if (i==data_size){
			printf("Too large data set! Need to increase data_size!\n");
			exit(1);
		}		
		data_size = i;
		
		
		
		entropy = calculate_entropy(data,data_size);
		expected_size = (int)(entropy * (double)(data_size) / 8.0);
		

		//printf("\n--------\ndata[0]=%d,data[1]=%d,data[data_size-1]=%d\n",data[0],data[1],data[data_size-1]);

		//make ranges for data
		//printf("start makeranges\n");
		//load count table file
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,T);
		FILE *fpDiff = fopen(tmpFileName, "r");
		if (fpDiff  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}
		
		absize(&alphabet_size,fpDiff);
		cm = (int*)malloc((alphabet_size+2)*sizeof(int));
		makeRanges(cm, alphabet_size, RANGE_SIZE_IN_BITS,fpDiff);
		fclose(fpDiff);
		
		
		//encode
		///////////end make ranges////////////
		
		///////////encode by blocks////////////
		if (data_size & BMASK){ //not a multiple of block size
			printf("File size=%d is a not a multiple of SetBlockSize=%f!\n",data_size,exp2(BlockSize));
			exit(1);
		}

		//write to a file
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,Type);
		strcat(tmpFileName,"Code");
		if (nproc>0){
			tmpint = sprintf(str,"%02d",iproc);
			strcat(tmpFileName,str);
		}
		FILE *fpCode = fopen(tmpFileName, "wb");
		if (fpCode  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}
		//write BlockAdd to a file
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,Type);
		strcat(tmpFileName,"BlockAdd");
		if (nproc>0){
			tmpint = sprintf(str,"%02d",iproc);
			strcat(tmpFileName,str);
		}
		FILE *fpBlockAdd = fopen(tmpFileName, "wb");
		if (fpBlockAdd  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}
		
		//output chrmadd storing the starting addr of each chrm
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,Type);
		strcat(tmpFileName,"ChrmAdd");
		if (nproc>0){
			tmpint = sprintf(str,"%02d",iproc);
			strcat(tmpFileName,str);
		}
		FILE *fpChrmAdd = fopen(tmpFileName, "w");
		if (fpChrmAdd  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}
		//first chrmadd=0
		fprintf(fpChrmAdd,"0 ");
		
		g_buffer = (unsigned char*)malloc(SetCodeSize);
		BlockAdd = 0; //address in bytes of a block in the code
		int ChrmCount = 0; //counter for chrm
		j = 0; //counter for block
//		int scount=0; //counter for SearchAdd
		current_byte = 0; 
		TMPWRITE = 0x00;
		TMPWRITE_COUNT = 7;	
		rm_encode = new RangeMapper(RANGE_SIZE_IN_BITS);
		
		//printf("last data[data_size-1] is %d, alphabet size is %d\n", data[data_size-1],alphabet_size);
		//for (i=0; i<40959; i++){
		for (i=0; i<data_size; i++){
			/*if (i>data_size-20 || !(i & SMASK)){
				printf("data[%d]= %d, BlockAdd %d\n",i,data[i],BlockAdd);
			}*/
			rm_encode->encodeRange(cm[data[i]], cm[data[i]+1]);//encode
/*			if (i==4096 || i==4095){
				printf("i=4096, LOW is %llX, HIGH is %llX, midpoint is %d\n",rm_encode->getLOW(),rm_encode->getHIGH(),rm_encode->getMidPoint());
				printf("tmpwrite_count=%hhX, TMPWRITE=%hhX \n", TMPWRITE_COUNT,TMPWRITE);
				printf("BlockAdd is %d\n",BlockAdd);
			}
*/			
			if (!((i+1) & BMASK)){ //a block is done
				/*if (pcount<10){
					printf("data_size %d, i %d, BlockAdd %d\n",data_size,i,BlockAdd);
					pcount++;
				}*/
				//old
				rm_encode->encodeRange(cm[alphabet_size], cm[alphabet_size+1]); //end of data marker
				rm_encode->flush();
				if (current_byte>SetCodeSize){
					printf("Code longer than expected! Increase SetCodeSize!\n");
					exit(1);
				}
				delete rm_encode;
				fwrite(g_buffer,1,current_byte,fpCode); //binary file write			
				//record BlockAdd
				fwrite(&BlockAdd, sizeof(int), 1, fpBlockAdd);
				//new
				rm_encode = new RangeMapper(RANGE_SIZE_IN_BITS);	
				BlockAdd += current_byte; //update next block address		
				current_byte = 0; 
				TMPWRITE = 0x00;
				TMPWRITE_COUNT = 7;
				j++;
				if (j==BlockAux[ChrmCount+1]){ //end of a chrm
					fprintf(fpChrmAdd, "%d ",BlockAdd);
					//printf("next chrm BlockAdd is %d",BlockAdd);
					ChrmCount++;
				}
				
			}
		}
		fclose(fpCode);
		free(data);
		//actual_size = BlockAdd;
		//record total size
		fwrite(&BlockAdd, sizeof(int), 1, fpBlockAdd);
		fclose(fpDiffSeqMatlab);
		fclose(fpBlockAdd);
		fclose(fpChrmAdd);
		free(g_buffer);
		free(cm);
		
		//end encoding
		clock_t end_encoding = clock();
		SPAMR((fpResult,"\n\n//combine/encoder_by_blocks\n//File is %s%s\n", FileName, Type));
		SPAMR((fpResult,"Time for encoding %2.3f sec.\n", (double)(end_encoding - start_encoding)/CLOCKS_PER_SEC));
		SPAMR((fpResult,"Diff Seq length is %d, entropy is %2.3f, non-zero-count alphabet_size is %d, RANGE_SIZE_IN_BITS is %d.\n",
		data_size, entropy, alphabet_size,RANGE_SIZE_IN_BITS));
		SPAMR((fpResult,"Expected size %d bytes, Actual size =BlockAdd=  %d bytes.\n", expected_size,BlockAdd));
		
		/*printf("\n\n//combine/encoder_by_blocks\n//File is %s%s\n", FileName, Type);
		printf("Time for encoding %2.3f sec.\n", (double)(end_encoding - start_encoding)/CLOCKS_PER_SEC);
		printf("Diff Seq length is %d, entropy is %2.3f, non-zero-count alphabet_size is %d, RANGE_SIZE_IN_BITS is %d.\n",
		data_size, entropy, alphabet_size,RANGE_SIZE_IN_BITS);
		printf("Expected size %d bytes, Actual size =BlockAdd=  %d bytes. Number of SearchAdd is %d\n", expected_size,BlockAdd,scount);
		*/
	}
	
	

	
	
	
	fclose(fpResult);

	return 0;
	
}
