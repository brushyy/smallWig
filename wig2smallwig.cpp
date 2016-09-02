/*
 * comopress a wig file 
 * created by Zhiying Wang
 */


#include "parameters.h"
#include "arithmetic.h"
#include "encoder_by_blocks.cpp"
#include "encoder.cpp"
#include "getseqblock.c"
#include "getmatlabseq.c"
#include "getseq.c"
#include "fmergeadd.c"
#include "fmergesub.c"
#include "fmerge.c"
#include "prll.c"
#include "fmergecount.c"
#include "lpaq1.cpp"


int main(int argc, char * argv[]){
	
	clock_t start = clock();
	
	//printf("compress start!\n");
	int m = 0; //indicate if -m is enabled
	int r = 0; //indicate if -r is enalbed 
	int p = 0; //indicate if -p is enabled
	//int f = 0; //indicate if -f is enabled
	//int opt = 0; //indicate which functions to use
	bool isOK = true; //indicate if argv format is OK
	int i,j;
	char N[10];
	int nproc;
	char **targv; //tmp argv to put to subsequent functions, if format correct, only stores two strings: path, filename
	char **targv2;
	targv = (char **)malloc(10 * sizeof(char *));
	targv2 = (char **)malloc(10 * sizeof(char *));
	for(i = 0; i < 10; i++){
		targv[i] = (char *)malloc(10000 * sizeof(char));
		targv2[i] = (char *)malloc(10000 * sizeof(char));
	}
	j=1;
	for (i = 1; i<10 && i < argc && isOK; i++){ //should argc<10
        if (strcmp(argv[i], "-m") == 0){  /* Process optional arguments. */
            m = 1;
            if (argc==5 && i+1 <= argc){
				i++;
				strcpy(N, argv[i]);
			}
			else
				isOK=false;
		}
		else if (strcmp(argv[i], "-r") == 0){
			r = 1;
			if ((argc==5 || argc==7 ) && i+1 <= argc){
				i++;
				BlockSize = atoi(argv[i]);
				//printf("blocksize %d\n",BlockSize);
			}
			else
				isOK=false;
		}
		else if (strcmp(argv[i], "-p") == 0){
			p = 1;
			if ((argc==7 ) && i+1 <= argc){
				i++;
				nproc = atoi(argv[i]);
			}
			else
				isOK=false;
		}
		/*else if (strcmp(argv[i],"-f") == 0){
			f = 1;
			if ((argc==7 || argc==9) && i+1 <= argc){
				i++;
				SetFloatSize = pow10(atoi(argv[i]));
			}
			else
				isOK=false;
		}*/
		else{
			strcpy(targv[j],argv[i]);
			j++;
		}
	}
	if (isOK){
		
		if (m==0 && r==0 && p==0 && argc==3){ //regular
			getseq(targv);
			encoder(targv);
			fmerge(0,targv,-1);
		}
		else if (m==1 && r==0 && p==0 ){ //context mixing
			getseq(targv);
	
			//make targv same as format for lpaq1
			strcpy(targv2[1],N);//N
			char FileName[1000];
			//Diff
			strcpy(targv2[2],targv[1]);
			strcat(targv2[2],"DiffSeqMatlab");//input
			strcpy(targv2[3],targv[1]);
			strcat(targv2[3],"DiffSeqMatlablpaq");//output
			lpaq1(4,targv2);
			//LenDiff
			strcpy(targv2[2],targv[1]);
			strcat(targv2[2],"LenDiffSeqMatlab");//input
			strcpy(targv2[3],targv[1]);
			strcat(targv2[3],"LenDiffSeqMatlablpaq");//output
			lpaq1(4,targv2);
			
			fmerge(1,targv,-1);
		}
		else if (m==0 && r==1 && p==0){ //random access
			//printf("random access!\n");
			getseqblock(targv,-1,-1);
			getmatlabseq(targv,-1,-1);
			encoder_by_blocks(targv,-1,-1);
			fmerge(2,targv,-1);
		}
		else if (m==0 && r==1 && p==1){ //random access & parallel
			prll(&nproc, targv,getseqblock);
			fmergecount(nproc,targv); //merge count files from different processes to get overall frequency
			prll(&nproc, targv,getmatlabseq);
			prll(&nproc, targv,encoder_by_blocks);
			fmerge(3,targv,nproc);
		}
		else
			isOK = false;
	}
	if (isOK==false)
		printf("To use:\nwig2smallwig InputFile OutputFile \n"
			"options: \n"
			"-m [N=0..9, uses 3 + 3*2^N MB memory, decompress should use same N] \n"
			"	context mixing \n"
			"-r [B, encode block size from 8 to 32]\n"
			"	random access and encode by blocks of size 2^B\n"
			"-p [total number of processes] \n"
			"	parallel realization, only available if -r is enabled and -m is disabled\n");
			//"-f [F=0,1,..., number of digits after decimal point for floating-point values]\n"
			//"   floating point values, available if -r is enabled and -m is disabled\n");
    else{
		clock_t end = clock();
		SPAM(("Compressing time in seconds\n%2.3f\n", (double)(end-start)/CLOCKS_PER_SEC));
	}
    for(i = 0; i < 10; i++){
		free(targv[i]); 
		free(targv2[i]);
	}
	free(targv);
	free(targv2);
	
	
           
	return 0;

}
