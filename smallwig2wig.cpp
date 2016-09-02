/*
 * decomopress to a wig file
 * created by Zhiying Wang
 */


#include "parameters.h"
#include "arithmetic.h"
#include "decoder_by_blocks.cpp"
#include "decoder.cpp"
#include "partial_decoder_by_blocks.cpp"
#include "to_wig.c"
#include "to_wig_by_blocks.c"
#include "to_wig_by_blocks_partial.c"
#include "to_wig_sub.c"
#include "fsplit.c"
#include "fmergesub.c"
#include "fremove.c"
#include "lpaq1.cpp"
#include "prll.c"



int main(int argc, char * argv[]){
	
	clock_t start = clock();
	
	//int m = 0; //indicate if -m is enabled
	int s = 0; //indicate if -s is enalbed 
	int p = 0; //indicate if -p is enabled
	int opt; //indicate which functions to use
	bool isOK = true; //indicate if argv format is OK
	int i,j;
	//char N[10];
	int nproc;
	int Qstart, Qend;
	char ChrmName[8];
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
        /*if (strcmp(argv[i], "-m") == 0){  // Process optional arguments. 
            m = 1;
            if (argc!=4)
				isOK=false;
		}
		else*/ 
		if (strcmp(argv[i], "-s") == 0){
			s = 1;
			if ((argc==7) && i+3 <= argc){
				i++;
				strcpy(ChrmName, argv[i]);
				i++;
				Qstart = atoi(argv[i]);
				i++;
				Qend = atoi(argv[i]);
			}
			else
				isOK=false;
		}
		else if (strcmp(argv[i], "-p") == 0){
			p = 1;
			if ((argc==5) && i+1 <= argc){
				i++;
				nproc = atoi(argv[i]);
			}
			else
				isOK=false;
		}		
		else{
			strcpy(targv[j],argv[i]);
			j++;
		}
	}
	if (j!=3){//targv should only left two arguments: path, filename
		isOK = false;
	}
	if (isOK){
		if (s==0 && p==0){ //decode whole file
			fsplit(&opt,targv);
			switch (opt){
				case 0: //regular
				{
					decoder(targv);
					to_wig(targv);
					fremove(opt,targv,-1);
					break;
				}
				case 1: //lpaq, or context mixing
				{
					//make targv same as format for lpaq1
					strcpy(targv2[1],"d");//decode option
					char FileName[1000];
					//Diff
					strcpy(targv2[3],targv[1]);
					strcat(targv2[3],"DiffSeqMatlabDecode");//output
					strcpy(targv2[2],targv[1]);
					strcat(targv2[2],"DiffSeqMatlablpaq");//input
					lpaq1(4,targv2);
					//LenDiff
					strcpy(targv2[3],targv[1]);
					strcat(targv2[3],"LenDiffSeqMatlabDecode");//output
					strcpy(targv2[2],targv[1]);
					strcat(targv2[2],"LenDiffSeqMatlablpaq");//input
					lpaq1(4,targv2);
					to_wig(targv);
					fremove(opt,targv,-1);
					break;
				}
				case 2: //encode by blocks
				{
					decoder_by_blocks(targv,-1,-1);
					to_wig_by_blocks(targv,-1,-1);
					/////////////////???????????????????????????fremove(opt,targv,-1);
					break;
				}
				default:
					printf("Input file is wrong! Please check again!\n");
					exit(1);
			}
			
		}
		else if ( s==1 && p==0){ //subsequence access
			fsplit(&opt,targv);
			switch (opt){
				case 0: //regular encoded
				{
					decoder(targv);
					to_wig_sub(targv,ChrmName,Qstart,Qend);
					fremove(opt,targv,-1);
					break;
				}
				case 1: //lpaq encoded
				{
					//make targv same as format for lpaq1
					strcpy(targv2[1],"d");//decode option
					char FileName[1000];
					//Diff
					strcpy(targv2[3],targv[1]);
					strcat(targv2[3],"DiffSeqMatlabDecode");//output
					strcpy(targv2[2],targv[1]);
					strcat(targv2[2],"DiffSeqMatlablpaq");//input
					lpaq1(4,targv2);
					//LenDiff
					strcpy(targv2[3],targv[1]);
					strcat(targv2[3],"LenDiffSeqMatlabDecode");//output
					strcpy(targv2[2],targv[1]);
					strcat(targv2[2],"LenDiffSeqMatlablpaq");//input
					lpaq1(4,targv2);
					to_wig_sub(targv,ChrmName,Qstart,Qend);
					fremove(opt,targv,-1);
					break;
				}
				case 2: //encoded by blocks
				{
					int Bstart;
					int Bend;
					partial_decoder_by_blocks(targv,ChrmName,Qstart,Qend, &Bstart, &Bend);
					to_wig_by_blocks_partial(targv,ChrmName,Qstart,Qend, &Bstart, &Bend);
					fremove(opt,targv,-1);
					break;
				}		
				default:
					printf("Input file is wrong! Please check again!\n");
					exit(1);
			}
		}
		else if ( s==0 && p==1){ //parallel
			fsplit(&opt,targv);
			if (opt==2){
				prll(&nproc, targv,decoder_by_blocks);
				prll(&nproc, targv,to_wig_by_blocks);
				fremove(3,targv,nproc);
			}	
			else{
				printf("Parallel option is not available for this file! opt=%d\n",opt);
				fremove(opt,targv,-1);
				isOK = false;
			}		
		}
		else
			isOK = false;
	}
	if (isOK==false)
		printf("To use:\n smallwig2wig [InputFile] [OutputFile] \n"
			"options: \n"
			"-s [ChrmName (e.g. chr1)] [Query Start (integer)] [Query End (integer)] \n"
			"	subsequence query \n"
			"-p [total number of processes]\n"
			"	parallel realization, only availabe if -s is disabled and -r is enabled in encoding\n");
    else{
		clock_t end = clock();
		SPAM(("Decompressing time in seconds\n%2.3f\n", (double)(end-start)/CLOCKS_PER_SEC));
		if (s==1)
			SPAM(("query time per bp %e\n", ((double)(end-start)/CLOCKS_PER_SEC)/(Qend-Qstart+1)));
	}
    for(i = 0; i < 10; i++){
		free(targv[i]);
		free(targv2[i]);
	}
	free(targv);
	free(targv2); 
	
	return 0;

}

