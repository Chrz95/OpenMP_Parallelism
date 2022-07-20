#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#define MAXCHAR 100

// Global Variables
double time_E = 0 ;
double time_F = 0;
long long int NumOfCells = 0 ;
int MaxValue ;
int Traceback_Steps = 0 ;
int cnt = 0 ;
pthread_mutex_t MutexMaxValue , mutexcnt;
int NumOfThreads ;
int * ScoringMatrix;
long int Start_Diag_I , Start_Diag_J ;
pthread_barrier_t func_barrier;
char* FirstString  ;
char* SecondString  ;
int length1  ;
int length2  ;
int gap  ;
int match  ;
int mismatch  ;
int D_Sz_all ; 
int Q_max ; 
int block_size = 1  ; 

struct CSM_Arguments
{
	int j  ;
	int NumOfDiagElements;
	int Start_Diag_I;
	int Start_Diag_J;
	int* StartPointI ;
	int* StartPointJ ;
	int MaxValue ;
	int Thread_ID ;
};

double gettime(void)
{
	struct timeval ttime;

	gettimeofday(&ttime, NULL);
	return ttime.tv_sec+ttime.tv_usec * 0.000001;
}

char* reverseConstString(char const* str)
{
    int n = strlen(str);
		char temp ;
    char *rev = (char*) malloc(n+1);

    strcpy(rev, str);

    for (int i=0, j=n-1; i<j; i++,j--)
		{
			temp = rev[i] ;
			rev[i] = rev[j];
			rev[j] = temp ;
		}

    // return pointer of reversed string
    return rev;
}

int totalMax (int A,int B,int C,int D)
{
	int max = 0 ;

	if (A > B)
		max = A ;
	else
		max = B ;

	if (C > max)
		max = C ;

	if (D > max)
		max = D ;

	return max ;
}

void *CalcScoringMatrix (void *Args)
{
	struct CSM_Arguments *my_data ;
	int Current_Diag_I ;
	int Current_Diag_J ;
	my_data = (struct CSM_Arguments *) Args ;

	// Passing the arguments
	int j = my_data->j ;
	int NumOfDiagElements = my_data->NumOfDiagElements;
	int* StartPointI = my_data->StartPointI;
	int* StartPointJ = my_data->StartPointJ;

	for (int i = 0 ; i < block_size ; i++)
	{
		Current_Diag_I = Start_Diag_I + j - 1 + i;
		Current_Diag_J = Start_Diag_J - j + 1 - i;

		if ((Current_Diag_I < 1) || (Current_Diag_I > length2) || (Current_Diag_J < 1) || (Current_Diag_J > length1)) break ;
		else
		{
			if (SecondString[Current_Diag_I-1] == FirstString[Current_Diag_J-1]) ScoringMatrix[Current_Diag_J + Current_Diag_I * (length1 + 1)] = totalMax(0,ScoringMatrix[Current_Diag_J + (Current_Diag_I-1) * (length1 + 1)] + gap,ScoringMatrix[(Current_Diag_J-1) + Current_Diag_I * (length1 + 1)] + gap,ScoringMatrix[(Current_Diag_J -1) + (Current_Diag_I-1) * (length1 + 1)] + match);
			else ScoringMatrix[Current_Diag_J + Current_Diag_I * (length1 + 1)] = totalMax(0,ScoringMatrix[Current_Diag_J + (Current_Diag_I-1) * (length1 + 1)] + gap,ScoringMatrix[(Current_Diag_J-1) + Current_Diag_I * (length1 + 1)] + gap,ScoringMatrix[(Current_Diag_J-1) + (Current_Diag_I-1) * (length1 + 1)] + mismatch);

			if ((ScoringMatrix[Current_Diag_J + Current_Diag_I * (length1 + 1)] > MaxValue) && (ScoringMatrix[Current_Diag_J + Current_Diag_I * (length1 + 1)] !=0))
			{
					pthread_mutex_lock (&mutexcnt);
					cnt = 0;
					pthread_mutex_unlock (&mutexcnt);

		     	    for (int k = 0 ; k < length2 + 1 ; k++) // Re Initialize the Start Point matrix
					{
							 StartPointI[k] = 0 ;
							 StartPointJ[k] = 0 ;
					}

					pthread_mutex_lock (&MutexMaxValue);
					MaxValue = ScoringMatrix [Current_Diag_J + Current_Diag_I * (length1 + 1)] ;
					//	 printf ("MaxValue is %d \n",MaxValue);
					pthread_mutex_unlock (&MutexMaxValue);
					StartPointI[cnt] = Current_Diag_I ;
					StartPointJ[cnt] = Current_Diag_J ;
			}
			else if (ScoringMatrix[Current_Diag_J + Current_Diag_I * (length1 + 1)] == MaxValue)
			{
					pthread_mutex_lock (&mutexcnt);
					cnt++;
					pthread_mutex_unlock (&mutexcnt);

					if (cnt < length2 + 1)
					{
							StartPointI[cnt] = Current_Diag_I ;
							StartPointJ[cnt] = Current_Diag_J ;
					}
			}		
		}		
  }
}

struct TB_Arguments
{
	int* StartPointI ;
	int* StartPointJ ;
	int i;
	char* Q_traceback;
	char* D_traceback;
	FILE *fp1 ;
};

void *Traceback (void *Args)
{
		struct TB_Arguments *my_data ;
		my_data = (struct TB_Arguments *) Args ;
		int i = my_data->i ;
		int* StartPointI = my_data->StartPointI;
		int* StartPointJ = my_data->StartPointJ;
		FILE * fp1 = my_data->fp1 ;

		int trace_cnt = 0 ;
		int traceback_length  = length2/100 + 2 ;
		int I = StartPointI[i]  ;
		int J = StartPointJ[i] ;
		int index = J + I * (length1 + 1) ;
		int index_hor = (J-1) + I * (length1 + 1) ;
		int index_vert =  J + (I-1) * (length1 + 1) ;
		int index_diag = (J-1) + (I-1) * (length1 + 1) ;

		char* Q_traceback = (char *) malloc(traceback_length*length1) ;
		char* D_traceback = (char *) malloc(traceback_length*length1) ;

		for (int j=0 ; j < traceback_length*length1 ; j++) // Initialize traceback matrices
		{
			Q_traceback[j] = '\0';
			D_traceback[j] = '\0';
		}
		while ((ScoringMatrix[index_diag] != 0) || (ScoringMatrix[index_hor] != 0) || (ScoringMatrix[index_vert] != 0)) // All surrounding elements must be zero to stop
		{
			index = J + I * (length1 + 1) ;
			index_hor = (J-1) + I * (length1 + 1) ;
			index_vert =  J + (I-1) * (length1 + 1) ;
			index_diag = (J-1) + (I-1) * (length1 + 1) ;

			if ((ScoringMatrix[index_diag] >= ScoringMatrix[index_hor]) && (ScoringMatrix[index_diag] >= ScoringMatrix[index_vert]))
			{
				Q_traceback[trace_cnt] = FirstString[J - 1] ;
				D_traceback[trace_cnt] = SecondString[I -  1] ;
				trace_cnt++;

				Traceback_Steps ++ ;
				I-- ;
				J-- ;
			}
			else if (ScoringMatrix[index_hor] >= ScoringMatrix[index_vert]) // Horizontal
			{
				Q_traceback[trace_cnt] = FirstString[J - 1] ;
				D_traceback[trace_cnt] = '-' ;
				trace_cnt++;
				Traceback_Steps++ ;
				J-- ;
			}
			else // if (ScoringMatrix[I -1][J] >= ScoringMatrix[I][J -1]) // Vertical
			{
				Q_traceback[trace_cnt] = '-' ;
				D_traceback[trace_cnt] =  SecondString[I - 1] ;
				trace_cnt++;
				Traceback_Steps++ ;
				I-- ;
			}
		}

		int Stop = StartPointI[i] - 1 ;
		int Start = I ;

		fprintf (fp1,"\nMatch %d [Score: %d,Start: %d,Stop: %d]\n",i+1,MaxValue,Start,Stop);
		fprintf (fp1,"D: %s \n",reverseConstString(D_traceback));
		fprintf (fp1,"Q: %s \n",reverseConstString(Q_traceback));

		free (Q_traceback);
		free (D_traceback);
	}

 int Smith_Waterman (char* output_filename,int match,int mismatch,int gap,FILE **fp)
 {
	 	length1 = Q_max; // Q
	 	length2 = D_Sz_all; // D
	 	int StartPointI[length2 + 1];
		int StartPointJ[length1 + 1];
		int index = 0;
		int index_hor = 0 ;
		int index_vert = 0 ;
		int index_diag = 0 ;
		FILE* fp1 = *fp ;
		double timeE_0 = gettime();
		double timeE_1 = 0 ;
		double timeF_0 , timeF_1 ;
		pthread_t threads[NumOfThreads];
		struct CSM_Arguments CSM_Arguments_Array[NumOfThreads]; // Thread arguments
		int rc ;

		// Scoring Matrix Calculation
		for (int i = 0 ; i <  length2 + 1 ; i++) // Initialize first column to 0s
		{
			index = i * (length1 + 1) ;
			ScoringMatrix[index] = 0 ;
		}

		for (int j = 0 ; j <  length1 + 1 ; j++) // Initialize first row to 0s
		{
			index = j ;
			ScoringMatrix[index] = 0 ;
		}

		int i, j ;
		long int NumOfRevDiags = length1 + length2 - 1 ;
		long int NumOfDiagElements ;

        for (i = 1; i <= NumOfRevDiags; i++)
        {
			// Find out how many elements the current diagonal has
			if ((i < length1 + 1) && ((i < length2 + 1))) NumOfDiagElements =  i; // Diagonal is increasing in size
			else if (i < length2 + 1) NumOfDiagElements = length1 ; // Diagonal remains the same size
			else NumOfDiagElements = 2 * (length1 + 1) - i + (length2 - length1) - 2; // Diagonal is decreasing in size

			// Find the position of the first element of the current diagonal
			if (i <= length1) // Start Point is in same line
			{
					Start_Diag_I = 1;
					Start_Diag_J = i;
			}
			else // Start Point is in same Column
			{
					Start_Diag_I = i - length1 + 1;
					Start_Diag_J = length1 ;
			}

			//	printf ("(Start_Diag_I,Start_Diag_J) = (%d,%d)", Start_Diag_I,Start_Diag_J) ;
			block_size = NumOfDiagElements / NumOfThreads ; 
			if (block_size == 0) block_size = 1 ; 
	
            for (j = 1 ; j <= NumOfDiagElements; j = j + block_size)
            {
					CSM_Arguments_Array[(j-1) % NumOfThreads].j = j ;
					CSM_Arguments_Array[(j-1) % NumOfThreads].NumOfDiagElements = NumOfDiagElements ;
					CSM_Arguments_Array[(j-1) % NumOfThreads].StartPointI = StartPointI ;
					CSM_Arguments_Array[(j-1) % NumOfThreads].StartPointJ = StartPointJ  ;
					int rc = pthread_create(&threads[(j-1) % NumOfThreads], NULL, &CalcScoringMatrix, (void *) &CSM_Arguments_Array[(j-1) % NumOfThreads]) ;
    			    if (rc)
					{
							printf("Problem \n");
							exit(-1);
					}
					if (pthread_join(threads[(j-1) % NumOfThreads], NULL) )
					{
			 				printf("error joining thread.");
			 				abort();
	 				}

            }	
	}
				
	cnt = 0 ;

	while (StartPointI[cnt] != 0) // Calculate number of Traceback Starting Points
	{
		//printf ("Staring Points are (%d,%d) \n",StartPointI[cnt],StartPointJ[cnt]);
		cnt ++ ;
	}

	timeE_1 = gettime();
	time_E += timeE_1 - timeE_0 ;

 	// Traceback
	timeF_0 = gettime();

	int ThreadsForTraceBack = NumOfThreads;
	if (cnt < NumOfThreads) ThreadsForTraceBack = cnt ;

	pthread_t threads1[NumOfThreads] ;

	struct TB_Arguments TB_Arguments_Array[ThreadsForTraceBack];

	for (int i=0 ; i < cnt ; i++)
	{
		TB_Arguments_Array[i % NumOfThreads].i = i ;
		TB_Arguments_Array[i % NumOfThreads].StartPointI = StartPointI ;
		TB_Arguments_Array[i % NumOfThreads].StartPointJ = StartPointJ  ;
		TB_Arguments_Array[i % NumOfThreads].fp1 = fp1  ;

		int rc = pthread_create(&threads1[i % NumOfThreads], NULL, &Traceback, (void *) &TB_Arguments_Array[i % NumOfThreads]) ;
	  	if (rc)
		{
				printf("Problem \n");
				exit(-1);
		}

		if (pthread_join(threads1[i % NumOfThreads], NULL) )
		{
	 			printf("error joining thread.");
	 			abort();
 		}
	}

	timeF_1 = gettime();
	time_F += timeF_1 - timeF_0 ;

	return Traceback_Steps ;
 }

 void Read_Dataset (char* output_filename,char* path,int match,int mismatch,int gap)
 {
 	FILE *fp;
    char* str= (char *) malloc(MAXCHAR);
    int pairs , Q_min  ;
    int count = 0 ;
    char* delim = (char *) malloc(2) ;
    char* token = (char *) malloc(MAXCHAR);
    int k = 0 ;
	char* Q_string ;
	char* D_string ;
    char ch ;
    char* ch_str = (char *) malloc (2) ;
	int Traceback_Steps = 0 ;

    strcpy(delim,":");

	fp = fopen(path,"r"); // Input file

    if (fp == NULL)
    {
        printf("Could not open file %s",path);
        return ;
    }

	FILE* fp1 = fopen(output_filename,"w"); // Output file

	if (fp1 == NULL)
	{
			printf("Could not open file %s",output_filename);
			return ;
	}

    while ((fgets(str, MAXCHAR, fp) != NULL) && (count <=3)) // For starters , read only the variable lines
    {
	   	// get the first token
   		token = strtok(str,delim);

  		// walk through other tokens
	   	if (token != NULL)  // Extract the number from each line
	   	{
	    	token = strtok(NULL,delim);
	    	if (count == 0) pairs = atoi (token) ;
				if (count == 1) Q_min = atoi (token) ;
				if (count == 2) Q_max = atoi (token) ;
				if (count == 3) D_Sz_all = atoi (token) ;
	   	}
	   	count++;
	}

	ScoringMatrix = (int *) malloc((D_Sz_all + 1) * (Q_max + 1) * sizeof(int));

	printf ("A) Total Number of Pairs : %d \n",pairs);
	free(delim);
	free(str);

	for (int i=0 ; i< pairs ; i++)
	{
		///////////////// Read a pair of strings /////////////////////////////////////
		Q_string = (char *) malloc(Q_max + 2); // Initialize Q_string
		D_string = (char *) malloc(D_Sz_all + 2); // Initialize D_string
		ch = fgetc(fp) ;

		while((ch != EOF) && (ch != 'D'))
		{
			ch = fgetc(fp) ;
			if ((ch != ':') && (ch != ' ') && (ch != '\n') && (ch != '\t')  && (ch  != EOF) && (ch  != 'Q') && (ch  != 'D'))
			{
				ch_str[0] = ch ;
				strcat(Q_string,ch_str);
			}
		}

		fseek (fp,-1,SEEK_CUR); // Move the pointer left of D

		while((ch != EOF) && (ch != 'Q'))
		{
			ch = fgetc(fp) ;
			if ((ch != ':') && (ch != ' ') && (ch != '\n') && (ch != '\t')  && (ch  != EOF) && (ch  != 'D') && (ch  != 'Q'))
			{
				ch_str[0] = ch ;
				strcat(D_string,ch_str);
			}
		}

		fseek (fp,-1,SEEK_CUR);

		fprintf(fp1,"\nPair %d : \n",i +1);

		NumOfCells += (strlen(Q_string))*(strlen(D_string));

		// Execute the Smith-Waterman algorith for each pair of strings
		FirstString = Q_string ;
		SecondString = D_string ;
		Traceback_Steps += Smith_Waterman (output_filename,match,mismatch,gap,&fp1);

		free (Q_string);
		free (D_string);
	}

	printf ("B) Total Number of Cells that take value : %lld \n",NumOfCells);
	printf ("C) Total Number of Traceback Steps : %d \n",Traceback_Steps);

	free (ch_str);
	fclose(fp1);
	fclose(fp);
	return;
}

void main(int argc, char *argv[])
{
	double time0 = gettime();
	double time1 = 0 ;
	double elapsed_time = 0 ;
	match = atoi(argv[3]);
    mismatch = atoi(argv[4]);
    gap = atoi(argv[5]);
	NumOfThreads = atoi(argv[6]);
    char* output_filename = (char *) malloc(MAXCHAR);
    char* path = (char *) malloc(MAXCHAR);
    strcpy(output_filename,"Report_");
    output_filename = strcat(output_filename,argv[1]);
	output_filename = strcat(output_filename,"_PTH_");
	output_filename = strcat(output_filename,argv[6]);
	output_filename = strcat(output_filename,".txt");
    path = argv[2];

	// Read the Dataset
    Read_Dataset(output_filename,path,match,mismatch,gap);

	// Print the performance indicators
	time1  = gettime();
	elapsed_time = time1 - time0 ;
	printf("D) Program Execution Time : %lf sec\n",elapsed_time);
	printf("E) Scoring Matrix Calculation Time : %lf sec\n",time_E);
	printf("F) Traceback Time : %lf sec\n",time_F);
	printf("G) CUPS (Program Execution Time) : %lf CUPS\n",NumOfCells/elapsed_time);
	printf("H) CUPS (Cell Calculation Time) : %lf CUPS\n",NumOfCells/time_E);
	free(ScoringMatrix) ;
	pthread_exit(NULL);
}
