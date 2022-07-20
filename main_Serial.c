#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#define MAXCHAR 100

// Global Variables
double time_E = 0 ;
double time_F = 0;
long long int NumOfCells = 0 ;

double gettime(void)
{
	struct timeval ttime;
	gettimeofday(&ttime, NULL);
	return ttime.tv_sec+ttime.tv_usec * 0.000001;
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

 int Smith_Waterman (char* FirstString,char* SecondString,char* output_filename,int match,int mismatch,int gap,FILE **fp)
 {
	 	int length1 = strlen(FirstString); // Q
	 	int length2 = strlen(SecondString); // D
	 	int StartPointI[length2 + 1];
		int StartPointJ[length1 + 1];
	 	int MaxValue = 0 , cnt = 0 ;
		int *ScoringMatrix = malloc((length2 + 1) * (length1 + 1) * sizeof(int));
		int index = 0;
		int index_hor = 0 ;
		int index_vert = 0 ;
		int index_diag = 0 ;
		FILE* fp1 = *fp ;
		double timeE_0 = gettime();
		double timeE_1 = 0 ;
		double timeF_0 , timeF_1 ;

		// printf ("The Scoring Matrix has a size of %d bytes \n",(length2 + 1) * (length1 + 1) * sizeof(int)) ;

		// Scoring Matrix Calculation

		for (int i = 0 ; i <  length2 + 1 ; i++) // Initialize first column to 0
		{
			index = i * (length1 + 1) ;
			ScoringMatrix[index] = 0 ;
		}

		for (int j = 0 ; j <  length1 + 1 ; j++) // Initialize first row to 0
		{
			index = j ;
			ScoringMatrix[index] = 0 ;
		}

		index = 0 ;

		for (int i = 1 ; i < length2 + 1  ; i++) // Calculating the marix
		{
				for (int j = 1 ; j < length1 + 1; j++)
				{
							index = j + i * (length1 + 1) ;
							index_hor = (j-1) + i * (length1 + 1) ;
							index_vert =  j + (i-1) * (length1 + 1) ;
							index_diag = (j-1) + (i-1) * (length1 + 1) ;

							if (SecondString[i-1] == FirstString[j - 1]) ScoringMatrix[index] = totalMax(0,ScoringMatrix[index_vert] + gap,ScoringMatrix[index_hor] + gap,ScoringMatrix[index_diag] + match);
							else ScoringMatrix[index] = totalMax(0,ScoringMatrix[index_vert] + gap,ScoringMatrix[index_hor] + gap,ScoringMatrix[index_diag] + mismatch);

							if ((ScoringMatrix[index] > MaxValue) && (ScoringMatrix[index] !=0)) // Each time a new max value is found 
							{
								cnt = 0;

								for (int k = 0 ; k < length2 + 1 ; k++) // Re Initialize the Start Point matrix
								{
										StartPointI[k] = 0 ;
										StartPointJ[k] = 0 ;
								}

								MaxValue = ScoringMatrix[index] ;
								StartPointI[cnt] = i ;
								StartPointJ[cnt] = j ;
							}
							else if (ScoringMatrix[index] == MaxValue) // If the same max value is found, place its coordinates in the StartPoint matricrs
							{
									cnt++ ;
									if (cnt < length2 + 1)
									{
											StartPointI[cnt] = i  ;
											StartPointJ[cnt] = j  ;
									}
							}
					}
			}

	timeE_1 = gettime();
	time_E += timeE_1 - timeE_0 ;

 	// Traceback

	timeF_0 = gettime();

	int I = 0 , J = 0 ;
	int trace_cnt = 0 ;
	int Start = 0 , Stop = 0 ;
	int traceback_length = 0 ;
	int Traceback_Steps = 0 ;

	traceback_length = length2/100 + 2 ; // Εmperical equation 

	char* Q_traceback = (char *) malloc(traceback_length*length1) ;
	char* D_traceback = (char *) malloc(traceback_length*length1) ;

	index = 0 ;
	index_hor = 0 ;
	index_vert =  0 ;
	index_diag = 0 ;

	for (int i=0 ; i <= cnt ; i++)
	{
		trace_cnt = 0 ;

		I = StartPointI[i]  ;
		J = StartPointJ[i] ;

		for (int i=0 ; i < traceback_length*length1 ; i++) // Initialize traceback matrices
		{
			Q_traceback[i] = '\0';
			D_traceback[i] = '\0';
		}

		index = J + I * (length1 + 1) ;
		index_hor = (J-1) + I * (length1 + 1) ;
		index_vert =  J + (I-1) * (length1 + 1) ;
		index_diag = (J-1) + (I-1) * (length1 + 1) ;

		while ((ScoringMatrix[index_diag] != 0) || (ScoringMatrix[index_hor] != 0) || (ScoringMatrix[index_vert] != 0)) // All surrounding elements must be zero to stop
		{
			index = J + I * (length1 + 1) ;
			index_hor = (J-1) + I * (length1 + 1) ;
			index_vert =  J + (I-1) * (length1 + 1) ;
			index_diag = (J-1) + (I-1) * (length1 + 1) ;

			if ((ScoringMatrix[index_diag] >= ScoringMatrix[index_hor]) && (ScoringMatrix[index_diag] >= ScoringMatrix[index_vert])) Diagonal
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

		Stop = StartPointI[i] - 1 ;
		Start = I ;

		// Write Match to file

		fprintf (fp1,"\nMatch %d [Score: %d,Start: %d,Stop: %d]\n",i+1,MaxValue,Start,Stop);
		fprintf (fp1,"D: %s \n",reverseConstString(D_traceback));
		fprintf (fp1,"Q: %s \n",reverseConstString(Q_traceback));
	}
		timeF_1 = gettime();
		time_F += timeF_1 - timeF_0 ;
		free(ScoringMatrix) ;
		free (Q_traceback);
		free (D_traceback);

	return Traceback_Steps;
 }

 void Read_Dataset (char* output_filename,char* path,int match,int mismatch,int gap)
 {
 		FILE *fp;
    char* str= (char *) malloc(MAXCHAR);
    int pairs , Q_min , Q_max , D_Sz_all;
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
		Traceback_Steps += Smith_Waterman (Q_string,D_string,output_filename,match,mismatch,gap,&fp1);

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

    // Assign the arguments
    char* output_filename = (char *) malloc(MAXCHAR);
    char* path = (char *) malloc(MAXCHAR);
    strcpy(output_filename,"Report_");
    output_filename = strcat(output_filename,argv[1]);
		output_filename = strcat(output_filename,".txt");
    path = argv[2];
    int match = atoi(argv[3]);
    int mismatch = atoi(argv[4]);
    int gap = atoi(argv[5]);

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
}
