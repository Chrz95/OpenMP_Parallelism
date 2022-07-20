# sudo ./run.sh -threads 4

############# - name #################
: '
if [ $1 = -name ]
then
ID=$2
fi

if [ $3 = -name ]
then
ID=$4
fi

if [ $5 = -name ]
then
ID=$6
fi

if [ $7 = -name ]
then
ID=$8
fi

if [ $9 = -name ]
then
ID=$10
fi

if [ $11 = -name ]
then
ID=$12
fi

#####################################

############# - input #################

if [ $1 = -input ]
then
PATH=$2
fi

if [ $3 = -input ]
then
PATH=$4
fi

if [ $5 = -input ]
then
PATH=$6
fi

if [ $7 = -input ]
then
PATH=$8
fi

if [ $9 = -input ]
then
PATH=$10
fi

if [ $11 = -input ]
then
PATH=$12
fi

#####################################

############# - match #################

if [ $1 = -match ]
then
INT1=$2
fi

if [ $3 = -match ]
then
INT1=$4
fi

if [ $5 = -match ]
then
INT1=$6
fi

if [ $7 = -match ]
then
INT1=$8
fi

if [ $9 = -match ]
then
INT1=$10
fi

if [ $11 = -match ]
then
INT1=$12
fi

#####################################

############# - mismatch #################

if [ $1 = -mismatch ]
then
INT2=$2
fi

if [ $3 = -mismatch ]
then
INT2=$4
fi

if [ $5 = -mismatch ]
then
INT2=$6
fi

if [ $7 = -mismatch ]
then
INT2=$8
fi

if [ $9 = -mismatch ]
then
INT2=$10
fi

if [ $11 = -mismatch ]
then
INT2=$12
fi
#####################################

############# - gap #################

if [ $1 = -gap ]
then
INT3=$2
fi

if [ $3 = -gap ]
then
INT3=$4
fi

if [ $5 = -gap ]
then
INT3=$6
fi

if [ $7 = -gap ]
then
INT3=$8
fi

if [ $9 = -gap ]
then
INT3=$10
fi

if [ $11 = -gap ]
then
INT3=$12
fi

#####################################

############# - threads #################

if [ $1 = -threads ]
then
INT4=$2
fi

if [ $3 = -threads ]
then
INT4=$4
fi

if [ $5 = -threads ]
then
INT4=$6
fi

if [ $7 = -threads ]
then
INT4=$8
fi

if [ $9 = -threads ]
then
INT4=$10
fi

if [ $11 = -threads ]
then
INT4=$12
fi
'

#####################################
echo   "\n \nExecuting Serial Program : \n"
gcc -o exec main_Serial.c

echo   "\nDataset 1 : \n"
./exec "1" "Datasets/D1.txt" 3 -1 -1
echo   "\nDataset 2 : \n"
./exec "2" "Datasets/D2.txt" 3 -1 -1
echo   "\nDataset 3 : \n"
./exec "3" "Datasets/D3.txt" 3 -1 -1
echo   "\nDataset 4 : \n"
./exec "4" "Datasets/D4.txt" 3 -1 -1
echo   "\nDataset 5 : \n"
./exec "5" "Datasets/D5.txt" 3 -1 -1
echo   "\nDataset 6 : \n"
./exec "6" "Datasets/D6.txt" 3 -1 -1
echo   "\nDataset 7 : \n"
./exec "7" "Datasets/D7.txt" 3 -1 -1
echo   "\nDataset 8 : \n"
./exec "8" "Datasets/D8.txt" 3 -1 -1
echo   "\nDataset 9 : \n"
./exec "9" "Datasets/D9.txt" 3 -1 -1
echo   "\nDataset 10 : \n"
./exec "10" "Datasets/D10.txt" 3 -1 -1


echo   "\n \nExecuting Parallel OpenMP Fine-Grained Program : \n"
gcc -fopenmp -o exec main_OMP_Fine.c

echo   "\nDataset 1 : \n"
./exec "1" "Datasets/D1.txt" 3 -1 -1 $2
echo   "\nDataset 2 : \n"
./exec "2" "Datasets/D2.txt" 3 -1 -1 $2
echo   "\nDataset 3 : \n"
./exec "3" "Datasets/D3.txt" 3 -1 -1 $2
echo   "\nDataset 4 : \n"
./exec "4" "Datasets/D4.txt" 3 -1 -1 $2
echo   "\nDataset 5 : \n"
./exec "5" "Datasets/D5.txt" 3 -1 -1 $2
echo   "\nDataset 6 : \n"
./exec "6" "Datasets/D6.txt" 3 -1 -1 $2
echo   "\nDataset 7 : \n"
./exec "7" "Datasets/D7.txt" 3 -1 -1 $2
echo   "\nDataset 8 : \n"
./exec "8" "Datasets/D8.txt" 3 -1 -1 $2
echo   "\nDataset 9 : \n"
./exec "9" "Datasets/D9.txt" 3 -1 -1 $2
echo   "\nDataset 10 : \n"
./exec "10" "Datasets/D10.txt" 3 -1 -1 $2


echo   "\n \nExecuting Parallel OpenMP Coarse-Grained Program : \n"
gcc -fopenmp -o exec main_OMP_Coarse.c

echo   "\nDataset 9 : \n"
./exec "9" "Datasets/D9.txt" 3 -1 -1 $2
echo   "\nDataset 10 : \n"
./exec "10" "Datasets/D10.txt" 3 -1 -1 $2


echo   "\nExecuting Parallel POSIX Coarse-Grained Program : \n"
gcc main_PTH_Coarse.c -o exec -lpthread

echo   "\nDataset 9 : \n"
./exec "9" "Datasets/D9.txt" 3 -1 -1 $2
echo   "\nDataset 10 : \n"
./exec "10" "Datasets/D10.txt" 3 -1 -1 $2

echo   "\n \nExecuting Parallel POSIX Fine-Grained Program : \n"
gcc main_PTH_Fine.c -o exec -lpthread

echo   "\nDataset 1 : \n"
./exec "1" "Datasets/D1.txt" 3 -1 -1 $2
echo   "\nDataset 2 : \n"
./exec "2" "Datasets/D2.txt" 3 -1 -1 $2
echo   "\nDataset 3 : \n"
./exec "3" "Datasets/D3.txt" 3 -1 -1 $2
echo   "\nDataset 4 : \n"
./exec "4" "Datasets/D4.txt" 3 -1 -1 $2
echo   "\nDataset 5 : \n"
./exec "5" "Datasets/D5.txt" 3 -1 -1 $2
echo   "\nDataset 6 : \n"
./exec "6" "Datasets/D6.txt" 3 -1 -1 $2
echo   "\nDataset 7 : \n"
./exec "7" "Datasets/D7.txt" 3 -1 -1 $2
echo   "\nDataset 8 : \n"
./exec "8" "Datasets/D8.txt" 3 -1 -1 $2
echo   "\nDataset 9 : \n"
./exec "9" "Datasets/D9.txt" 3 -1 -1 $2
echo   "\nDataset 10 : \n"
./exec "10" "Datasets/D10.txt" 3 -1 -1 $2

rm exec
