/*
Ruslan Abramov 306847393
Ravid Anbary 203373360
*/
#include "Header.h"
#include <mpi.h>
#include <omp.h>
#include <time.h>


const int RGB = 256;
const int MASTER = 0;
const int SLAVE = 1;

const char* ERROR_MSG_1 = "\n Memory could not be allocated --> aborting";

void makeIntArray(int** array, int size);
/***********************Master related******************************/
void mastersWork(int* randomNumbers, int* histogram, MPI_Status * status);
void initRandomArray(int* array, int arraySize);
void histogramOMP(int* array, int size, int* histogram);
void addHistograms(int *mastersHistogram, int *slavesHistogram, int size);
void printHistogram(int *histogram, int size);

/**************************slave related*********************************/
void slaveWork(int *randomNumbers, int randomNumbersSize, int *histogram, int histogramSize, MPI_Status *status);
void computeWithCuda(int* randomNumbersArray, int *histogram, int randomArraySize, int histogramSize);

int main(int argc, char *argv[])
{
	
	int   myid, numprocs;
	int *randomNumbers = NULL;
	int histogram[RGB] = { 0 };
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	if (myid == MASTER)
	{
		mastersWork(randomNumbers, histogram, &status);
		printHistogram(histogram , RGB);
	}
	if (myid == SLAVE){
		slaveWork(randomNumbers, RAN_ARRAY_SIZE / 2, histogram, RGB, &status);
	}

	MPI_Finalize();
	free(randomNumbers);
	
	return 0;
}

void initRandomArray(int* array, int arraySize)
{
	srand((unsigned int)time(NULL));   // should only be called once
#pragma omp parallel for
	for (int index = 0; index < RAN_ARRAY_SIZE; index++)
		array[index] = rand() % RGB;      //  a pseudo-random integer between 0 and 255
}

void histogramOMP(int* array, int size, int* histogram)
{
	
#pragma omp parallel 
	{	//each thread updates his histogram
		int privateHistoram[RGB] = { 0 };
#pragma omp for
		for (int index = 0; index < size; index++)
			privateHistoram[array[index]]++;

#pragma omp critical
		for (int index = 0; index < RGB; index++)
			histogram[index] += privateHistoram[index];
		
	}
	
}

void mastersWork(int* randomNumbers, int* histogram, MPI_Status * status)
{
	int slavesHistogram[RGB];
	makeIntArray(&randomNumbers, RAN_ARRAY_SIZE);
	initRandomArray(randomNumbers, RAN_ARRAY_SIZE);
	MPI_Send(&randomNumbers[RAN_ARRAY_SIZE / 2], RAN_ARRAY_SIZE / 2, MPI_INT, SLAVE, 0, MPI_COMM_WORLD);
	histogramOMP(randomNumbers, RAN_ARRAY_SIZE / 2, histogram);
	MPI_Recv(slavesHistogram, RGB, MPI_INT, SLAVE, 0, MPI_COMM_WORLD, status);
	addHistograms(histogram, slavesHistogram, RGB);
}

void addHistograms(int *mastersHistogram, int *slavesHistogram, int size)
{
#pragma omp parallel for
	for (int index = 0; index < size; index++)
		mastersHistogram[index] += slavesHistogram[index];
}

void makeIntArray(int** array, int size){
	*array = (int*)malloc(size * sizeof(int));
	if (array == NULL){
		/* Memory could not be allocated*/
		printf("%s", ERROR_MSG_1);
		MPI_Finalize();
		exit(1);
	}
}

void computeWithCuda(int* randomNumbersArray, int *histogram, int randomArraySize, int histogramSize)
{
	cudaError_t  cudaStatus = histogramWithCuda(randomNumbersArray, histogram, randomArraySize, histogramSize);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "histogramWithCuda failed!");
			MPI_Finalize();
			exit(1);
		}

		cudaStatus = cudaDeviceReset();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceReset failed!");
			MPI_Finalize();
			exit(1);
		}
}

void printHistogram(int *histogram, int size)
{
	printf("\n*****************************************\n");
	int sum = 0;
	for (int index = 0; index < size; index++)
	{
		printf("(%d) --> %d\n", index + 1, histogram[index]);
		sum += histogram[index];
	}
	printf("*******total numbers (%d) \n", sum);
	printf("*******total numbers (%d) - 8 * 1024 * 1024 (8MB)= %d \n", sum , sum - RAN_ARRAY_SIZE);
}

void slaveWork(int *randomNumbers , int randomNumbersSize , int *histogram , int histogramSize , MPI_Status *status)
{
	makeIntArray(&randomNumbers, randomNumbersSize);
	MPI_Recv(randomNumbers, randomNumbersSize, MPI_INT, MASTER, 0, MPI_COMM_WORLD, status);
	//histogram using cuda
	computeWithCuda(randomNumbers, histogram, randomNumbersSize, histogramSize);
	MPI_Send(histogram, histogramSize, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
}