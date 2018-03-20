/*
Ruslan Abramov 306847393
Ravid Anbary 203373360
*/
#include "Header.h"

const int numberOfThreads = 512;
const int numberOfBlocks = RAN_ARRAY_SIZE / numberOfThreads;

__global__ void firstKernel(int* randomNumbersArray, int *histogram, int randomArraySize)
{
	int arrayIndex = threadIdx.x + blockIdx.x * blockDim.x;
	//There is no actual error with atomicAdd it's just an Intellisense warning caused with vs and cuda intigration
	//we are using atomic add in  order to avoid race conditions
	if (arrayIndex < randomArraySize)
		atomicAdd( &histogram[randomNumbersArray[arrayIndex]],1);
}




// Helper function for using CUDA to add vectors in parallel.
cudaError_t histogramWithCuda(int* randomNumbersArray, int *histogram, int randomArraySize, int histogramSize)
{
	int *dev_randomNumbersArray = 0;
	int *dev_histogram = 0;

	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	int device;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	// Allocate GPU buffers for two vectors (one input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_randomNumbersArray, randomArraySize * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_histogram, histogramSize * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}



	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_randomNumbersArray, randomNumbersArray, randomArraySize * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}



	// Launch a kernel on the GPU with 512 threads for each element.
	//There is no actual error with kernel call it's just an Intellisense warning caused with vs and cuda intigration
	firstKernel << <numberOfBlocks, numberOfThreads >> > (dev_randomNumbersArray, dev_histogram, randomArraySize);
	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}

	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(histogram, dev_histogram, histogramSize * sizeof(int), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}



Error:
	cudaFree(dev_randomNumbersArray);
	cudaFree(dev_histogram);

	return cudaStatus;
}

