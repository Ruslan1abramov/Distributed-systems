/*
	Ruslan Abramov 306847393
	Ravid Anbary 203373360
*/
#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
const int RAN_ARRAY_SIZE = 1024 * 1024 * 8;

cudaError_t  histogramWithCuda(int* randomNumbersArray, int *histogram, int randomArraySize, int histogramSize);