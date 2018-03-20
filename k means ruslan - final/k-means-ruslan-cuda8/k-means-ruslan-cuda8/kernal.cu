
#include "DataTypes.h"
#include "cudaHeader.h"
#include <stdio.h>
#include <math.h>


__device__ double my_atomicAdd(double* address, double val);
__device__ double distanceBetween2Points(point2D p1, point2D p2);
__device__ int theSameCluster(clusteredPoint p1, clusteredPoint p2);


__device__ double my_atomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull =
		(unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,
			__double_as_longlong(val +
				__longlong_as_double(assumed)));
	} while (assumed != old);
	return __longlong_as_double(old);
}

__device__ double distanceBetween2Points(point2D p1, point2D p2)
{
	return powf(p1.x - p2.x, 2) + powf(p1.y - p2.y, 2);
}

__device__ int theSameCluster(clusteredPoint p1, clusteredPoint p2)
{
	return p1.myClusterId == p2.myClusterId ? 1 : NULL;
}

//kernals

__global__ void  calculateDistanceFromCenterOfMass(clusteredPoint* clusterPointArray, int numberOfPoints, double *distanceFromCenterOfMass, point2D centerOfMass)
{
	int arrayIndex = threadIdx.x + blockIdx.x * blockDim.x;
	if (arrayIndex < numberOfPoints)
	{
		distanceFromCenterOfMass[arrayIndex] = distanceBetween2Points(clusterPointArray[arrayIndex].point.initCoords, centerOfMass);
	}

}

//calculates the new point position
__global__ void pointMovmentInTime(clusteredPoint* clusterPointArray, int numberOfPoints, double time)
{
	int arrayIndex = threadIdx.x + blockIdx.x * blockDim.x;
	if (arrayIndex < numberOfPoints)
	{
		clusterPointArray[arrayIndex].point.initCoords.x += time * clusterPointArray[arrayIndex].point.initCoords.x;
		clusterPointArray[arrayIndex].point.initCoords.y += time * clusterPointArray[arrayIndex].point.initCoords.y;
	}
}
//calculates an array of farest distances between the points and a cluster
__global__ void distanceFromCluster(clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int cluster)
{
	int arrayIndex = threadIdx.x + blockIdx.x * blockDim.x;
	if (arrayIndex < numberOfPoints)
	{
		//distance between point and it's cluster
		double myCluster_distance = cluster == 0 ? -1 : distanceBetween2Points(clusterPointArray[arrayIndex].point.initCoords, clusterPointArray[arrayIndex].myCentroid);

		//distance between point and a cluster
		double distance = distanceBetween2Points(clusterPointArray[arrayIndex].point.initCoords, clustersArray[cluster].myPoint);
		//the point is closer to another cluster
		if (distance < myCluster_distance || myCluster_distance  < 0)
		{

			clusterPointArray[arrayIndex].myCentroid.x = clustersArray[cluster].myPoint.x;
			clusterPointArray[arrayIndex].myCentroid.y = clustersArray[cluster].myPoint.y;
			clusterPointArray[arrayIndex].myClusterId = cluster;
			//we remeber the cluster index
		}
	}
}
//adds each cluster it's new point
__global__ void addPointToCluster(clusteredPoint* clusterPointArray, int numberOfPoints, newClusterCenter *tempClusterCenter)
{
	int arrayIndex = threadIdx.x + blockIdx.x * blockDim.x;
	if (arrayIndex < numberOfPoints)
	{

		//this is the critical section where we add the coordinates and decide on the 2 far most points
		my_atomicAdd(&(tempClusterCenter[clusterPointArray[arrayIndex].myClusterId].accumulatedDistance.x), clusterPointArray[arrayIndex].point.initCoords.x);
		my_atomicAdd(&(tempClusterCenter[clusterPointArray[arrayIndex].myClusterId].accumulatedDistance.y), clusterPointArray[arrayIndex].point.initCoords.y);
		atomicAdd(&(tempClusterCenter[clusterPointArray[arrayIndex].myClusterId].numberOfPointsInCluster), 1);
	}
}


//calcuates the biggest distance between points in the same cluster
__global__ void distanceBetweenPointsInCluster(clusteredPoint* clusterPointArray, int numberOfPoints, double *biggestDistance, int pointIndex)
{
	int arrayIndex = threadIdx.x + blockIdx.x * blockDim.x;
	if (arrayIndex < numberOfPoints)
	{
		if (theSameCluster(clusterPointArray[arrayIndex], clusterPointArray[pointIndex]) != NULL)
		{
			double distance = distanceBetween2Points(clusterPointArray[arrayIndex].point.initCoords, clusterPointArray[pointIndex].point.initCoords);
			biggestDistance[arrayIndex] = distance > biggestDistance[arrayIndex] ? distance : biggestDistance[arrayIndex];
		}
	}
}


cudaError_t compute_distanceFromCenterOfMass(clusteredPoint* clusterPointArray, int numberOfPoints, point2D centerOfMass, double *distanceFromCenterOfMass)
{
	clusteredPoint* dev_clusterPointArray = NULL;
	double* dev_distanceFromCenterOfMass = NULL;
	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	// Allocate GPU buffers  .
	
	cudaStatus = cudaMalloc((void**)&dev_clusterPointArray, numberOfPoints * sizeof(clusteredPoint));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	

	cudaStatus = cudaMalloc((void**)&dev_distanceFromCenterOfMass, numberOfPoints * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	
	// Copy input vectors from host memory to GPU buffers.

	cudaStatus = cudaMemcpy(dev_clusterPointArray, clusterPointArray, numberOfPoints * sizeof(clusteredPoint), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	

	calculateDistanceFromCenterOfMass << <numberOfThreads, numberOfPoints / numberOfThreads + 1 >> >(dev_clusterPointArray, numberOfPoints, dev_distanceFromCenterOfMass, centerOfMass);

	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(distanceFromCenterOfMass, dev_distanceFromCenterOfMass, numberOfPoints * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(clusterPointArray, dev_clusterPointArray, numberOfPoints * sizeof(clusteredPoint), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}


Error:
	cudaFree(dev_clusterPointArray);
	cudaFree(dev_distanceFromCenterOfMass);
	return cudaStatus;
}



cudaError_t computeNewClusterCenter(clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int numberOfClusters, newClusterCenter *tempClusterCenter )
{
	clusteredPoint* dev_clusterPointArray = NULL;
	cluster* dev_clustersArray = NULL;
	newClusterCenter* dev_tempClusterCenter = NULL;

	cudaError_t cudaStatus;
	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	// Allocate GPU buffers  .
	cudaStatus = cudaMalloc((void**)&dev_clusterPointArray, numberOfPoints * sizeof(clusteredPoint));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_clustersArray, numberOfClusters * sizeof(cluster));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_tempClusterCenter, numberOfClusters * sizeof(newClusterCenter));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}


	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_clusterPointArray, clusterPointArray, numberOfPoints * sizeof(clusteredPoint), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_clustersArray, clustersArray, numberOfClusters * sizeof(cluster), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	//cheking n points with k clusters
	for (int counter = 0; counter < numberOfClusters; counter++)
	{
		distanceFromCluster << <numberOfThreads, numberOfPoints / numberOfThreads + 1 >> >(dev_clusterPointArray, numberOfPoints, dev_clustersArray, counter);
		cudaStatus = cudaDeviceSynchronize();
	}
	//adding the points to their cluster
	//this is a critical area so we use locks
	addPointToCluster << <numberOfThreads, numberOfPoints / numberOfThreads + 1 >> >(dev_clusterPointArray, numberOfPoints, dev_tempClusterCenter);
	cudaStatus = cudaDeviceSynchronize();
	
	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(tempClusterCenter, dev_tempClusterCenter, numberOfClusters * sizeof(newClusterCenter), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(clusterPointArray, dev_clusterPointArray, numberOfPoints * sizeof(clusteredPoint), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}



Error:
	cudaFree(dev_clusterPointArray);
	cudaFree(dev_clustersArray);
	cudaFree(dev_tempClusterCenter);

	return cudaStatus;
}

cudaError_t movePointInTime(clusteredPoint* clusterPointArray, int numberOfPoints, double time)
{
	clusteredPoint* dev_clusterPointArray = NULL;

	cudaError_t cudaStatus;
	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	// Allocate GPU buffers  .
	cudaStatus = cudaMalloc((void**)&dev_clusterPointArray, numberOfPoints * sizeof(clusteredPoint));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}


	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_clusterPointArray, clusterPointArray, numberOfPoints * sizeof(clusteredPoint), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}


	//moving the points
	pointMovmentInTime << <numberOfThreads, numberOfPoints / numberOfThreads + 1 >> >(dev_clusterPointArray, numberOfPoints, time);


	// Copy output vector from GPU buffer to host memory.

	cudaStatus = cudaMemcpy(clusterPointArray, dev_clusterPointArray, numberOfPoints * sizeof(clusteredPoint), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
Error:
	cudaFree(dev_clusterPointArray);


	return cudaStatus;
}

/*this method found out to really slow --> therfore i have decided to not use it but i still want to know why isn't it good*/
cudaError_t computeDistanceBetweenPointsInNewClusterCenter(clusteredPoint* clusterPointArray, int numberOfPoints, double *maxDistance)
{
	clusteredPoint* dev_clusterPointArray = NULL;
	double* dev_maxDistance = NULL;


	cudaError_t cudaStatus;
	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	// Allocate GPU buffers  .
	cudaStatus = cudaMalloc((void**)&dev_clusterPointArray, numberOfPoints * sizeof(clusteredPoint));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_maxDistance, numberOfPoints * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_clusterPointArray, clusterPointArray, numberOfPoints * sizeof(clusteredPoint), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	for (int counter = 0; counter < numberOfPoints; counter++)
	{
		distanceBetweenPointsInCluster << < numberOfThreads, numberOfPoints / numberOfThreads + 1 >> > (dev_clusterPointArray, numberOfPoints, dev_maxDistance, counter);
		printf("%d\n", counter);
		fflush(stdout);
	}


	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(maxDistance, dev_maxDistance, numberOfPoints * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
Error:
	cudaFree(dev_clusterPointArray);
	cudaFree(dev_maxDistance);

	return cudaStatus;
}
