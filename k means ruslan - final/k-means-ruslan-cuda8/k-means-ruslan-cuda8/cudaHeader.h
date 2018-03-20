#include "cuda_runtime.h"
#include "device_launch_parameters.h"


const int numberOfThreads = 512;

cudaError_t compute_distanceFromCenterOfMass(clusteredPoint* clusterPointArray, int numberOfPoints, point2D centerOfMass, double *distanceFromCenterOfMass);
cudaError_t computeNewClusterCenter(clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int numberOfClusters, newClusterCenter *tempClusterCenter);
cudaError_t computeDistanceBetweenPointsInNewClusterCenter(clusteredPoint* clusterPointArray, int numberOfPoints, double *maxDistance);
cudaError_t movePointInTime(clusteredPoint* clusterPointArray, int numberOfPoints, double time);
