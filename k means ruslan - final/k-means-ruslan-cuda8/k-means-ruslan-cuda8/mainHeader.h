
#include "DataTypes.h"
#include "globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

/**************K Mean*****************/
/**************init*******************/
void initClusters(clusteredPoint* clusterPoint, int numberOfPoints, cluster* clustersArray, int numberOfClusters);
void findMaxDistanceIndex(clusteredPoint* clusterPointArray , int numberOfPoints, double *distanceFromCenterOfMass, int *maxIndex, int *pointsAdded, int clusterToAdd);
void addNewCluster(clusteredPoint* clusterPointArray, int numberOfPoints, double *distanceFromCenterOfMass, cluster* clustersArray, int clusterToAdd, int *pointsAdded);
void recalculateCenterOfMass(cluster* clustersArray, int numberOfClusters, point2D *centerOfMass);
void compute_distanceFromCenterOfMass_with_cuda(clusteredPoint* clusterPointArray, int numberOfPoints, point2D centerOfMass, double *distanceFromCenterOfMass);
int pointAllreadyAdded(clusteredPoint* clusterPointArray, int *pointsAdded, int point, int clusterToAdd);
void computeFirstClusters(clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int numberOfClusters, int *pointsAdded);
void makeFirstRandomCluster(clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int *pointsAdded);
double distanceBetween(point2D point1, point2D point2);
double calculateClusterQualitiy(cluster* clustersArray);

/************itirations*************/
newClusterCenter* k_means_iteration(clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int numberOfClusters, int *flag_point_moved_between_clusters);
int didAPointMoved(clusteredPoint* newClusterPointArray, clusteredPoint* oldClusterPointArray, int numberOfPoints);

void k_means(int myid, clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int numberOfClusters, double* tempQualitiy, double* time, int numprocs, MPI_Status* statusP, MPI_Datatype mpiCluster, MPI_Datatype mpiNewClusterCenter, MPI_Datatype mpiClusterdPoint);
void addClusterData(newClusterCenter* mastersNewClusters, newClusterCenter* temp, int numberOfClusters);
void calculateNewClusterArray(newClusterCenter* mastersNewClusters, cluster* clustersArray, int numberOfClusters, clusteredPoint* points, double* distance, int numberOfPoints);
void masterCalculatesNewQualitiy(cluster* clustersArray, double *qualitiy);
void checkTermination(double tempQualitiy, int counter, double time, int* keepIteratingP);
void movePoints(clusteredPoint* clusterPointArray, int numberOfPoints, double time);

void masterCreateNewClusterCentersWithDiameter(newClusterCenter* mastersNewClusters, clusteredPoint* mastersPoints, cluster* clustersArray, int numberOfClusters, int numprocs, int masterNumberOfPoints, int slaveNumberOfPoints, MPI_Status* statusP, MPI_Datatype mpiNewClusterCenter, MPI_Datatype mpiClusterdPoint);
void setDiameter(cluster* clustersArray, int numberOfClusters, clusteredPoint* points, double* distance, int numberOfPoints);

void masterCreateNewClusterCentersNoDiameter(newClusterCenter* mastersNewClusters, cluster* clustersArray, int numberOfClusters, int numprocs, MPI_Status* statusP, MPI_Datatype mpiNewClusterCenter);
void computeDistanceBetweenPointsInNewClusterCenter_with_cuda(clusteredPoint* tempPoints, int numberOfPoints, double *distance);
/*****init and final*****/
void startMpi(int *argc, char **argv[], int* myidP, int* numprocsP, MPI_Datatype *point2DMPIType, MPI_Datatype *velocity2DMPIType, MPI_Datatype *movingPoint2DMPIType, MPI_Datatype *clusteredPoint, MPI_Datatype *newClusterCenter, MPI_Datatype *cluster);
void finalize(clusteredPoint* clusterPointArray, cluster* clustersArray, MPI_Datatype *point2DMPIType, MPI_Datatype *velocity2DMPIType, MPI_Datatype *movingPoint2DMPIType, MPI_Datatype *clusteredPoint, MPI_Datatype *newClusterCenter, MPI_Datatype *cluster);
/******memory******/
void makeClusteredPoint(clusteredPoint** clusterPointArrayPointer, int size);
void makeClustersArray(cluster** clustersArrayPointer, int size);
void makeDoubleArray(double** array, int size);
void makenewClusterCenterArray(newClusterCenter** array, int size);
void makeIntArray(int** array, int size);

void copyOldClusterArray(clusteredPoint** oldClusterPointArrayClusters, clusteredPoint* clusterPointArray, int numberOfPoints);
void shiftArray(clusteredPoint* points, int offset, int numberOfPoints);
void copyMastersPoints(clusteredPoint* temp, clusteredPoint* mester, int numberOfPoints);
/***************mpi process commuincation*****************/
void masterSendPointsToSlaves(clusteredPoint* clusterPointArray, int numberOfPoints, int numprocs, MPI_Datatype mpiClusterdPoint, int *masterProcessNumberOfPoints);
void slavesRecivePointsFromMaster(clusteredPoint* clusterPointArray, int chunckSize, MPI_Datatype mpiClusterdPoint, MPI_Status *statusP);
void slavesSend_newClusterCenter(newClusterCenter* newClusters, int numberOfClusters, MPI_Datatype mpiNewClusterCenter);
void slavesSend_Points(clusteredPoint* clusterPointArray, int numberOfPoints, MPI_Datatype mpiClusterdPoint);///
void broadCastInitialDataToSlaves(int *numberOfPointsForProcess, MPI_Comm comm);
void SendClustersToSlave(cluster* clustersArray, int numberOfClusters, MPI_Datatype mpiCluster);
void broadCastSlavesLoopVar(int* keepIterating, int* updatedCounter, double* updatedTime);
/*********error handling**********/
void ERROR(const char* MSG);
/*************file handling******************/
void scanGlobalsFromFile(FILE * file, clusteredPoint** clusterPointArrayPointer, cluster** clustersArrayPointer);
void scanMovingPoint(FILE * file, float* x, float* y, float* xVelocity, float* yVelocity);
void readPointsFromFile(char* filePath, clusteredPoint** clusterPointArrayPointer, cluster** clustersArrayPointer);
void writeToFile(char* outputFilePath, cluster* clustersArray, int numberOfClusters, double* tempQualitiy, double time, double startTime);
void writeClusterCenters(FILE* file, cluster* clustersArray, int numberOfClusters);
/***************mpi data types************/
void createPoint2DMpiDataType(MPI_Datatype *point2DMPIType);
void createVelocity2DMpiDataType(MPI_Datatype *velocity2DMPIType);
void createMovingPoint2DMpiDataType(MPI_Datatype *movingPoint2DMPIType, MPI_Datatype point2DMPIType, MPI_Datatype velocity2DMPIType);
void createClusteredPointMpiDataType(MPI_Datatype *clusteredPoint, MPI_Datatype movingPoint2DMPIType, MPI_Datatype point2DMPIType);
void createNewClusterCenterMpiDataType(MPI_Datatype *newClusterCenter, MPI_Datatype point2DMPIType);
void createClusterMpiDataType(MPI_Datatype *cluster, MPI_Datatype point2DMPIType);