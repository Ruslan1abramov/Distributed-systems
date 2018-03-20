#include "mainHeader.h"
#include "cudaHeader.h"


char* filePath = "C:\\Users\\afeka\\Desktop\\input_1.txt";
char* outfilePath = "C:\\Users\\afeka\\Desktop\\randomPoints_res.txt";


void main(int argc, char *argv[])
{

	clusteredPoint* clusterPointArray = NULL;
	cluster* clustersArray = NULL;
	double calculatedQualitiy = 0;
	int myid, numprocs;
	int numberOfPointsForProcess;
	double time = 0;
	/*stratring MPI*/
	MPI_Datatype point2DMPIType, velocity2DMPIType, movingPoint2DMPIType, clusteredPointsMPI, newClusterCenterMPI, clusterMPI;
	MPI_Status statusP;
	startMpi(&argc, &argv, &myid, &numprocs, &point2DMPIType, &velocity2DMPIType, &movingPoint2DMPIType, &clusteredPointsMPI, &newClusterCenterMPI, &clusterMPI);

	double sTime = 0;

	if (myid == MASTER)
	{
		sTime = MPI_Wtime();
		readPointsFromFile(filePath, &clusterPointArray, &clustersArray);
		initClusters(clusterPointArray, N, clustersArray, K);
		numberOfPointsForProcess = N / numprocs;
	}
	broadCastInitialDataToSlaves(&numberOfPointsForProcess, MPI_COMM_WORLD);

	if (myid == MASTER)
	{
		masterSendPointsToSlaves(clusterPointArray, N, numprocs, clusteredPointsMPI, &numberOfPointsForProcess);
	}
	else
	{
		makeClusteredPoint(&clusterPointArray, numberOfPointsForProcess);
		slavesRecivePointsFromMaster(clusterPointArray, numberOfPointsForProcess, clusteredPointsMPI, &statusP);
		makeClustersArray(&clustersArray, K);
	}
	
	k_means(myid, clusterPointArray, numberOfPointsForProcess, clustersArray, K, &calculatedQualitiy, &time, numprocs, &statusP, clusterMPI, newClusterCenterMPI, clusteredPointsMPI);
	if (myid == MASTER)
	{
		writeToFile(outfilePath, clustersArray, K, &calculatedQualitiy, time - dT, sTime);
	}

	finalize(clusterPointArray, clustersArray, &point2DMPIType, &velocity2DMPIType, &movingPoint2DMPIType, &clusteredPointsMPI, &newClusterCenterMPI, &clusterMPI);
	printf("Process # %d ----> DONE\n" , myid);
	fflush(stdout);
}

/*creating an array of clusteredPoint*/
void makeClusteredPoint(clusteredPoint** clusterPointArrayPointer, int size) {
	*clusterPointArrayPointer = (clusteredPoint *)malloc(size * sizeof(clusteredPoint));
	if (*clusterPointArrayPointer == NULL) {
		/* Memory could not be allocated*/
		ERROR(ERROR_MSG_1);
	}
}
/*creating an array of cluster*/
void  makeClustersArray(cluster** clustersArrayPointer, int size) {
	*clustersArrayPointer = (cluster *)malloc(size * sizeof(cluster));
	if (*clustersArrayPointer == NULL) {
		/* Memory could not be allocated*/
		ERROR(ERROR_MSG_1);
	}
}
/*creating an array of doubles*/
void makeDoubleArray(double** array, int size) {
	*array = (double*)malloc(size * sizeof(double));
	if (array == NULL) {
		/* Memory could not be allocated*/
		ERROR(ERROR_MSG_1);
	}
}

/*creating an array of Int*/
void makeIntArray(int** array, int size) {
	*array = (int*)malloc(size * sizeof(int));
	if (array == NULL) {
		/* Memory could not be allocated*/
		ERROR(ERROR_MSG_1);
	}
}

/*copying the old cluster Id*/
/*this will be used in order to decide if there was movement of points between clusters*/
void copyOldClusterArray(clusteredPoint** oldClusterPointArrayClusters, clusteredPoint* clusterPointArray, int numberOfPoints)
{
	makeClusteredPoint(oldClusterPointArrayClusters, numberOfPoints);
	for (int counter = 0; counter < numberOfPoints; counter++)
	{
		(*oldClusterPointArrayClusters)[counter].myClusterId = clusterPointArray[counter].myClusterId;
	}

}
void scanGlobalsFromFile(FILE * file, clusteredPoint** clusterPointArrayPointer, cluster** clustersArrayPointer)
{
	/*file is build in a way that the first line contains 6 information verbials
	all the other lines contain information need to create a movingPoint in each line*/
	//reading number of points
	fscanf_s((FILE*)file, "%d", &N);
	//creating a clusteredPoint array size of N
	makeClusteredPoint(clusterPointArrayPointer, N);
	//reading number of clusters
	fscanf_s((FILE*)file, "%d", &K);
	//creating a centroid array size of K
	makeClustersArray(clustersArrayPointer, K);
	//reading the end of time interval[0, T]
	fscanf_s((FILE*)file, "%lf", &T);
	//reading the number that defines a moments
	fscanf_s((FILE*)file, "%lf", &dT);
	//reading the number of maximum iterations
	fscanf_s((FILE*)file, "%d", &LIMIT);
	//reading the minimum cluster quality
	fscanf_s((FILE*)file, "%lf", &QM);
}

void scanMovingPoint(FILE * file, float* x, float* y, float* xVelocity, float* yVelocity) {
	/*reading a point from the file*/
	fscanf_s((FILE*)file, "%f", x);
	fscanf_s((FILE*)file, "%f", y);
	fscanf_s((FILE*)file, "%f", xVelocity);
	fscanf_s((FILE*)file, "%f", yVelocity);
}
void readPointsFromFile(char* filePath, clusteredPoint** clusterPointArrayPointer, cluster** clustersArrayPointer) {
	FILE * file;
	fopen_s(&file, filePath, "r");
	if (file == NULL)
	{
		printf(ERROR_MSG_2);
		exit(1);
	}
	if (file) {
		scanGlobalsFromFile(file, clusterPointArrayPointer, clustersArrayPointer);

		for (int i = 0; i < N; i++) {
			float x = 0, y = 0;
			float xVelocity = 0, yVelocity = 0;
			scanMovingPoint(file, &x, &y, &xVelocity, &yVelocity);
			clusteredPoint newClusteredPoint;
			newClusteredPoint.point.initCoords.x = x;
			newClusteredPoint.point.initCoords.y = y;
			newClusteredPoint.point.velocity.xVelocity = xVelocity;
			newClusteredPoint.point.velocity.yVelocity = yVelocity;
			newClusteredPoint.myCentroid.x = NULL;
			newClusteredPoint.myCentroid.y = NULL;
			newClusteredPoint.myClusterId = -1;

			(*clusterPointArrayPointer)[i] = //placing a new clusteredPoint
				newClusteredPoint;
		}
		fclose(file);
	}
	//complaxity of O(N)
}

void initClusters(clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int numberOfClusters) {
	/*we choose the first centroids in the following way:
	1. first we choose a random moving point
	2. the next centroid will be the the farthest point from the centroid
	3. the next centroid will be the farthest point from both of the centroids
	4. and so on until we calculate k centroids
	this way we provide a better prediction and thus more likely to convarge faster
	*/
	int *pointsAdded = NULL;
	makeIntArray(&pointsAdded, numberOfClusters);

	makeFirstRandomCluster(clusterPointArray, numberOfPoints, clustersArray, pointsAdded);
	computeFirstClusters(clusterPointArray, numberOfPoints, clustersArray, numberOfClusters, pointsAdded);

	free(pointsAdded);
}


/*choosing one random point as first cluster*/
void makeFirstRandomCluster(clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int *pointsAdded) {
	srand((unsigned int)time(NULL));   // should only be called once
	int pointIndex = rand() % numberOfPoints;
	clusteredPoint cls = clusterPointArray[pointIndex]; //  a pseudo-random integer between 0 and numberOfPoints
	cluster firstCluster;
	firstCluster.myPoint = cls.point.initCoords;
	firstCluster.diamtre = 0;
	clustersArray[0] = firstCluster;
	pointsAdded[0] = pointIndex;
}
/*computing the first clusters*/
void computeFirstClusters(clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int numberOfClusters, int *pointsAdded) {
	double *distanceFromCenterOfMass = NULL;
	makeDoubleArray(&distanceFromCenterOfMass, numberOfPoints);
	point2D centerOfMass = clustersArray[0].myPoint;
	//each iteration we recalculate the cluster center of mass and we pick the farest point from it as the next cluster
	for (int counter = 1; counter < numberOfClusters; counter++)
	{
		recalculateCenterOfMass(clustersArray, counter, &centerOfMass);
		compute_distanceFromCenterOfMass_with_cuda(clusterPointArray, numberOfPoints, centerOfMass, distanceFromCenterOfMass);
		addNewCluster(clusterPointArray, numberOfPoints, distanceFromCenterOfMass, clustersArray, counter, pointsAdded);
	}

	free(distanceFromCenterOfMass);
	//complexity of O(K^2*N) --> i use cuda in order to speed things up
	//this is k^2 because in each iteration I scan the cluster already being added in order to check if the point I want to add was already chosen
}
/*calculating with cuda the distance of each point from center of mass*/
void compute_distanceFromCenterOfMass_with_cuda(clusteredPoint* clusterPointArray, int numberOfPoints, point2D centerOfMass, double *distanceFromCenterOfMass)
{
	cudaError_t  cudaStatus = compute_distanceFromCenterOfMass(clusterPointArray, numberOfPoints, centerOfMass, distanceFromCenterOfMass);
	if (cudaStatus != cudaSuccess)
		ERROR(MY_CUDA_ERROR_MSG_1);

	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess)
		ERROR(MY_CUDA_ERROR_DeviceReset);
}
void addNewCluster(clusteredPoint* clusterPointArray, int numberOfPoints, double *distanceFromCenterOfMass, cluster* clustersArray, int clusterToAdd, int *pointsAdded)
{
	int pointToAdd = 0;
	findMaxDistanceIndex(clusterPointArray , numberOfPoints, distanceFromCenterOfMass, &pointToAdd, pointsAdded, clusterToAdd);
	clustersArray[clusterToAdd].myPoint = clusterPointArray[pointToAdd].point.initCoords;
	clustersArray[clusterToAdd].diamtre = 0;
	pointsAdded[clusterToAdd] = pointToAdd;
}
void findMaxDistanceIndex(clusteredPoint* clusterPointArray,int numberOfPoints, double *distanceFromCenterOfMass, int *maxIndex, int *pointsAdded, int clusterToAdd)
{
	*maxIndex = 0;
	for (int counter = 0; counter < numberOfPoints; counter++)
		if (distanceFromCenterOfMass[*maxIndex] < distanceFromCenterOfMass[counter] && pointAllreadyAdded(clusterPointArray, pointsAdded, counter, clusterToAdd) == 0)
			*maxIndex = counter;
}
/*this method is used in order to prevent from choosing the same points as clusters*/
int pointAllreadyAdded(clusteredPoint* clusterPointArray, int *pointsAdded, int point, int clusterToAdd)
{
	int beenAdded = 0;
	for (int counter = 0; counter < clusterToAdd; counter++)
	{
		if (clusterPointArray[pointsAdded[counter]].point.initCoords.x == clusterPointArray[point].point.initCoords.x && clusterPointArray[pointsAdded[counter]].point.initCoords.y == clusterPointArray[point].point.initCoords.y)
		{
			beenAdded = 1;
			return beenAdded;
		}
	}
	return beenAdded;
}

void recalculateCenterOfMass(cluster* clustersArray, int numberOfClusters, point2D *centerOfMass)
{
	double x = 0, y = 0;
	for (int counter = 0; counter < numberOfClusters; counter++)
	{
		x += clustersArray[counter].myPoint.x;
		y += clustersArray[counter].myPoint.y;
	}
	*centerOfMass = { x / numberOfClusters, y / numberOfClusters };

}
/*calculating distance between 2 points*/
double distanceBetween(point2D point1, point2D point2) {
	return (point1.x - point2.x)*(point1.x - point2.x) + (point1.y - point2.y)*(point1.y - point2.y);
}

double calculateClusterQualitiy(cluster* clustersArray) {
	double* clusterQuality = NULL;
	makeDoubleArray(&clusterQuality, K);
	//inint cluster quality
#pragma omp parallel for
	for (int counter = 0; counter < K; counter++)
		clusterQuality[counter] = 0;
#pragma omp parallel for
	for (int counter = 0; counter < K; counter++) {
		//calculating each clusters quality
		double diamtre = clustersArray[counter].diamtre;
		for (int i = 0; i < K; i++) {
			//for each cluster other then ourself
			if (i != counter) {
				clusterQuality[counter] += diamtre / sqrt(distanceBetween(clustersArray[counter].myPoint, clustersArray[i].myPoint));
			}
		}
		/*
		some complexity notes this method is done in O(k^2) using openMp we speed things up by *numberOfThreads at best
		it is possible to use cuda in this section but because K << N there won't be any added value for cuda
		*/
	}
	double q = 0;
	for (int counter = 0; counter < K; counter++)
		q += clusterQuality[counter];
	free(clusterQuality);
	return q / (K*(K - 1));
}


void ERROR(const char* MSG) {
	printf("%s", MSG);
	MPI_Finalize();
	exit(1);
}

void startMpi(int *argc, char **argv[], int* myidP, int* numprocsP, MPI_Datatype *point2DMPIType, MPI_Datatype *velocity2DMPIType, MPI_Datatype *movingPoint2DMPIType, MPI_Datatype *clusteredPoint, MPI_Datatype *newClusterCenter, MPI_Datatype *cluster) {
	MPI_Init(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, myidP);
	MPI_Comm_size(MPI_COMM_WORLD, numprocsP);
	createPoint2DMpiDataType(point2DMPIType);
	createVelocity2DMpiDataType(velocity2DMPIType);
	createMovingPoint2DMpiDataType(movingPoint2DMPIType, *point2DMPIType, *velocity2DMPIType);
	createClusteredPointMpiDataType(clusteredPoint, *movingPoint2DMPIType, *point2DMPIType);
	createNewClusterCenterMpiDataType(newClusterCenter, *point2DMPIType);
	createClusterMpiDataType(cluster, *point2DMPIType);

}
void broadCastInitialDataToSlaves(int *numberOfPointsForProcess, MPI_Comm comm)
{
	MPI_Bcast(&N, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&K, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&dT, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(numberOfPointsForProcess, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
}

void masterSendPointsToSlaves(clusteredPoint* clusterPointArray, int numberOfPoints, int numprocs, MPI_Datatype mpiClusterdPoint, int *masterProcessNumberOfPoints) {
	//sends each slave N/number of process points
	//the master will also calculate N/number of process points + N - number of process points*N/number of process points
	int offset = numberOfPoints - numprocs * (numberOfPoints / numprocs);
	int chunckSize = numberOfPoints / numprocs;
	*masterProcessNumberOfPoints = offset + *masterProcessNumberOfPoints;
	for (int counter = 1; counter < numprocs; counter++)
		MPI_Send(&clusterPointArray[offset + counter*chunckSize], chunckSize, mpiClusterdPoint, counter, 0, MPI_COMM_WORLD);
}

void slavesRecivePointsFromMaster(clusteredPoint* clusterPointArray, int chunckSize, MPI_Datatype mpiClusterdPoint, MPI_Status *statusP) {
	MPI_Recv(clusterPointArray, chunckSize, mpiClusterdPoint, MASTER, 0, MPI_COMM_WORLD, statusP);
}

void SendClustersToSlave(cluster* clustersArray, int numberOfClusters, MPI_Datatype mpiCluster) {
	//broadcasting the clusters from the master to everyone
	MPI_Bcast(clustersArray, numberOfClusters, mpiCluster, MASTER, MPI_COMM_WORLD);
}
/*checking if a point switched cluster*/
int didAPointMoved(clusteredPoint* newClusterPointArray, clusteredPoint* oldClusterPointArray, int numberOfPoints)
{
	int flag_moved = 0;
	for (int counter = 0; counter < numberOfPoints; counter++)
		if (newClusterPointArray[counter].myClusterId != oldClusterPointArray[counter].myClusterId)
			flag_moved++;
	return flag_moved;
}
newClusterCenter* k_means_iteration(clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int numberOfClusters, int *flag_point_moved_between_clusters) {
	
	//copying the old cluster id
	clusteredPoint* oldClusterPointArray = NULL;
	copyOldClusterArray(&oldClusterPointArray, clusterPointArray , numberOfPoints);
	
	newClusterCenter* temp = NULL;
	makenewClusterCenterArray(&temp, numberOfClusters);
	//using cuda to calculate next clusters
	cudaError_t  cudaStatus = computeNewClusterCenter(clusterPointArray, numberOfPoints, clustersArray, numberOfClusters, temp );
	if (cudaStatus != cudaSuccess)
		ERROR(MY_CUDA_ERROR_MSG_1);

	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess)
		ERROR(MY_CUDA_ERROR_DeviceReset);

	//Checking if points moved between clusters
	*flag_point_moved_between_clusters = didAPointMoved(clusterPointArray, oldClusterPointArray, numberOfPoints);
	free(oldClusterPointArray);
	return temp;
	//O(K*N) using cuda we can reduce it, each point will search with K threads
}

void slavesSend_newClusterCenter(newClusterCenter* newClusters, int numberOfClusters, MPI_Datatype mpiNewClusterCenter) {
	//slaves send their new cluster centers
	MPI_Send(newClusters, numberOfClusters, mpiNewClusterCenter, MASTER, 0, MPI_COMM_WORLD);
	//O(k)
}

void slavesSend_Points(clusteredPoint* clusterPointArray, int numberOfPoints, MPI_Datatype mpiClusterdPoint)
{
	//slaves send their points
	MPI_Send(clusterPointArray, numberOfPoints, mpiClusterdPoint, MASTER, 0, MPI_COMM_WORLD);
	//O(n)
}


void masterCreateNewClusterCentersWithDiameter(newClusterCenter* mastersNewClusters, clusteredPoint* mastersPoints, cluster* clustersArray, int numberOfClusters, int numprocs, int masterNumberOfPoints, int slaveNumberOfPoints, MPI_Status* statusP, MPI_Datatype mpiNewClusterCenter, MPI_Datatype mpiClusterdPoint) {
	//master receives from slaves their new cluster centers and calculates the new cluster centers
	newClusterCenter* temp = NULL;
	makenewClusterCenterArray(&temp, numberOfClusters);
	for (int counter = 1; counter < numprocs; counter++)
	{
		MPI_Recv(temp, numberOfClusters, mpiNewClusterCenter, counter, 0, MPI_COMM_WORLD, statusP);
		addClusterData(mastersNewClusters, temp, numberOfClusters);
	}

	//gathering the points from the slaves and creating an array of points in order to calculate the diameter
	clusteredPoint* tempPoints = NULL;
	makeClusteredPoint(&tempPoints, N);
	//master recivies from slaves their points
	for (int counter = 1; counter < numprocs; counter++)
	{
		MPI_Recv(tempPoints, slaveNumberOfPoints, mpiClusterdPoint, counter, 0, MPI_COMM_WORLD, statusP);
		int offset = (counter - 1) * slaveNumberOfPoints + masterNumberOfPoints;
		shiftArray(tempPoints, offset, slaveNumberOfPoints);
	}
	copyMastersPoints(tempPoints, mastersPoints, masterNumberOfPoints);

	//calculates an array contains  the square distance between each point and its far most point in the cluster
	double *distance = NULL;
	makeDoubleArray(&distance, N);
	computeDistanceBetweenPointsInNewClusterCenter_with_cuda(tempPoints, N, distance);
	calculateNewClusterArray(mastersNewClusters, clustersArray, numberOfClusters, tempPoints, distance, N);

	free(temp);
	free(tempPoints);
	free(distance);
}

/*this method call cuda that computes the distance between each point and it far most Point inside the cluster -- it had bad time results*/
/*i've updated this method to work woth omp it is much faster this way*/
void computeDistanceBetweenPointsInNewClusterCenter_with_cuda(clusteredPoint* tempPoints , int numberOfPoints , double *distances)
{
	double st = MPI_Wtime();
	printf("Staring to calculate distances ---\n");
	fflush(stdout);
	
#pragma omp parallel for
	for (int out_counter = 0; out_counter < numberOfPoints; out_counter++)
	{
		for (int innerCounter = out_counter + 1; innerCounter < numberOfPoints; innerCounter++)
		{
			if (tempPoints[out_counter].myClusterId == tempPoints[innerCounter].myClusterId )
			{
				double distance = distanceBetween(tempPoints[out_counter].point.initCoords, tempPoints[innerCounter].point.initCoords);
				distances[out_counter] = distance > distances[out_counter] ? distance : distances[out_counter];
			}
		}
	}
	printf("It took %lf to calculate distances ---\n" , MPI_Wtime() - st);
	fflush(stdout);
	/*
	cudaError_t  cudaStatus = computeDistanceBetweenPointsInNewClusterCenter(tempPoints, N, distances);
	//O(n^2) - using cuda will reduce it
	//I must admit that while i use cuda for this part it is still really slow
	if (cudaStatus != cudaSuccess)
		ERROR(MY_CUDA_ERROR_MSG_1);

	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess)
		ERROR(MY_CUDA_ERROR_DeviceReset);
		*/
}
void masterCreateNewClusterCentersNoDiameter(newClusterCenter* mastersNewClusters, cluster* clustersArray, int numberOfClusters, int numprocs, MPI_Status* statusP, MPI_Datatype mpiNewClusterCenter) {
	//master recives from slaves their new cluster centers and calculates the new cluster centers
	newClusterCenter* temp = NULL;
	makenewClusterCenterArray(&temp, numberOfClusters);
	for (int counter = 1; counter < numprocs; counter++)
	{
		MPI_Recv(temp, numberOfClusters, mpiNewClusterCenter, counter, 0, MPI_COMM_WORLD, statusP);
		addClusterData(mastersNewClusters, temp, numberOfClusters);
	}

#pragma omp parallel for
	for (int clusterCounter = 0; clusterCounter < numberOfClusters; clusterCounter++) {

		//setting the new center coords
		clustersArray[clusterCounter].myPoint.x = mastersNewClusters[clusterCounter].accumulatedDistance.x / mastersNewClusters[clusterCounter].numberOfPointsInCluster;
		clustersArray[clusterCounter].myPoint.y = mastersNewClusters[clusterCounter].accumulatedDistance.y / mastersNewClusters[clusterCounter].numberOfPointsInCluster;
	}
	free(temp);
	
}

void shiftArray(clusteredPoint* points, int offset, int numberOfPoints)
{
	for (int counter = 0; counter < numberOfPoints; counter++)
		points[offset + counter] = points[counter];
}
void copyMastersPoints(clusteredPoint* temp, clusteredPoint* mester, int numberOfPoints)
{
	for (int counter = 0; counter < numberOfPoints; counter++)
		temp[counter] = mester[counter];
}
void calculateNewClusterArray(newClusterCenter* mastersNewClusters, cluster* clustersArray, int numberOfClusters, clusteredPoint* points, double* distance, int numberOfPoints)
{
	setDiameter(clustersArray, numberOfClusters, points, distance, numberOfPoints);
	//O(NK) 
#pragma omp parallel for
	for (int clusterCounter = 0; clusterCounter < numberOfClusters; clusterCounter++) {

		//setting the new center coords
		clustersArray[clusterCounter].myPoint.x = mastersNewClusters[clusterCounter].accumulatedDistance.x / mastersNewClusters[clusterCounter].numberOfPointsInCluster;
		clustersArray[clusterCounter].myPoint.y = mastersNewClusters[clusterCounter].accumulatedDistance.y / mastersNewClusters[clusterCounter].numberOfPointsInCluster;
	}

	
}

void setDiameter(cluster* clustersArray, int numberOfClusters , clusteredPoint* points, double* distance, int numberOfPoints )
{
	for (int counter = 0; counter < numberOfPoints; counter++)
	{

		point2D clusterCoords = points[counter].myCentroid;
		double newDiamtre = sqrt(distance[counter]);
		for (int clusterCounter = 0; clusterCounter < numberOfClusters; clusterCounter++)
		{
			if (clusterCoords.x == clustersArray[clusterCounter].myPoint.x && clusterCoords.y == clustersArray[clusterCounter].myPoint.y)
				clustersArray[clusterCounter].diamtre = newDiamtre > clustersArray[clusterCounter].diamtre ? newDiamtre : clustersArray[clusterCounter].diamtre;
		}

	}
}
void addClusterData(newClusterCenter* mastersNewClusters, newClusterCenter* temp, int numberOfClusters)
{
#pragma omp parallel for
	for (int clusterCounter = 0; clusterCounter < numberOfClusters; clusterCounter++)
	{
		//adding the distances
		mastersNewClusters[clusterCounter].accumulatedDistance.x += temp[clusterCounter].accumulatedDistance.x;
		mastersNewClusters[clusterCounter].accumulatedDistance.y += temp[clusterCounter].accumulatedDistance.y;
		mastersNewClusters[clusterCounter].numberOfPointsInCluster += temp[clusterCounter].numberOfPointsInCluster;

	}
}
void makenewClusterCenterArray(newClusterCenter** array, int size) {
	*array = (newClusterCenter*)malloc(size * sizeof(newClusterCenter));
	if (array == NULL) {
		/* Memory could not be allocated*/
		ERROR(ERROR_MSG_1);
	}
}
void masterCalculatesNewQualitiy(cluster* clustersArray, double *qualitiy) {
	*qualitiy = calculateClusterQualitiy(clustersArray);
}
/*main algorithm*/
void k_means(int myid, clusteredPoint* clusterPointArray, int numberOfPoints, cluster* clustersArray, int numberOfClusters, double* tempQuality, double* time, int numprocs, MPI_Status* statusP, MPI_Datatype mpiCluster, MPI_Datatype mpiNewClusterCenter, MPI_Datatype mpiClusterdPoint)
{
	int keepIterating = DONT_STOP;
	int counter = 0;
	*time = 0;
	while (keepIterating == DONT_STOP)
	{
		//at the start of each iteration the master will send the clusters to the slaves
		SendClustersToSlave(clustersArray, numberOfClusters, mpiCluster);
		int flag_did_a_point_moved = 0;
		//pairing points to clusters -- checking how many points moved between clusters -- calculating the new cluster partial data
		newClusterCenter* newClusters = k_means_iteration(clusterPointArray, numberOfPoints, clustersArray, numberOfClusters, &flag_did_a_point_moved);
		if (myid != MASTER)
		{
			//each slave sends to master how many points switched cluster
			MPI_Send(&flag_did_a_point_moved, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
		}
		else
		{
			//the master receives from slaves how many points moved
			int temp = 0;
			for (int counter = 1; counter < numprocs; counter++)
			{
				MPI_Recv(&temp, 1, MPI_INT, counter, 0, MPI_COMM_WORLD, statusP);
				flag_did_a_point_moved += temp;
			}
		}
		//letting the slaves know if points moved between clusters
		MPI_Bcast(&flag_did_a_point_moved, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
		if (myid != MASTER)
		{
			slavesSend_newClusterCenter(newClusters, numberOfClusters, mpiNewClusterCenter);
			//if points didn't move between clusters than the slave needs to send its points to master
			//because the master needs to evaluates the clusters diameter
			if (flag_did_a_point_moved == 0 || counter + 1 == LIMIT)
				slavesSend_Points(clusterPointArray, numberOfPoints, mpiClusterdPoint);
		}
		else
		{
			if (flag_did_a_point_moved == 0 || counter + 1 == LIMIT)
			{	//point stopped moving between clusters we need to create the clusters with diameter
				masterCreateNewClusterCentersWithDiameter(newClusters, clusterPointArray, clustersArray, numberOfClusters, numprocs, numberOfPoints, N / numprocs, statusP, mpiNewClusterCenter, mpiClusterdPoint);
				masterCalculatesNewQualitiy(clustersArray, tempQuality);
				//resetting  the counter and moving time
				counter = 0;
				*time += dT;
				checkTermination(*tempQuality, counter, *time, &keepIterating);
			}
			else
			{	//points are still moving between clusters there is no need to set a diameter for the clusters
				masterCreateNewClusterCentersNoDiameter(newClusters, clustersArray, numberOfClusters, numprocs, statusP, mpiNewClusterCenter);
				counter++;
				keepIterating = DONT_STOP ;
			}
			printf("counter -----> %d points moved ---> %d quality --->%f \n", counter, flag_did_a_point_moved, *tempQuality);
			fflush(stdout);
		}
		broadCastSlavesLoopVar(&keepIterating, &counter, time);
		if (keepIterating == DONT_STOP && flag_did_a_point_moved == 0)
			movePoints(clusterPointArray, numberOfPoints, dT);

		free(newClusters);
	}

}
/*moving points in time using cuda*/
void movePoints(clusteredPoint* clusterPointArray , int numberOfPoints , double time)
{
	cudaError_t cudaStatus = movePointInTime(clusterPointArray, numberOfPoints, time);
	if (cudaStatus != cudaSuccess)
		ERROR(MY_CUDA_ERROR_MSG_1);

	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess)
		ERROR(MY_CUDA_ERROR_DeviceReset);
}

void checkTermination(double tempQualitiy, int counter, double time, int* keepIteratingP)
{
	*keepIteratingP = (tempQualitiy <= QM || counter >= LIMIT || time >= T )? STOP : DONT_STOP;
}

void broadCastSlavesLoopVar(int* keepIterating, int* updatedCounter, double* updatedTime)
{
	//sending the loop varibales to the slaves this way we keep the loop in sync
	MPI_Bcast(keepIterating, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(updatedCounter, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(updatedTime, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

}


void writeToFile(char* outputFilePath, cluster* clustersArray, int numberOfClusters, double* tempQualitiy, double time , double startTime)
{  //writing the ouputFile
	FILE *file = NULL;
	fopen_s(&file ,outputFilePath, "w");
	if (file == NULL)
		ERROR(ERROR_MSG_2);
	char *outmsg = "Centers of the clusters: \n";
	fprintf(file, "Ruslan Abramov - 306847393 \n");
	fprintf(file, "Time it took- %lf seconds \n", MPI_Wtime() - startTime);
	fprintf(file, "First occurrence at t = %f with q = %f \n", time, *tempQualitiy);
	fprintf(file, "%s", outmsg);
	writeClusterCenters(file, clustersArray, numberOfClusters);
	fclose(file);
}

void writeClusterCenters(FILE* file, cluster* clustersArray, int numberOfClusters) {
	//writing the cluster centers
	for (int counter = 0; counter < numberOfClusters; counter++)
		fprintf(file, "%f\t%f\n", clustersArray[counter].myPoint.x, clustersArray[counter].myPoint.y);
}

void createPoint2DMpiDataType(MPI_Datatype *point2DMPIType)
{

	point2D recivedPoint = { 0, 0 };
	MPI_Aint disp[2];
	int blocklen[2] = { 1, 1 };
	MPI_Datatype type[2] = { MPI_DOUBLE, MPI_DOUBLE };
	disp[0] = (char *)&recivedPoint.x - (char *)&recivedPoint;
	disp[1] = (char *)&recivedPoint.y - (char *)&recivedPoint;


	MPI_Type_create_struct(2, blocklen, disp, type, point2DMPIType);
	MPI_Type_commit(point2DMPIType);
}

void createVelocity2DMpiDataType(MPI_Datatype *velocity2DMPIType)
{
	MPI_Datatype type[2] = { MPI_DOUBLE, MPI_DOUBLE };
	int blocklen[2] = { 1, 1 };
	MPI_Aint disp[2] = { 0, sizeof(double) };
	MPI_Type_create_struct(2, blocklen, disp, type, velocity2DMPIType);
	MPI_Type_commit(velocity2DMPIType);
}


void createMovingPoint2DMpiDataType(MPI_Datatype *movingPoint2DMPIType, MPI_Datatype point2DMPIType, MPI_Datatype velocity2DMPIType)
{

	MPI_Datatype type[2] = { point2DMPIType, velocity2DMPIType };
	int blocklen[2] = { 1, 1 };
	MPI_Aint disp[2] = { 0, sizeof(point2D) };
	MPI_Type_create_struct(2, blocklen, disp, type, movingPoint2DMPIType);
	MPI_Type_commit(movingPoint2DMPIType);

}

void createClusteredPointMpiDataType(MPI_Datatype *clusteredPoints, MPI_Datatype movingPoint2DMPIType, MPI_Datatype point2DMPIType)
{

	MPI_Datatype type[3] = { movingPoint2DMPIType, point2DMPIType, MPI_INT };
	int blocklen[3] = { 1, 1, 1 };
	MPI_Aint disp[3] = { 0, sizeof(movingPoint2D) , sizeof(movingPoint2D)+ sizeof(point2D) };
	MPI_Type_create_struct(3, blocklen, disp, type, clusteredPoints);
	MPI_Type_commit(clusteredPoints);

}
void createNewClusterCenterMpiDataType(MPI_Datatype *newClusterCenter, MPI_Datatype point2DMPIType)
{
	MPI_Datatype type[4] = { point2DMPIType, point2DMPIType, MPI_DOUBLE, MPI_INT };
	int blocklen[4] = { 1, 1, 1, 1 };
	MPI_Aint disp[4] = { 0, sizeof(point2D), 2 * sizeof(point2D), 2 * sizeof(point2D) + sizeof(double) };
	MPI_Type_create_struct(4, blocklen, disp, type, newClusterCenter);
	MPI_Type_commit(newClusterCenter);
}

void createClusterMpiDataType(MPI_Datatype *cluster, MPI_Datatype point2DMPIType)
{
	MPI_Datatype type[2] = { point2DMPIType, MPI_DOUBLE };
	int blocklen[2] = { 1, 1 };
	MPI_Aint disp[2] = { 0, sizeof(point2D) };
	MPI_Type_create_struct(2, blocklen, disp, type, cluster);
	MPI_Type_commit(cluster);
}

void finalize(clusteredPoint* clusterPointArray, cluster* clustersArray, MPI_Datatype *point2DMPIType, MPI_Datatype *velocity2DMPIType, MPI_Datatype *movingPoint2DMPIType, MPI_Datatype *clusteredPoints, MPI_Datatype *newClusterCenter, MPI_Datatype *cluster)
{
	MPI_Type_free(point2DMPIType);
	MPI_Type_free(velocity2DMPIType);
	MPI_Type_free(movingPoint2DMPIType);
	MPI_Type_free(clusteredPoints);
	MPI_Type_free(newClusterCenter);
	MPI_Type_free(cluster);
	free(clusterPointArray);
	free(clustersArray);
	MPI_Finalize();
}
