#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define FILEPATH "matrix.txt"
/********global consts ************/
const int NOT_A_SQAURE_NUMBER	= - 1;
const int UP					= 0;
const int DOWN					= 1;
const int LEFT					= 2;
const int RIGHT					= 3;
const int NO_ERROR				= 0;
const int THERE_IS_AN_ERROR		= 1;
const int MASTER				= 0;
const int NOT_SORTED			= 0;
const int SORTED				= 4;
const int EVEN					= 0;
const int ODD					= 1;
const int BIGGER				= 1;
const int EQUAL					= -1;
const int EDGE					= 1;

/* point struct */
typedef struct {
	double x; // x coor
	double y; // y coor
	double distanceFromP1; // distance from {10 ,20}
}point;
const point P1 = { 10, 20 , 0 };

/*******Error Msg********/
const char* ERROR_MSG_1 = "\n Memory could not be allocated --> aborting";
const char* ERROR_MSG_2 = "\n Matrix size must be a square number --> aborting";
const char* ERROR_MSG_3 = "\n Number of Process must be equal to matrix size --> aborting";

/********Point related************/ 
double distanceFromP(point p1, double x, double y);
void createPointMpiDataType(MPI_Datatype *PointMPIType);
/********Matrix related*************/
point** readPointsFromFile(char* filePath, int* matrixSize);
point** makePointMat(int size);
void printMat(point** pointMatrix , int size , int order);
void freeMat(point** pointMatrix);
int isColEdge(int col, int size);
int isRowEdge(int row, int size);
/********Sort related*********/
int isMyDistanceBigger(point myPoint, point myNeighbor);
void oddEvenSort(int myid, int numprocs, int *nbrs, int * coords,
					point *myPoint, MPI_Datatype PointMPIType, MPI_Comm cartcomm, MPI_Status *status,
																						int colRowSize);
void innerSwitch(point *myPoint, int *evenOrOddStep, int *nbrs,
					int *didIChanged, int oddOrEvenPos, int rowsOrCols,
						MPI_Datatype PointMPIType, MPI_Status *status, MPI_Comm cartcomm,
															int isOnEdge, int isMyRowOddOrEven);
void sortEvenStepRowsOrCols(point *myPoint, int* nbrs,
								int* didIChanged, int oddOrEvenPos, int rowsOrColls,
									MPI_Comm cartcomm, MPI_Datatype PointMPIType, MPI_Status *status,
																					int isMyRowOddOrEven );
void sortOddStepRowsOrCols(point *myPoint, int* nbrs,
							int* didIChanged, int oddOrEvenPos, int rowsOrColls,
								MPI_Comm cartcomm, MPI_Datatype PointMPIType, MPI_Status *status,
																	int isOnEdge, int isMyRowOddOrEven);
void sortOddRow(int isBigger, int isEdge, int oddOrEvenPos, int oddOrEvenStep, int *toChange);
void chekingOddStep(int *didIChanged, int isOnEdge, int isMyDisBig, int oddOrEvenPos);
void isSortedMatrix(int *isSorted, int *didIChanged, int myid, int numprocs, MPI_Comm cartcomm, MPI_Status *status);
int myNbr(int oddOrEvenStep, int oddOrEvenPos, int *nbrs, int rowsOrColls);

/*******file and init related*********/
void masterReadAndPrint(point ***matrixP, int *error, char *filePath, int *size);
void masterFinal(point ***matrixP, int *size);
void sendRecivePoints(point** matrix, MPI_Datatype PointMPIType, MPI_Comm comm, MPI_Status *status, point *myPoint, int numprocs, int myid, int size);
/*********sending the sorted points to the master*********/
void recivePointsFromSlaves(point*** matrix, MPI_Datatype PointMPIType, MPI_Comm comm, MPI_Status *status, point *myPoint, int numprocs, int myid, int size);

/***********MPI Cartesian tooplogy*********/
void cartesianToop(int numOfDims, int *dimSize, int *isCyclicDim, int reOrder, MPI_Comm *cartcomm, int myid, int *myCoords, int *myNbrs);

int main(int argc, char *argv[]){
	char* filePath = FILEPATH;
	int size, error = 0;
	point **matrix = NULL, myPoint = { 0, 0, 0 };
	int   myid, numprocs;
	MPI_Comm cartcomm;
	MPI_Status status;
	//Mpi init for point struct
	MPI_Datatype PointMPIType;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	// Create MPI user data type for Point
	createPointMpiDataType(&PointMPIType);

	//master reads the data
	if (myid == MASTER)
		masterReadAndPrint(&matrix, &error, filePath, &size);
	//sending to all process the size of matrix and if there is an error
	MPI_Bcast(&size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&error, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	if (numprocs != size * size) {
		error = THERE_IS_AN_ERROR;
		if (myid == MASTER) {
			printf("%s", ERROR_MSG_3);
			printf("\n%d", size);
		}
	}
	
	if (error == NO_ERROR)
	{
		sendRecivePoints(matrix, PointMPIType, MPI_COMM_WORLD, &status, &myPoint, numprocs, myid, size);
		int dims[2] = { size, size }, isCyclic[2] = { 1, 1 };
		int coords[2], nbrs[4];
		//create the cartesian virtual topology, get rank, coordinates, neighbor ranks
		cartesianToop(2, dims, isCyclic, 0, &cartcomm, myid, coords, nbrs);

		oddEvenSort(myid, numprocs, nbrs, coords, &myPoint, PointMPIType, cartcomm, &status , size);
		
		//send my new point and coords to the master
		if (myid != MASTER){
			MPI_Send(&myPoint, 1, PointMPIType, MASTER, 0, cartcomm);
			MPI_Send(coords, 2, MPI_INT, MASTER, 0, cartcomm);
		}
		//get the sorted points
		if (myid == MASTER){
			recivePointsFromSlaves(&matrix, PointMPIType, cartcomm, &status, &myPoint, numprocs, myid, size);
			masterFinal(&matrix, &size);
		}
	}

	MPI_Type_free(&PointMPIType);
	MPI_Finalize();
}

//reading from the file
point** readPointsFromFile(char* filePath, int* matrixSize){
	int size;
	point** pointMatrix = NULL;
	FILE * file;
	fopen_s(&file,filePath, "r");
	if (file){

		//first we read n^2 which is the amount of points to read
		//and therefore the size of our matrix
		fscanf_s((FILE*)file, "%d", &size);
		//normalizing size and checking if size = n^2
		size = (int)sqrt((float)size)*(int)sqrt((float)size) == size ? 
			(int)sqrt((float)size) : -1;
		if (size == NOT_A_SQAURE_NUMBER || size % 2 == ODD){
			//the input is not valid
			printf("%s", ERROR_MSG_2);
			return NULL; 
		}
		//creating a point matrix size of n*n
		pointMatrix = makePointMat(size);

		//reading n*n points
		for (int i = 0, rows = 0, cols; i < size*size;){
			//points are saved in file in this way:
			//"x y\n"
			float x;
			float y;
			fscanf_s((FILE*)file, "%f", &x);
			fscanf_s((FILE*)file, "%f", &y) ;
			cols = i%size;
			pointMatrix[rows][cols] = //placing a new point
				point{ x, y, distanceFromP(P1, x, y) };
			i++;
			//after cols == n we go to the next row
			rows = i%size == 0 ? rows + 1 : rows;
		}
		fclose(file);
	}
	*matrixSize = size;
	return pointMatrix;
}

point** makePointMat(int size){
	point** pointMatrix = (point **)malloc(size * sizeof(point*));
	//pointer to an array of [lines] pointers
	if (pointMatrix == NULL) {
		/* Memory could not be allocated*/
		printf("%s", ERROR_MSG_1);
		return NULL;
	}

	for (int i = 0; i < size; ++i){
		pointMatrix[i] = (point *)malloc(size * sizeof(point));
		if (pointMatrix[i] == NULL) {
			/* Memory could not be allocated*/
			printf("%s", ERROR_MSG_1);
			return NULL;
		}
	}
	return pointMatrix;
}
double distanceFromP(point p1, double x, double y)
{	//the distance between points : ((x1-x2)^2 + (y1 - y2)^2)^0.5
	return sqrt(pow(p1.x - x, 2) + pow(p1.y - y, 2));
}

//printing the matrix
void printMat(point** pointMatrix , int size, int order){
	for (int rows = 0; rows <  size; rows++){
		for (int cols = 0; cols < size; cols++){
			int colPos = order == 1 && rows % 2 == 1 ? size - 1 - cols : cols;
			printf("(%.3f, %.3f --> %.3f)\n", 
				pointMatrix[rows][colPos].x, pointMatrix[rows][colPos].y, pointMatrix[rows][colPos].distanceFromP1);
		}

	}
}

//releasing memory
void freeMat(point** pointMatrix){
	for (int rows = 0; rows < sizeof pointMatrix / sizeof pointMatrix[0]; rows++)
		free(pointMatrix[rows]);
	free(pointMatrix);
}

//return 1 if mydistance is bigger then my neighbor
int isMyDistanceBigger(point myPoint, point myNeighbor){
	if (myPoint.distanceFromP1 > myNeighbor.distanceFromP1)
		return BIGGER;
	else
		return myPoint.distanceFromP1 == myNeighbor.distanceFromP1 ? EQUAL : 0;
}


void masterReadAndPrint(point ***matrixP, int *error, char *filePath, int *size){
	*matrixP = readPointsFromFile(filePath, size);
	*error = *matrixP == NULL ? THERE_IS_AN_ERROR : NO_ERROR;
	if (*error == NO_ERROR){
		printMat(*matrixP, *size , 0);
	}
}

void masterFinal(point ***matrixP , int *size){
	printf("The sorted points : \n");
	printMat(*matrixP, *size , 1);
	freeMat(*matrixP);
}

void sendRecivePoints(point** matrix, MPI_Datatype PointMPIType, MPI_Comm comm, MPI_Status *status, point *myPoint, int numprocs, int myid , int size){
	if (numprocs == size * size) { // checking if number of process is equal to the matrix size
		if (myid == MASTER) 
		{ //the master sends each point to its process
			for (int counter = 1, row = 0; counter < numprocs;) {
				MPI_Send(&matrix[row][counter % size], 1, PointMPIType, counter, 0, comm);
				counter++;
				row = counter % size == 0 ? row + 1 : row;
			}
			*myPoint = matrix[0][0];
		}
		else 
		{	//each slave receives a point from the master
			MPI_Recv(myPoint, 1, PointMPIType, MASTER, 0, comm, status);
		}
	}
}

void recivePointsFromSlaves(point*** matrix, MPI_Datatype PointMPIType, MPI_Comm comm, MPI_Status *status,
	point *myPoint, int numprocs, int myid, int size)
{
	point tempP;

	int coords[2];
	(*matrix)[0][0] = *myPoint;
	for (int counter = 1; counter < numprocs; counter++){
		MPI_Recv(&tempP, 1, PointMPIType, MPI_ANY_SOURCE, 0, comm, status);
		MPI_Recv(&coords, 2, MPI_INT, (*status).MPI_SOURCE, 0, comm, status);
		(*matrix)[coords[0]][coords[1]] = tempP;
	}
}

void createPointMpiDataType(MPI_Datatype *PointMPIType){
	MPI_Datatype type[3] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
	int blocklen[3] = { 1, 1, 1 };
	MPI_Aint disp[3] = { 0, sizeof(double), sizeof(double) * 2 };
	MPI_Type_create_struct(3, blocklen, disp, type, PointMPIType);
	MPI_Type_commit(PointMPIType);
}


void oddEvenSort(int myid, int numprocs, int *nbrs, int * coords, point *myPoint,
					MPI_Datatype PointMPIType, MPI_Comm cartcomm, MPI_Status *status, int colRowSize){
	int isSorted = NOT_SORTED;
	int iteration = 0, didIChanged = 0;
	int isMyRowOddOrEven = coords[0] % 2 == 0 ? EVEN : ODD, isMyColOddOrEven = coords[1] % 2 == 0 ? EVEN : ODD;
	//cheking if coords on the edge for odd steps
	int isColOnEdge = isColEdge(coords[1], colRowSize), isRowOnEdge = isRowEdge(coords[0], colRowSize);
	//even or odd step
	int workingOnRows = EVEN, workingOnCols = EVEN;
	//sorting
	while (isSorted != SORTED )
	{
		if (iteration % 2 == 0)		//working on rows switching cols
			innerSwitch(myPoint, &workingOnRows, nbrs, &didIChanged, isMyColOddOrEven,
						(iteration + 1) % 2, PointMPIType, status, cartcomm, isColOnEdge, isMyRowOddOrEven);

		else		 //working on col switching rows
			innerSwitch(myPoint, &workingOnCols, nbrs, &didIChanged, isMyRowOddOrEven,
						(iteration + 1) % 2, PointMPIType, status, cartcomm, isRowOnEdge, isMyRowOddOrEven);
		
		iteration++;
		isSortedMatrix(&isSorted, &didIChanged, myid, numprocs, cartcomm, status);
	}
	if (myid == MASTER){
		printf("\n****************************\nIt took %d iterations to sort the matrix\n", iteration);
		printf("Which is log%d^2 + c = %f + %f \n", colRowSize, 2 * (log((double)colRowSize) / log(2.0)), iteration - 2 * log((double)colRowSize) / log(2.0));
		printf("In this case the complaxity needs to be O(logn) + 1 where n is the size of matrix\n");
	}
}




void isSortedMatrix(int *isSorted , int *didIChanged ,int myid, int numprocs , MPI_Comm cartcomm , MPI_Status *status)
{
	//telling the master if we made a change
	if (myid != MASTER){
		MPI_Send(didIChanged, 1, MPI_INT, MASTER, 0, cartcomm);
	}
	//at the end of each iteration we need to check if the matrix is sorted
	else{
		int tempDidSomebodyChanged;
		for (int counter = 1; counter < numprocs; counter++){
			MPI_Recv(&tempDidSomebodyChanged, 1, MPI_INT, MPI_ANY_SOURCE, 0, cartcomm, status);
			*didIChanged += tempDidSomebodyChanged;
			}
		//we need that all 4 sort oporations will have no changes
		*isSorted = *didIChanged == 0 ? *isSorted + 1 : NOT_SORTED;
		}
	MPI_Bcast(isSorted, 1, MPI_INT, MASTER, cartcomm);
}


void innerSwitch(point *myPoint, int *evenOrOddStep, int *nbrs,
					int *didIChanged, int oddOrEvenPos, int rowsOrCols,
						MPI_Datatype PointMPIType, MPI_Status *status, MPI_Comm cartcomm,
														int isOnEdge, int isMyRowOddOrEven)
{	//evenStep
	if (*evenOrOddStep == EVEN){
		sortEvenStepRowsOrCols(myPoint, nbrs, didIChanged, oddOrEvenPos, rowsOrCols, cartcomm, PointMPIType, status, isMyRowOddOrEven);
		*evenOrOddStep = ODD;
	}//oddStep
	else{
		sortOddStepRowsOrCols(myPoint, nbrs, didIChanged, oddOrEvenPos, rowsOrCols, cartcomm, PointMPIType, status, isOnEdge, isMyRowOddOrEven);
		*evenOrOddStep = EVEN;
	}
}

void sortEvenStepRowsOrCols(point *myPoint, int* nbrs, int* didIChanged,
								int oddOrEvenPos, int rowsOrColls, MPI_Comm cartcomm,
									MPI_Datatype PointMPIType, MPI_Status *status, int isMyRowOddOrEven)
{
	point neighboursPoint;
	int myNbrs = myNbr(EVEN, oddOrEvenPos, nbrs, rowsOrColls);
	MPI_Send(myPoint, 1, PointMPIType, myNbrs, 0, cartcomm);
	MPI_Recv(&neighboursPoint, 1, PointMPIType, myNbrs, 0, cartcomm, status);
	int isMyDisBig = isMyDistanceBigger(*myPoint, neighboursPoint);
	//sorting cols during an EVEN step when the col is in an odd row
	if (isMyDisBig == EQUAL)
		*didIChanged = 0;
	else if (isMyRowOddOrEven == ODD && rowsOrColls == 1)
		sortOddRow(isMyDisBig, 0, oddOrEvenPos, EVEN, didIChanged);
	/*
	*sorting cols during an even step when the col is on an even row
	*sorting rows during an even step
	*/
	else
	{
		if (oddOrEvenPos == EVEN && isMyDisBig == BIGGER || oddOrEvenPos == ODD && isMyDisBig != BIGGER)
			*didIChanged = 1;
		else
			*didIChanged = 0;
	}
	*myPoint = *didIChanged == 1 ? neighboursPoint : *myPoint;
}
void sortOddStepRowsOrCols(point *myPoint, int* nbrs, int* didIChanged,
								int oddOrEvenPos, int rowsOrColls, MPI_Comm cartcomm,
									MPI_Datatype PointMPIType, MPI_Status *status, int isOnEdge, int isMyRowOddOrEven )
{
	point neighboursPoint;
	int myNbrs = myNbr(ODD, oddOrEvenPos, nbrs, rowsOrColls);
	MPI_Send(myPoint, 1, PointMPIType, myNbrs,0, cartcomm);
	MPI_Recv(&neighboursPoint, 1, PointMPIType, myNbrs,0, cartcomm, status);
	int isMyDisBig = isMyDistanceBigger(*myPoint, neighboursPoint);
	if (isMyDisBig == EQUAL)
		*didIChanged = 0;
	//sorting cols during an odd step when the col is in an odd row
	else if (isMyRowOddOrEven == ODD && rowsOrColls == 1)
		sortOddRow(isMyDisBig, isOnEdge, oddOrEvenPos, ODD, didIChanged);
	/*
	*sorting cols during an odd step when the col is on an even row
	*sorting rows during an odd step
	*/
	else  
		chekingOddStep(didIChanged, isOnEdge, isMyDisBig, oddOrEvenPos);

	*myPoint = *didIChanged == 1 ? neighboursPoint : *myPoint;
	
}

int isColEdge(int col, int size)
{
	return col % size == 0 || col % size == size - 1 ? EDGE : 0;
}

int isRowEdge(int row, int size)
{
	return row  == 0 || row == size - 1 ? EDGE : 0;
}

void sortOddRow(int isBigger, int isEdge, int oddOrEvenPos ,int oddOrEvenStep, int *toChange) {
	if (oddOrEvenStep == EVEN)
	{
		if (oddOrEvenPos == EVEN && isBigger != BIGGER || oddOrEvenPos == ODD && isBigger == BIGGER)
			*toChange = 1;
		else
			*toChange = 0;
	}
	else{
		if (isEdge == EDGE  && (oddOrEvenPos == EVEN && isBigger != BIGGER || oddOrEvenPos == ODD && isBigger == BIGGER))
			*toChange = 1;
		else if (isEdge != EDGE && (oddOrEvenPos == EVEN && isBigger == BIGGER || oddOrEvenPos == ODD && isBigger != BIGGER))
			*toChange = 1;
		else
			*toChange = 0;
	}
	
}

void chekingOddStep(int *didIChanged, int isOnEdge, int isMyDisBig, int oddOrEvenPos)
{
	if (isOnEdge == EDGE && (oddOrEvenPos == ODD && isMyDisBig != BIGGER || oddOrEvenPos == EVEN && isMyDisBig == BIGGER))
		*didIChanged = 1;
	else if (isOnEdge != EDGE && (oddOrEvenPos == ODD && isMyDisBig == BIGGER || oddOrEvenPos == EVEN && isMyDisBig != BIGGER))
		*didIChanged = 1;
	else
		*didIChanged = 0;
}


int myNbr(int oddOrEvenStep ,int oddOrEvenPos, int *nbrs, int rowsOrColls){
	if (oddOrEvenStep == EVEN)
		return oddOrEvenPos == EVEN ? nbrs[rowsOrColls == 1 ? RIGHT : DOWN] : nbrs[rowsOrColls == 1 ? LEFT : UP];
	else
		return oddOrEvenPos == ODD ? nbrs[rowsOrColls == 1 ? RIGHT : DOWN] : nbrs[rowsOrColls == 1 ? LEFT : UP];
}

void cartesianToop(int numOfDims, int *dimSize, int *isCyclicDim, int reOrder, MPI_Comm *cartcomm,int myid, int *myCoords, int *myNbrs)
{
	MPI_Cart_create(MPI_COMM_WORLD, numOfDims, dimSize, isCyclicDim, 0, cartcomm);
	MPI_Comm_rank(*cartcomm, &myid);
	MPI_Cart_coords(*cartcomm, myid, 2, myCoords);
	MPI_Cart_shift(*cartcomm, 0, 1, &myNbrs[UP], &myNbrs[DOWN]);
	MPI_Cart_shift(*cartcomm, 1, 1, &myNbrs[LEFT], &myNbrs[RIGHT]);
}