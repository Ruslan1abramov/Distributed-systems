//This header contatins the structs being used in the k - means algorithm

/* 2D point */
typedef struct {
	double x; // x coords
	double y; // y coords
} point2D;

/* velocity in 2D */
typedef struct {
	double xVelocity; // x velocity
	double yVelocity; // y velocity
} velocity2D;

/* moving point struct */
typedef struct {
	point2D initCoords; //starting point
	velocity2D velocity;
} movingPoint2D;

/* moving point assigned to a cluster*/
typedef struct {
	movingPoint2D point;
	point2D myCentroid;
	int myClusterId;
} clusteredPoint;

/* the slave process needs to calculate a part of the new clusters center*/
typedef struct{
	point2D myPoint;
	point2D accumulatedDistance; //each point assigned to a cluster will be add
	double diamtre;
	int numberOfPointsInCluster;
	//int id;
	/* when the master process will receive this struct then it can easly determine the new clusters center*/
} newClusterCenter;

typedef struct{
	point2D myPoint;
	double diamtre;
	//int id;
} cluster;