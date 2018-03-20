int	N = 0; // number of points
int	K = 0; // number of clusters to find
int	LIMIT = 0;// the maximum number of iterations for K - MEAN algorithm.
double QM = 0;// quality measure to stop
double T = 0;// defines the end of time interval[0, T]
double dT = 0;// defines moments t = n*dT, n = { 0, 1, 2, …, T / dT } for which calculate the clusters and the quality

const int MASTER = 0;
const int STOP = 0;
const int DONT_STOP = 1;

const char* OUTPUTFILENAME = "distrbuted_k_means_output_ruslan_abramov_306847393.txt";
/*******Error Msg********/
const char* ERROR_MSG_1 = "\n Memory could not be allocated --> aborting";
const char* ERROR_MSG_2 = "Error opening file!\n";
const char* MY_CUDA_ERROR_MSG_1 = "computeFirstClustersWithCuda failed!";
const char* MY_CUDA_ERROR_DeviceReset = "cudaDeviceReset failed!";