This is our share sort implementation

Assumptions:
1. The matrix size in the file must be an even square number (4 , 16 , 36 , 64 , 100 , etc)
2. The matrix has the same amount of rows and cols
3. The first number in the file is the matrix size (n^2)
4. The next n^2 items in the file are tuples representing a point in the Euclid space
5. The number of processes must be n^2
6. The points are sorted by their distance from the point (10 , 20)
7. The file we read needs to be located in mpich2\bin where the mpiexec.exe is located
	because we are using mpiexec.exe his folder is our current working directory
8. The file name is matrix.txt