Ruslan Abramov 306847393
Ravid Anbary 203373360

Some notes:
1. our array size is 8Mb we did it this size so that it will divide equally into thread blocks in cuda
2. we used 8Kb of cuda blocks each block contains 512 threads
3. avoiding race conditions in cuda was made by atomic add function
4. avoiding race conditions in openMP was made by doing some sort of reduction,
   each thread has it's own RGB size array,
   after the calculation, each thread copies its array into the shared RGB size array inside of a critical section,
   in this way, we managed to avoid blocking the 4Mb random numbers array and managed to block only a small array.
5. The last thing to note is the 2 Intellisense error in the cuda code, 
   we are using VS 2013 it's a known issue that calling kernel function causes Intellisense error,
   and for some reason, the atomic add function inside the kernel creates Intellisense error. Never the last the code compiles without any warnings or errors.
6. Thank you! 