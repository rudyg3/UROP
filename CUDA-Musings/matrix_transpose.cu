#include <stdio.h>
#include "../../common/book.h"

const int ROWS = 3000;
const int COLS = 2000;


// Arguments are pointers to 
__global__ void transpose(int (*a)[COLS], int (*t)[ROWS]) {
  int cur_row = blockIdx.y;
  int cur_col = blockIdx.x;
  t[cur_row][cur_col] = a[cur_col][cur_row];
}


int main( void ) {
 
  // Declare arrays to be used on host (allocate memory off heap to prevent segfault)
  int *A = (int *)malloc(ROWS * COLS * sizeof(int));
  int *T = (int *)malloc(ROWS * COLS * sizeof(int));

  // Declare pointers to be used/evaluated on device
  // Pointers must be of type array containing number of elements in first element of 2D array
  int (*Aptr)[COLS];
  int (*Tptr)[ROWS];

  dim3 grid(ROWS,COLS);

  // Allocate space for the pointers on the device
  HANDLE_ERROR( cudaMalloc( (void**)&Aptr, ROWS * COLS * sizeof(int) ) ); 
  HANDLE_ERROR( cudaMalloc( (void**)&Tptr, ROWS * COLS * sizeof(int) ) );

  // Fill up matrix A
  printf("Filling A matrix...\n");
  int fill = 0;
  for(int j = 0; j < ROWS; j++) {
    for(int i = 0; i < COLS; i++) {
      A[i+j*COLS] = fill;
      fill++;
      //printf("A[%d][%d] = %d\n",j,i,A[i+j*COLS]);
	}
  }

  // Point data stored in A on host to device with memory address: Aptr
  HANDLE_ERROR( cudaMemcpy( Aptr, A, ROWS * COLS * sizeof(int), cudaMemcpyHostToDevice ) );

  printf("Transposing matrix...\n");
  transpose<<<grid,1>>>(Aptr, Tptr); // Two-dimensional grid of blocks simulates matrix

  // Copy back transposed matrix (Tptr) to host
  HANDLE_ERROR( cudaMemcpy( T, Tptr, ROWS * COLS * sizeof(int), cudaMemcpyDeviceToHost ) );

  // Check to see if transposition worked
  int success = 1;
  for(int j = 0; j < COLS; j++) {
    for(int i = 0; i < ROWS; i++) {      
      //printf("T[%d][%d] = %d\n",j,i,T[i+j*ROWS]);
      if (T[i+j*ROWS] != A[j+i*COLS])
      	 success = 0;
    }
  }

  if (success) {
     printf("Success!\n");
  }

  else {
     printf("Failure\n");
  }

 
  cudaFree( Aptr );
  cudaFree( Tptr );
  free( A );
  free( T );

  return 0;
}
