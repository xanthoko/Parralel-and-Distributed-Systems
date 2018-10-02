#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "pagerank.h"
#define max_e 0.00001

int N, max_it=1, k, VERBOSE = 0, threads;
int *ingoing;
float *b, *page_rank;

int main(int argc, char **argv){

	if (argc < 3){
		printf("something bad happened.\n");
		exit(1);
	}

	N = atoi(argv[1]);  // number of websites
	max_it = atoi(argv[2]);  // number of maximum iterations
	float const_term = (1-0.85)/N;  // constant that will fill the b matrix
	threads = atoi(argv[3]);  // number of threads in case of parallel executing

	// --------------allocate memory for required arrays--------------
	b = (float *)malloc(N * sizeof(float));
	page_rank = (float *)malloc(N * sizeof(float));
	ingoing = (int *)malloc(N * sizeof(int)); 
	float **adj_matrix = (float **)malloc(N * sizeof(float *));

	arrays_init(const_term);

	// read ingoing matrix
	FILE * fp;
	fp = fopen ("data/ingoing.txt", "r");

	for (int i=0; i<N; i++){ 
		int k = fscanf(fp, "%d", &ingoing[i]);
		adj_matrix[i] = malloc(2 * ingoing[i] * sizeof(float));
	}
	fclose(fp);

	// read adjacency matrix
	fp = fopen ("data/sparse.txt", "r");

	for (int i=0; i<N; i++){ 
		for (int j=0; j<2*ingoing[i]; j++){
		  int k = fscanf(fp, "%f", &adj_matrix[i][j]);
		}
	}
	fclose(fp);

	// perform the gauss seidel method to solve the pagrank problem
	// clock_t is used to count the time taken to complete the solving
	clock_t time;
	time = clock(); 
	gs(adj_matrix);
	time = clock() - time; 
	double time_taken = ((double)time)/CLOCKS_PER_SEC; // in seconds
	printf("TIME TAKEN: %f\n", time_taken);	

	return 0;
}


// initiliaze arrays. Initial values of pagerank cells are set to 1/(number of websites)
void arrays_init(float b_val){
	for (int i=0;i<N;i++){
		b[i] = b_val;
		page_rank[i] = 1.0/N;
	}
}


// The Gauss Seidel algorithm implementation.
void gs(float** array){
	int i, row, col, j, it;
	float row_sum, error, temp, max_error;

	for (it=0;it<max_it;it++){
		#pragma omp parallel for num_threads(threads)
		for(i=0;i<N;i++){
			row_sum = 0;
			row = i;
			for(j=0;j<ingoing[i];j++){
				col = array[i][2*j];
				row_sum += array[row][2*j+1] * page_rank[col];
			}
			temp = (b[i] - row_sum);
			error = fabs(page_rank[i] - temp);
			if (error>max_e){
				max_error = error;
			} 
			page_rank[i] = temp;
		}
		// printf("iteration: %d error: %f\n", it+1, max_error);
	}
}
