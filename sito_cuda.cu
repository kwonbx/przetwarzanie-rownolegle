#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cuda_runtime.h>

#define N 17 
#define EPS 1e-6
#define MAX_G6_LEN 32
#define BATCH_SIZE 100000

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", \
                    __FILE__, __LINE__, err, cudaGetErrorString(err), #call); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

__device__ void getEigenvalues_dev(double m[N][N], int n, double evals[N]) {
    for(int it=0; it<500; it++) {
        double max_v = 0;
        int p=0, q=1;
        for(int i=0; i<n; i++) {
            for(int j=i+1; j<n; j++) {
                if(fabs(m[i][j]) > max_v) { 
                    max_v = fabs(m[i][j]); 
                    p=i; q=j; 
                }
            }
        }
        if(max_v < 1e-12) break;

        double theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
        double t;
        if (fabs(theta) < 1e10) {
            t = (theta >= 0) ? 1.0/(theta + sqrt(1.0+theta*theta)) : -1.0/(-theta + sqrt(1.0+theta*theta));
        } else {
            t = 0.5 / theta;
        }
        
        double c = 1.0 / sqrt(1.0+t*t);
        double s = t*c;
        double tau = s/(1.0+c);
        double apq = m[p][q];
        
        m[p][p] -= t*apq; 
        m[q][q] += t*apq; 
        m[p][q] = m[q][p] = 0;
        
        for(int i=0; i<n; i++) {
            if(i!=p && i!=q) {
                double aip = m[i][p], aiq = m[i][q];
                m[i][p] = m[p][i] = aip - s*(aiq + aip*tau);
                m[i][q] = m[q][i] = aiq + s*(aip - aiq*tau);
            }
        }
    }
    for(int i=0; i<n; i++) evals[i] = m[i][i];
}

__device__ int decodeGraph6_dev(const char* g6, double adj[N][N]) {
    if (g6[0] < 63 || g6[0] > 126) return 0;
    int n = g6[0] - 63;

    for(int i=0; i<n; i++) 
        for(int j=0; j<n; j++) adj[i][j] = 0.0;

    int char_ptr = 1;
    int bit_ptr = 0;
    int current_val = 0;

    for (int j = 1; j < n; j++) {
        for (int i = 0; i < j; i++) {
            if (bit_ptr == 0) {
                current_val = g6[char_ptr++] - 63;
                bit_ptr = 6;
            }
            if ((current_val >> (bit_ptr - 1)) & 1) {
                adj[i][j] = adj[j][i] = 1.0;
            }
            bit_ptr--;
        }
    }
    return n;
}

__global__ void processGraphsKernel(const char* g6_batch, int* results, int count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < count) {
        char local_g6[MAX_G6_LEN];
        double adj[N][N];
        double evals[N];
        
        int base_idx = idx * MAX_G6_LEN;
        for(int i = 0; i < MAX_G6_LEN; i++) {
            local_g6[i] = g6_batch[base_idx + i];
            if(local_g6[i] == '\0') break;
        }

        int n = decodeGraph6_dev(local_g6, adj);
        
        int isIntegral = 0;
        if (n > 0) {
            getEigenvalues_dev(adj, n, evals);
            isIntegral = 1;
            for(int i=0; i<n; i++) {
                double rounded = round(evals[i]);
                if(fabs(evals[i] - rounded) > EPS) {
                    isIntegral = 0;
                    break;
                }
            }
        }

        results[idx] = isIntegral;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        return 1;
    }

    FILE *f = fopen(argv[1], "r");
    if (!f) return 1;

    long max_lines = -1;
    if (argc >= 3) {
        max_lines = atol(argv[2]);
        if (max_lines <= 0) {
            fclose(f);
            return 1;
        }
    }

    char *h_g6_batch = (char*)malloc(BATCH_SIZE * MAX_G6_LEN * sizeof(char));
    int  *h_results  = (int*)malloc(BATCH_SIZE * sizeof(int));

    char *d_g6_batch;
    int  *d_results;
    CUDA_CHECK(cudaMalloc((void**)&d_g6_batch, BATCH_SIZE * MAX_G6_LEN * sizeof(char)));
    CUDA_CHECK(cudaMalloc((void**)&d_results, BATCH_SIZE * sizeof(int)));

    char line[1024];
    long lines_read = 0;
    int id = 1;
    int current_batch_size = 0;

    clock_t start_time = clock();
    double start_time_sec = (double)start_time / CLOCKS_PER_SEC;

    int threadsPerBlock = 256;

    while (max_lines == -1 || lines_read < max_lines) {
        if (fscanf(f, "%1023s", line) != 1) break;
        
        lines_read++;
        
        strncpy(&h_g6_batch[current_batch_size * MAX_G6_LEN], line, MAX_G6_LEN - 1);
        h_g6_batch[current_batch_size * MAX_G6_LEN + MAX_G6_LEN - 1] = '\0';
        
        current_batch_size++;

        if (current_batch_size == BATCH_SIZE) {
            CUDA_CHECK(cudaMemcpy(d_g6_batch, h_g6_batch, current_batch_size * MAX_G6_LEN * sizeof(char), cudaMemcpyHostToDevice));
            
            int blocksPerGrid = (current_batch_size + threadsPerBlock - 1) / threadsPerBlock;
            
            processGraphsKernel<<<blocksPerGrid, threadsPerBlock>>>(d_g6_batch, d_results, current_batch_size);
            CUDA_CHECK(cudaGetLastError());
            
            CUDA_CHECK(cudaMemcpy(h_results, d_results, current_batch_size * sizeof(int), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaDeviceSynchronize());

            clock_t current_time = clock();
            double time_found = (double)(current_time - start_time) / CLOCKS_PER_SEC;

            for (int i = 0; i < current_batch_size; i++) {
                if (h_results[i]) {
                    printf("%d: %s (czas pobrania paczki: %f s)\n", id, &h_g6_batch[i * MAX_G6_LEN], time_found);
                    fflush(stdout);
                    id++;
                }
            }
            current_batch_size = 0;
        }
    }

    if (current_batch_size > 0) {
        CUDA_CHECK(cudaMemcpy(d_g6_batch, h_g6_batch, current_batch_size * MAX_G6_LEN * sizeof(char), cudaMemcpyHostToDevice));
        
        int blocksPerGrid = (current_batch_size + threadsPerBlock - 1) / threadsPerBlock;
        processGraphsKernel<<<blocksPerGrid, threadsPerBlock>>>(d_g6_batch, d_results, current_batch_size);
        CUDA_CHECK(cudaMemcpy(h_results, d_results, current_batch_size * sizeof(int), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaDeviceSynchronize());

        clock_t current_time = clock();
        double time_found = (double)(current_time - start_time) / CLOCKS_PER_SEC;

        for (int i = 0; i < current_batch_size; i++) {
            if (h_results[i]) {
                printf("%d: %s (czas pobrania paczki: %f s)\n", id, &h_g6_batch[i * MAX_G6_LEN], time_found);
                fflush(stdout);
            }
            id++;
        }
    }

    clock_t end_time = clock();
    double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("Czas dzialania dla %ld grafow: %f s\n", lines_read, time_spent);

    free(h_g6_batch);
    free(h_results);
    cudaFree(d_g6_batch);
    cudaFree(d_results);
    fclose(f);

    return 0;
}