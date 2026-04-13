#include <iostream>
#include <cuda_runtime.h>

__global__ void cuda_hello() {
    printf("Hello World z GPU! Watek: %d\n", threadIdx.x);
}

int main() {
    cuda_hello<<<1, 5>>>();
    cudaDeviceSynchronize();
    return 0;
}