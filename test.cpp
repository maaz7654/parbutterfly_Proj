#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    int rank, size;
    
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    
    // Get the rank of the process and the total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Print the rank and total number of processes
    printf("Hello from rank %d of %d\n", rank, size);

    // Finalize the MPI environment
    MPI_Finalize();
    
    return 0;
}
