#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
    // constants
    const int SIZE = MPI_MAX_PROCESSOR_NAME;

    // init MPI
    MPI_Init(&argc, &argv);

    // get processor name
    char cpu_name[SIZE];    // string for processor name
    int cpu_name_len;       // processor name length
    MPI_Get_processor_name(cpu_name, &cpu_name_len);

    // get number of processes
    int world_size;         // number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // get the rank of the process
    int prc_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &prc_rank);

    cout << "CPU: " << cpu_name << ", Processes #: " << world_size << ", Process rank: " << prc_rank << endl; 

    // clear MPI resources
    MPI_Finalize();
    
    return 0;
}//main

