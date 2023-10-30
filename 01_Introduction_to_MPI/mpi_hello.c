#include <mpi.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
    /* constants */
    const int SIZE =  MPI_MAX_PROCESSOR_NAME;

    /* variables */
    int world_size;         /* number of processes */
    int prc_rank;           /* rank of the process */
    int cpu_name_len;       /* processor name length */
    char cpu_name[SIZE];    /* string for processor name */

    /* init */
    world_size = 0;
    prc_rank = 0;
    cpu_name_len = 0;

    /* init MPI */
    MPI_Init(&argc, &argv);

    /* get processor name */
    MPI_Get_processor_name(cpu_name, &cpu_name_len);

    /* get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /* get the rank of the process */
    MPI_Comm_rank(MPI_COMM_WORLD, &prc_rank);

    printf("CPU: %s, Processes #: %d, Process rank: %d\n", cpu_name, world_size, prc_rank);

    /* clear MPI resources */
    MPI_Finalize();

    return 0;
}//main
