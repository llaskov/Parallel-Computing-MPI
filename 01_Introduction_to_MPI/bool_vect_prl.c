#include <mpi.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
    /* variables */
    int world_size;         /* number of processes */
    int prc_rank;           /* rank of the process */
    int i;                  /* loop counter */
    int vect_size;          /* number of elements of the array */
    int weight;             /* accumulated weight in the process */
    int part_size;          /* number of elements per process */
    int start;              /* starting index for the current process */
    int end;                /* ending index for the current process */

    /* init */
    world_size = 0;
    prc_rank = 0;
    unsigned int vect[] = 
        {0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1};
    i = 0;
    weight = 0;
    vect_size = sizeof(vect) / sizeof(vect[0]);
    part_size = 0;
    start = 0;
    end = 0;

    /* init MPI */
    MPI_Init(&argc, &argv);

    /* get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /* get the rank of the process */
    MPI_Comm_rank(MPI_COMM_WORLD, &prc_rank);

    part_size = vect_size / world_size;
    start = prc_rank * part_size;
    end = start + part_size - 1;

    for (i = start; i <= end; i++)
    {
        if (vect[i])
        {
            weight++;
        }
    }

    printf("Process %d counted %d weight\n", prc_rank, weight);

    /* clear MPI resources */
    MPI_Finalize();

    return 0;
}/* main */
