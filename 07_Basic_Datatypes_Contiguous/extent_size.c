#include <mpi.h>
#include <stdio.h>

#define NUMB_DT 2

void printDataType(int curr_rank)
{
    int dt_indx;            /* index of the datatype */
    int size;               /* size of the datatype */
    char* name;             /* string for the datatype name */ 
    MPI_Aint lb;            /* lower boundary */
    MPI_Aint ex;            /* extent */
    MPI_Datatype dt;        /* MPI datatype */

    dt_indx = 0;
    size = 0;
    lb = 0;
    ex = 0;
    
    dt_indx = curr_rank % NUMB_DT;

    switch(dt_indx)
    {
        case 0:
            name = "MPI_CHAR";
            dt = MPI_CHAR;            
            break;
        case 1:
            name = "MPI_INT";
            dt = MPI_INT;
            break;
        default:
            name = "Unknown datatype index";
    }

    MPI_Type_get_extent(
        dt,
        &lb,
        &ex
    );
    MPI_Type_size(dt, &size);
    printf("Process %d: %s lower bound %lu, extent %lu, size %d bytes\n", curr_rank, name, lb, ex, size);
}/* printDataType */

int main(int argc, char* argv[])
{
    /* MPI variables */
    int numb_prcs;          /* number of processes */
    int curr_rank;          /* rank of the current process */

    /* init */
    numb_prcs = 0;
    curr_rank = 0;
    
    /* init MPI */
    MPI_Init(&argc, &argv);

    /* get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &numb_prcs);

    /* get the rank of the current process */
    MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);

    printDataType(curr_rank);    
    
    /* clear MPI resources */
    MPI_Finalize();
    
    return 0;
}/* main */
