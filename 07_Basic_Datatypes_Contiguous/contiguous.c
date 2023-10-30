#include <mpi.h>
#include <stdio.h>

#define NUMB_PRCS 2
#define DT_SIZE 2

void exchange(int curr_rank)
{
    enum Roles{SEND, RECV};

    char buff[DT_SIZE];

    MPI_Datatype dt;

    MPI_Type_contiguous(DT_SIZE, MPI_CHAR, &dt);
    MPI_Type_commit(&dt);

    switch (curr_rank)
    {
        case SEND:
            buff[0] = 'a';
            buff[1] = '1';
            printf("Process %d sends: %c, %c\n", curr_rank, buff[0], buff[1]);
            MPI_Send(&buff, 1, dt, RECV, 0, MPI_COMM_WORLD);
            break;
        case RECV:
            MPI_Recv(&buff, 1, dt, SEND, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Process %d receives: %c, %c\n", curr_rank, buff[0], buff[1]);
            break;
        default:
            fprintf(stderr, "Wrong role code in exchange().\n");
	        MPI_Abort(MPI_COMM_WORLD, 1);

    }

    MPI_Type_free(&dt);
}/* exchange */

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

    if (numb_prcs != NUMB_PRCS)
    {
        fprintf(stderr, "Usage: run with exactly two processes.\n");
	    MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    exchange(curr_rank);
    
    /* clear MPI resources */
    MPI_Finalize();
    
    return 0;
}/* main */

