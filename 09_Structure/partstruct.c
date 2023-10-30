#include <mpi.h>
#include <stdio.h>
#include <stddef.h>

#define NUMB_BLCK 2

typedef struct
{
    double xcoord;
    double ycoord;
    int label;
} Part;

int main(int argc, char* argv[])
{
    /* MPI variables */
    int numb_prcs;                  /* number of processes */
    int curr_rank;                  /* rank of the current process */

    MPI_Datatype dt_part;           /* store the MPI structure datatype */
    MPI_Datatype types[NUMB_BLCK];  /* array of data types compose the structure */
    MPI_Aint dsplc[NUMB_BLCK];      /* array of displacements: offsets from beginning of structure */
    int block_count[NUMB_BLCK];     /* array with number of elements for each type */

    /* algorithm variables */
    int i;
    Part send;
    Part rcve;

    /* init */
    i = 0;
    numb_prcs = 0;
    curr_rank = 0;
    send.xcoord = 0.1;
    send.ycoord = 0.2;
    send.label = 1;

    /* init MPI */
    MPI_Init(&argc, &argv);

    /* get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &numb_prcs);

    /* get the rank of the current process */
    MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);

    /* init types */
    types[0] = MPI_DOUBLE;
    types[1] = MPI_INT;

    /* init number of element per block */
    block_count[0] = 2;
    block_count[1] = 1;

    /* init displacement */
    dsplc[0] = offsetof(Part, xcoord);
    dsplc[1] = offsetof(Part, label);

    /* build the new structure datatype */
    MPI_Type_create_struct(
            NUMB_BLCK,      /* number of blocks */
            block_count,    /* number of elements in each block */
            dsplc,          /* displacements */
            types,          /* data types */
            &dt_part        /* handle of the datatype */
            );

    /* commit datatype to the system */
    MPI_Type_commit(&dt_part);

    switch (curr_rank)
    {
        case 0:
            MPI_Send(&send, 1, dt_part, 1, 0, MPI_COMM_WORLD);
            break;
        case 1:
            MPI_Recv(&rcve, 1, dt_part, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Received particle: (%f, %f), %d\n", rcve.xcoord, rcve.ycoord, rcve.label);
            break;
        default:
            fprintf(stderr, "Wrong role code in exchange().\n");
	        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* free the datatype from the system */
    MPI_Type_free(&dt_part);

    /* clear MPI resources */
    MPI_Finalize();
    
    return 0;
}/* main */
