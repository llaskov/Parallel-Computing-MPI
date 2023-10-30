#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#define INT_IN_FILE 256
#define BYTE_IN_FILE INT_IN_FILE * sizeof(int)

/* print array */
void printArr(int* arr, int size, int curr_rank)
{
    int i;

    printf("Rank %d:\n", curr_rank);
    for (i = 0; i < size; i++)
    {
        printf("%d ", arr[i]);
    }
    printf("\n");
}/* printArr */

/* read  */

int main(int argc, char* argv[])
{
    /* MPI variables */
    MPI_File mpi_file;      /* MPI file handler */
    MPI_Status stat;        /* MPI status */ 
    int numb_prcs;          /* number of processes */
    int curr_rank;          /* rank of the current process */

    /* algorithm variables */
    int size;               /* number of bytes in the buffer */
    int numb;               /* number of integers in the buffer */
    int *buff;              /* buffer to store the elements from the ASCII file */

    /* init */
    numb_prcs = 0;
    curr_rank = 0;
    size = 0;
    buff = NULL;
    
    /* init MPI */
    MPI_Init(&argc, &argv);

    /* get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &numb_prcs);

    /* get the rank of the current process */
    MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);

    /* size in bytes, numb in integers, allocate memory */
    size = BYTE_IN_FILE / numb_prcs;
    buff = (int*) malloc(size);
    numb = INT_IN_FILE / numb_prcs;

    /* open file and read in the buffer */
    MPI_File_open(
        MPI_COMM_WORLD,     /* communicator */
        "data",             /* file name */
        MPI_MODE_RDONLY,    /* modes */
        MPI_INFO_NULL,      /* handle to the info object */
        &mpi_file           /* file handle */
    );
    MPI_File_seek(
        mpi_file,           /* file handle */
        curr_rank * size,   /* offset of the file */
        MPI_SEEK_SET        /* offset must be calculated from the head of the file */
    );
    MPI_File_read(
        mpi_file,           /* file handle */
        buff,               /* buffer to write */
        numb,               /* number of elements */
        MPI_INT,            /* data type of each buffer element */
        &stat               /* status object */
    );

    /* close the MPI file */
    MPI_File_close(&mpi_file);

    /* print the buffer */
    printArr(buff, numb, curr_rank);

    /* free the buffer */
    free(buff);

    /* clear MPI resources */
    MPI_Finalize();
    
    return 0;
}/* main */  
