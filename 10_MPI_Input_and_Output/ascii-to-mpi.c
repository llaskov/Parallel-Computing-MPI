#include <mpi.h>
#include <stdio.h>

#define NUMB_PRCS 1
#define BUFF_CPCT 256

/* read the ASCII file and store the contents in the buffer array */
int                     /* return number of elements in the buffer */ 
    readAscii(
        int buff[]      /* buffer to store the input from the file */
    )
{
    /* variables */
    int i;              /* row index */
    int j;              /* column index */
    int k;              /* buffer index */
    FILE* in_file;      /* input .txt file */

    /* init */
    i = 0;
    j = 0;
    k = 0;
    in_file = NULL;

    /* open input .txt file */
    in_file = fopen("in.txt", "r");
    if (in_file == NULL)
    {
        fprintf(stderr, "Error opening in.txt.\n");
	    MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* read number of rows and store them on the first position */
    fscanf(in_file, "%d", &buff[k++]);

    /* read the matrix of 0 and 1 */
    for (i = 0; i < buff[0]; i++)
    {
        for (j = 0; j < buff[0]; j++)
        {
            fscanf(in_file, "%1d", &buff[k++]);
        }
    }

    /* close the ASCII file */
    fclose(in_file);

    return k;
}/* readAscii */

/* output to MPI data file */
void wfile(
        int buff[]      /* buffer to store in the data file */
    )
{
    /* variables */
    MPI_File mpi_file;      /* MPI file handler */
    MPI_Status stat;        /* MPI status */

    /* open output MPI data file */
    if (MPI_File_open(
                MPI_COMM_SELF,                      /* communicator: open a file independently of other processes */
                "data",                             /* file name */
                MPI_MODE_CREATE | MPI_MODE_WRONLY,  /* modes */
                MPI_INFO_NULL,                      /* handler to the info object */
                &mpi_file                           /* file handler */
        )
    )
    {    
        fprintf(stderr, "Error creating the MPI file.\n");
	    MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    /* set the view for the process */
    MPI_File_set_view(
            mpi_file,           /* file handler */
            0,                  /* displacement */
            MPI_INT,            /* elementary data type */ 
            MPI_INT,            /* file type */
            "native",           /* data representation */
            MPI_INFO_NULL       /* handler to the info object */
            );

    /* write MPI file to disk */
    MPI_File_write(
            mpi_file,           /* file handler */
            buff,               /* buffer to write */
            BUFF_CPCT,          /* number of elements */
            MPI_INT,            /* data type of each buffer element */
            &stat               /* status object */
            );

    /* close the MPI file */
    MPI_File_close(&mpi_file);
}/* wfile */

int main(int argc, char* argv[])
{
    /* MPI variables */
    int numb_prcs;          /* number of processes */
    int curr_rank;          /* rank of the current process */

    /* algorithm variables */
    int size;               /* number of elements in the buffer */
    int buff[BUFF_CPCT];    /* buffer to store the elements from the ASCII file */

    /* init */
    numb_prcs = 0;
    curr_rank = 0;
    size = 0;
    
    /* init MPI */
    MPI_Init(&argc, &argv);

    /* get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &numb_prcs);

    /* get the rank of the current process */
    MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);

    if (numb_prcs != NUMB_PRCS)
    {
        fprintf(stderr, "Usage: run with exactly one process.\n");
	    MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* read the ASCII file */
    size = readAscii(buff);
    printf("Size: %d, Rows: %d\n", size, buff[0]);

    /* write the MPI file */
    wfile(buff);

    /* clear MPI resources */
    MPI_Finalize();
    
    return 0;
}/* main */
