#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define ROOT 0              /* rank of the root process */
#define NUMB_PRCS 2         /* number of allowed processes */
#define NUMB_BLCK 4         /* number of blocks in the structure */
#define NUMB_CARG 8         /* number of allowed command line arguments */
#define STR_CPCT 16         /* string capacity */

/* command line arguments for Mandelbrot program */
typedef struct
{
    char display[STR_CPCT]; /* name of the display */
    int max_iter;           /* maximum number of iterations */
    double xmin;            /* x-coordinate of lower left corner */
    double ymin;            /* y-coordinate of lower left corner */
    double xmax;            /* x-coordinate of upper right corner */
    double ymax;            /* y-coordinate of upper right corner */
    int wdth;               /* display width in pixels */
    int hght;               /* display height in pixels */
} CmdlArg;

/* print the CmdArg structure in the standard output */
void printArgs(const CmdlArg* ptr_cmd_arg)
{
    printf("Name of the display: %s\n", ptr_cmd_arg->display);
    printf("Maximum number of iterations: %d\n", ptr_cmd_arg->max_iter);
    printf("Coordinates of the screen: (%f, %f), (%f, %f)\n", 
            ptr_cmd_arg->xmin, 
            ptr_cmd_arg->ymin,
            ptr_cmd_arg->xmax, 
            ptr_cmd_arg->ymax
            );
    printf("Display size in pizels: (%d, %d)\n", ptr_cmd_arg->wdth, ptr_cmd_arg->hght);
}/* printCmdArg */

/* read command line arguments and store them into the structure */
void readArgs(
    int argc,               /* number of command line arguments */ 
    char* argv[],           /* command line arguments */
    CmdlArg* ptr_cmd_arg    /* pointer to the structure */
    )
{
    /* handle command line arguments */
    if (argc != NUMB_CARG + 1)
    {
        fprintf(stderr, "Usage: %s <display> <max iter> <xmin> <xmax> <ymin> <ymax> <width> <height>\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* store command line arguments into the structure */
    strcpy(ptr_cmd_arg->display, argv[1]);
    ptr_cmd_arg->max_iter = atoi(argv[2]);
    ptr_cmd_arg->xmin = atof(argv[3]);
    ptr_cmd_arg->xmax = atof(argv[4]);
    ptr_cmd_arg->ymin = atof(argv[5]);
    ptr_cmd_arg->ymax = atof(argv[6]);
    ptr_cmd_arg->wdth = atoi(argv[7]);
    ptr_cmd_arg->hght = atoi(argv[8]);
}/* readArgs */

/* create the structure datatype */
MPI_Datatype createStruct(const CmdlArg* ptr_cmd_arg)
{
    MPI_Datatype dt_cmd;            /* store the MPI structure datatype */
    MPI_Datatype types[NUMB_BLCK];  /* array of data types compose the structure */
    MPI_Aint dsplc[NUMB_BLCK];      /* array of displacements: offsets from beginning of structure */
    int block_count[NUMB_BLCK] = {  /* array with number of elements for each type */
        STR_CPCT,   /* size of the string */
        1,          /* max_iter is one int */
        4,          /* four floating point coordinates */
        2           /* two integers for display size in pixels */
    };

    /* init types */
    types[0] = MPI_CHAR;
    types[1] = MPI_INT;
    types[2] = MPI_DOUBLE;
    types[3] = MPI_INT;

    /* init displacement */
    dsplc[0] = offsetof(CmdlArg, display);
    dsplc[1] = offsetof(CmdlArg, max_iter);
    dsplc[2] = offsetof(CmdlArg, xmin);
    dsplc[3] = offsetof(CmdlArg, wdth);

    /* build the new structure datatype */
    MPI_Type_create_struct(
            NUMB_BLCK,      /* number of blocks */
            block_count,    /* number of elements in each block */
            dsplc,          /* displacements */
            types,          /* data types */
            &dt_cmd         /* handle of the datatype */
            );

    /* commit datatype to the system */
    MPI_Type_commit(&dt_cmd);

    return dt_cmd;
}/* createStruct */

/* communicate data */
void exchange(int curr_rank, CmdlArg* ptr_cmd_arg)
{
    enum Roles{SEND, RECV};

    MPI_Datatype dt_cmd;            /* store the MPI structure datatype */

    /* create the structure datatype */
    dt_cmd = createStruct(ptr_cmd_arg);

    switch (curr_rank)
    {
        case SEND:            
            printf("Process %d sends:\n", curr_rank);
            printArgs(ptr_cmd_arg);
            MPI_Send(ptr_cmd_arg, 1, dt_cmd, RECV, 0, MPI_COMM_WORLD);
            break;
        case RECV:
            MPI_Recv(ptr_cmd_arg, 1, dt_cmd, SEND, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Process %d receives:\n", curr_rank);
            printArgs(ptr_cmd_arg);
            break;
        default:
            fprintf(stderr, "Wrong role code in exchange().\n");
	        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* free the datatype from the system */
    MPI_Type_free(&dt_cmd);
}/* exchange */

int main(int argc, char* argv[])
{
    /* MPI variables */
    int numb_prcs;          /* number of processes */
    int curr_rank;          /* rank of the current process */

    /* algorithm variables */
    CmdlArg cmd_arg;

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

     /* root process reads command line arguments */
    if (curr_rank == ROOT)
    {
        /* handle command line arguments */
        readArgs(argc, argv, &cmd_arg);
    }

    /* communicate data */
    exchange(curr_rank, &cmd_arg);
    
    /* clear MPI resources */
    MPI_Finalize();
    
    return 0;
}/* main */
