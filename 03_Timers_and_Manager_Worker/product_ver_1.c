#include <mpi.h>
#include <stdio.h>

/* capacity of arrays */
#define MAX_ROWS 1000
#define MAX_COLS 1000

/* number of elements used */
#define SIZE_ROWS 100
#define SIZE_COLS 100

#define MIN(x, y) (x < y ? x : y)

/* initialize the vector with ones for the example */
void initVect(double vect[], int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        vect[i] = 1;
    }
}/* initVect */

/* initialize the matrix for the example */
void initMatr(double matr[MAX_ROWS][MAX_COLS], int rows, int cols)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            matr[i][j] = j + 1;
        }
    }
}/* initMatr */

/* print vector */
void printVect(double vect[], int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        printf("%.2f ", vect[i]);
    }
    printf("\n");
}/* print */

/* scalar product of two vectors */
double scalarProd(double frst[], double scnd[], int size)
{
    int i;
    int result;

    result = 0;
    for (i = 0; i < size; i++)
    {
        result += frst[i] * scnd[i];
    }

    return result;
}/* scalarProd */

int main(int argc, char* argv[])
{
    /* MPI variables */
    int world_size;                 /* number of processes */
    int prc_rank;                   /* rank of the process */
    int manager;                    /* rank of the manager */
    MPI_Status status;              /* status structure; information about received message */

    /* algorithm variables */
    int i, j;                       /* loop counters */
    int rows, cols;                 /* number of rows and columns */
    int numb_sent;                  /* number of sent rows */
    int more;                       /* more iterations */
    double a[MAX_ROWS][MAX_COLS];   /* matrix to multiply */
    double b[MAX_COLS];             /* vector to multiply */
    double c[MAX_ROWS];             /* resulting vector */
    double buff[MAX_COLS];          /* buffer for each sent row */
    double ans;                     /* an element of the result sent by a worker */

    /* init */
    world_size = 0;
    prc_rank = 0;
    manager = 0;
    more = 1;

    i = 0;
    j = 0;
    rows = SIZE_ROWS;
    cols = SIZE_COLS;
    numb_sent = 0;
    ans = 0.0;

    /* init MPI */
    MPI_Init(&argc, &argv);

    /* get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /* get the rank of the process */
    MPI_Comm_rank(MPI_COMM_WORLD, &prc_rank);

    if (prc_rank == manager)    /* manager */
    {
        /* manager initializes and then dispatches */

        /* initialize matrix and vector for the example */
        initMatr(a, rows, cols);
        initVect(b, cols);
        
        /* send b to each worker process */
        MPI_Bcast(
            &b,             /* data to be communicated */
            cols,           /* number of items */
            MPI_DOUBLE,     /* data type of communicated data */
            manager,        /* rank of the sender process */
            MPI_COMM_WORLD  /* communicator */
            );

        /* send a row of the matrix to each worker process */
        /* tag with row index */
        for (i = 1; i < MIN(world_size, rows); i++)
        {
            /* copy the row into the buffer */
            for (j = 0; j < cols; j++)
            {
                buff[j] = a[i - 1][j];
            }

            /* blocking send */
            MPI_Send(
                    &buff,          /* address of send buffer */
                    cols,           /* number of elements */
                    MPI_DOUBLE,     /* data type of communicated data */
                    i,              /* rank of the destination process */
                    i - 1,          /* message tag with the index of the row */
                    MPI_COMM_WORLD  /* communicator */
                    ); 

            /* count number of sent rows */
            numb_sent++;  
        }
        for (i = 0; i < rows; i++)
        {
            MPI_Recv(
                    &ans,           /* address of receive buffer */
                    1,              /* number of elements */
                    MPI_DOUBLE,     /* data type of communicated data */
                    MPI_ANY_SOURCE, /* rank of the source */
                    MPI_ANY_TAG,    /* message tag */
                    MPI_COMM_WORLD, /* communicator */
                    &status         /* status object */
                    );

            c[status.MPI_TAG] = ans;

            if (numb_sent < rows)   /* number of rows sent is less than number of rows of a */
            {
                /* copy the row into the buffer */
                for (j = 0; j < cols; j++)
                {
                    buff[j] = a[numb_sent][j];
                }

                /* blocking send */
                MPI_Send(
                        &buff,              /* address of send buffer */
                        cols,               /* number of elements */
                        MPI_DOUBLE,         /* data type of communicated data */
                        status.MPI_SOURCE,  /* rank of the destination process */
                        numb_sent,          /* message tag with the index of the row */
                        MPI_COMM_WORLD      /* communicator */
                        );

                /* count number of sent rows */
                numb_sent++;
            }
            else    /* all rows has been sent */
            {
                 printf("Tell to stop.\n");
                /* signal sender that there are no rows left */

                /* blocking send of done message */
                MPI_Send(
                        MPI_BOTTOM,         /* address of send buffer */
                        0,                  /* number of elements */
                        MPI_DOUBLE,         /* data type of communicated data */
                        status.MPI_SOURCE,  /* rank of the destination process */
                        MAX_ROWS,           /* message tag  */
                        MPI_COMM_WORLD      /* communicator */
                        );
            }  
        }           
    }
    else    /* worker */
    {
        /* workers receive the vector b */
        MPI_Bcast(
                &b,             /* data to be communicated */
                cols,           /* number of items */
                MPI_DOUBLE,     /* data type of communicated data */
                manager,        /* rank of the sender process */
                MPI_COMM_WORLD  /* communicator */
                );

        if (prc_rank < rows)
        {
            while (more)
            {
                MPI_Recv(
                    &buff,          /* address of receive buffer */
                    cols,           /* number of elements */
                    MPI_DOUBLE,     /* data type of communicated data */
                    manager,        /* rank of the source */
                    MPI_ANY_TAG,    /* message tag */
                    MPI_COMM_WORLD, /* communicator */
                    &status         /* status object */
                    );

                if (status.MPI_TAG < MAX_ROWS)
                {
                    /* calculate dot product between matrix row and column vector */
                    ans = scalarProd(buff, b, cols);

                    /* blocking send */
                    MPI_Send(
                            &ans,           /* address of send buffer */
                            1,              /* number of elements */
                            MPI_DOUBLE,     /* data type of communicated data */
                            manager,        /* rank of the destination process */
                            status.MPI_TAG, /* message tag with the index of the row */
                            MPI_COMM_WORLD  /* communicator */
                            );
                }
                else
                {
                    printf("I have to stop.\n");
                    more = 0; /* stop the loop */
                }
            }
        }
    }

    /* print the result in the manager process */
    if (prc_rank == manager)    /* manager */
    {
        printVect(c, rows);    
    }

    /* clear MPI resources */
    MPI_Finalize();

    return 0;
}/* main */
