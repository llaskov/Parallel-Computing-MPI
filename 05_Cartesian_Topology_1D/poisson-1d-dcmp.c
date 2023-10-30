#include <mpi.h>
#include <stdio.h>

/* constants */
#define MAXN 16                 /* size of the problem */
#define CPCT MAXN + 2           /* arrays capacity */
#define MAXI 10000              /* maximum number of iterations */
#define EPS 1.0e-5              /* accuracy */

/* macro to find minimum */
#define MIN(x, y) (x < y ? x : y)

/* compute start and end index in 1D decomposition for each process */
void decomp1d(int numb_intr, int numb_prcs, int curr_rank, int* start, int* end)
{
    int div;                    /* int division number rows by number of processes */
    int rem;                    /* integer division reminder to correct */

    div = 0;
    rem = 0;

    div = numb_intr / numb_prcs;    /* number of rows per process */
    *start = curr_rank * div + 1;   /* starting index of the current process */
    rem = numb_intr % numb_prcs;    /* number of rows may not be evenly divisible by the processes */
    *start += MIN(curr_rank, rem);  /* correct starting index */
    if (curr_rank < rem)
    {
        div++;
    }
    *end = *start + div - 1;
    if (*end > numb_intr || curr_rank == numb_prcs - 1)
    {
        *end = numb_intr;
    }
}/* decomp1d */

/* initialize the right-hand side f, and initial solution guess a and b */
/* f is initialized all with zeroes */
/* a and b: 0 inside, 1 left boundary, 1 upper boundary, 0 right boundary, 0 down boundary */
/* note that upper left [0][0] and upper right corners [0][numb_intr + 1] are both zero */
void init1d(
        double a[][CPCT],       /* initial solution guess */ 
        double b[][CPCT],       /* solution estimation after iteration */
        double f[][CPCT],       /* right-hand side */ 
        int numb_intr,          /* number of intervals */ 
        int start,              /* start row of the matrix decomposition */               
        int end                 /* end row of the matrix decomposition */
        )
{
    int i;
    int j;

    i = 0;
    j = 0;

    /* values inside the domain */
    for (i = start - 1; i <= end + 1; i++)
    {
        for (j = 0; j <= numb_intr + 1; j++)
        {
            a[i][j] = 0.0;
            b[i][j] = 0.0;
            f[i][j] = 0.0;
        }
    }

    /* boundary conditions */
    for (i = start; i <= end; i++)
    {
        a[i][0] = 1.0;
        b[i][0] = 1.0;
        a[i][numb_intr + 1] = 0.0;
        b[i][numb_intr + 1] = 0.0;
    }

    if (start == 1)
    {
        for (i = 1; i <= numb_intr; i++)
        {
            a[0][i] = 1.0;
            b[0][i] = 1.0;
        }
    }
}/* init1d */

/* exchange ghost points */
void exchng1d(
        double a[][CPCT],       /* initial solution guess */
        int numb_intr,          /* number of intervals */
        int start,              /* start row of the matrix decomposition */               
        int end,                /* end row of the matrix decomposition */
        MPI_Comm comm1d,        /* 1D Cartesian communicator */
        int belw_rank,          /* rank of the neighbor below */
        int abov_rank           /* rank of the neighbor above */
        )
{
    MPI_Sendrecv(
            a[end],                 /* initial address of send buffer */
            numb_intr,              /* number of elements to send */
            MPI_DOUBLE_PRECISION,   /* data type of send elements */
            abov_rank,              /* destination rank */
            0,                      /* send tag */
            a[start - 1],           /* initial address of receive buffer */
            numb_intr,              /* maximum number of elements to receive */
            MPI_DOUBLE_PRECISION,   /* data type of receive elements */
            belw_rank,              /* source rank */
            0,                      /* receive tag */
            comm1d,                 /* communicator */
            MPI_STATUS_IGNORE       /* status object: ignore it */            
            );
    MPI_Sendrecv(
            a[start],               /* initial address of send buffer */
            numb_intr,              /* number of elements to send */
            MPI_DOUBLE_PRECISION,   /* data type of send elements */
            belw_rank,              /* destination rank */
            1,                      /* send tag */
            a[end + 1],             /* initial address of receive buffer */
            numb_intr,              /* maximum number of elements to receive */
            MPI_DOUBLE_PRECISION,   /* data type of receive elements */
            abov_rank,              /* source rank */
            1,                      /* receive tag */
            comm1d,                 /* communicator */
            MPI_STATUS_IGNORE       /* status object: ignore it */            
            );
}/* exchng1d */

/* Jacobi sweep for 1D decomposition */
/* a contains initial estimation, b contains the result */
void jsweep1d(
        double a[][CPCT],       /* initial solution guess */ 
        double b[][CPCT],       /* solution estimation after iteration */
        double f[][CPCT],       /* right-hand side */ 
        int numb_intr,          /* number of intervals */ 
        int start,              /* start row of the matrix decomposition */               
        int end                 /* end row of the matrix decomposition */
        )
{
    int i;
    int j;
    double h;

    i = 0;
    j = 0;
    h = 0.0;

    h = 1.0 / (numb_intr + 1.0);    /* size of the interval in the domain */

    /* perform the sweep */
    for (i = start; i <= end; i++)
    {
        for (j = 1; j <= numb_intr; j++)
        {
            b[i][j] =
                0.25 
                * (
                    a[i][j - 1]
                    + a[i + 1][j]
                    + a[i - 1][j]
                    + a[i][j + 1]                
                    -
                    h * h * f[i][j]
                );
        }
    }
}/* jsweep1d */

/* difference between two arrays */
double diff(
        double a[][CPCT],       /* initial solution guess */ 
        double b[][CPCT],       /* solution estimation after iteration */
        int numb_intr,          /* number of intervals */ 
        int start,              /* start row of the matrix decomposition */               
        int end                 /* end row of the matrix decomposition */
        )
{
    int i;
    int j;
    double result;
    double local_diff;

    i = 0;
    j = 0;
    result = 0.0;
    local_diff = 0.0;

    for (i = start; i <= end; i++)
    {
        for (j = 1; j <= numb_intr; j++)
        {
            local_diff = a[i][j] - b[i][j];
            result += local_diff * local_diff; 
        }
    }

    return result;
}/* diff */

void printArr(
        double arr[][CPCT],
        int numb_intr,          /* number of intervals */ 
        int start,              /* start row of the matrix decomposition */               
        int end                 /* end row of the matrix decomposition */
        )
{
    int i;
    int j;

    i = 0;
    j = 0;

    for (i = start; i <= end; i++)
    {
        for (j = 0; j < numb_intr; j++)
        {
            printf("%f ", arr[i][j]);
        }
        printf("\n");
    }
}/* printArr */

int main(int argc, char* argv[])
{
    /* MPI variables */
    MPI_Comm comm1d;        /* 1D Cartesian communicator */
    int numb_prcs;          /* number of processes */
    int curr_rank;          /* rank of the current process */
    int belw_rank;          /* rank of the neighbor below */
    int abov_rank;          /* rank of the neighbor above */
    int dime_prcs[1];       /* array with number of processes in each dimension */
    int dime_peri[1];       /* logical array, periodic or not in each dimension */
    double time_start;      /* starting time */
    double time_end;        /* ending time */

    /* algorithm variables */
    int i;                  /* loop counter */
    int more;               /* more iterations or the solution has converged */
    int numb_intr;          /* number of intervals to divide the domain / size of the problem */
    int start;              /* start row of the matrix decomposition */
    int end;                /* end row of the matrix decomposition */
    double diff_work;
    double diff_norm;
    double a[CPCT][CPCT];   /* initial solution guess */
    double b[CPCT][CPCT];   /* solution estimation after iteration */
    double f[CPCT][CPCT];   /* right-hand side */

    /* init */
    numb_prcs = 0;
    curr_rank = 0;
    belw_rank = 0;
    abov_rank = 0;
    time_start = 0.0;
    time_end = 0.0;

    i = 0;
    more = 0;
    numb_intr = 0;
    start = 0;
    end = 0;
    diff_work = 0.0;
    diff_norm = 0.0;
    
    /* init MPI */
    MPI_Init(&argc, &argv);

    /* get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &numb_prcs);

    /* get the rank of the current process */
    MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);

    /* set number of intervals / size of the problem in process with rank 0 */
    if (curr_rank == 0)
    {
        numb_intr = MAXN;
    }
    MPI_Bcast(
        &numb_intr,             /* data to be communicated */
        1,                      /* number of items */
        MPI_INTEGER,            /* data type of communicated data */
        0,                      /* rank of the sender process */
        MPI_COMM_WORLD          /* communicator */
        );

    dime_prcs[0] = numb_prcs;
    dime_peri[0] = 0;

    /* decomposition of the domain with new Cartesian communicator */
    MPI_Cart_create(
        MPI_COMM_WORLD,         /* input communicator */
        1,                      /* number of dimensions of Cartesian grid */
        dime_prcs,              /* number of processes in each dimension */
        dime_peri,              /* whether the grid is periodic */
        1,                      /* true to reorder processes */
        &comm1d                 /* the new Cartesian communicator */
        );

    /* get the position of the current process in the Cartesian communicator */
    MPI_Comm_rank(comm1d, &curr_rank);

    /* get the neighbors of the current process in the Cartesian communicator */
    /* by shifted source and destination ranks */
    MPI_Cart_shift(
        comm1d,                 /* Cartesian communicator */
        0,                      /* coordinate dimension of shift */
        1,                      /* displacement > 0 up, < 0 down */
        &belw_rank,             /* rank of the source process */
        &abov_rank              /* rank of the destination process */
        );

    /* compute start and end indexes for each process */
    decomp1d(numb_intr, numb_prcs, curr_rank, &start, &end);
    printf("Process %d, start % d, end %d\n", curr_rank, start, end);

    /* initialize the right-hand side f, and initial solution guess a and b */
    init1d(a, b, f, numb_intr, start, end);

    /* synchronization between MPI processes */
    MPI_Barrier(MPI_COMM_WORLD);    /* completes after all group members have entered the barrier */

    time_start = MPI_Wtime();       /* starting time */

    /* iterative process to converge */
    more = 1;
    for (i = 0; more && i < MAXI; i++)
    {
        exchng1d(a, numb_intr, start, end, comm1d, belw_rank, abov_rank);
        jsweep1d(a, b, f, numb_intr, start, end);
        exchng1d(b, numb_intr, start, end, comm1d, belw_rank, abov_rank);
        jsweep1d(b, a, f, numb_intr, start, end);

        diff_work = diff(a, b, numb_intr, start, end);
        MPI_Allreduce(
                &diff_work,             /* address of sending buffer */
                &diff_norm,             /* address of receive buffer */
                1,                      /* number of elements */
                MPI_DOUBLE_PRECISION,   /* data type of communicated data */
                MPI_SUM,                /* operation applied */
                comm1d                  /* communicator */
                );

        /* print out difference between a and b matrices */
        /* printf("diff_norm: %f\n", diff_norm); */

        /* if difference is less than a small constant, converged */
        if (diff_norm < EPS)
        {
            more = 0;
        }
    }

    time_end = MPI_Wtime();         /* ending time */

    if (curr_rank == 0)
    {
        if (i >= MAXI)
        {
            printf("Failed to converge.\n");
        }
        else
        {
            printf("Converged: %d iterations, %f sec.\n", 2 * i, time_end - time_start);
        }
    }

    /* clear MPI resources */
    MPI_Finalize();
    
    return 0;
}/* main */
