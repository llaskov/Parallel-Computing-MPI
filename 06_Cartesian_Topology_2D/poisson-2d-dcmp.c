#include <mpi.h>
#include <stdio.h>

/* constants */
#define MAXN 16                 /* size of the problem */
#define CPCT MAXN + 2           /* arrays capacity */
#define MAXI 10000              /* maximum number of iterations */
#define EPS 1.0e-5              /* accuracy */
#define NDIM 2                  /* Cartesian decomposition number of dimensions */

/* macro to find minimum */
#define MIN(x, y) (x < y ? x : y)

/* window for each process in matrix coordinates */
struct PrcWin
{
    int sta_row;                /* start row of the matrix decomposition */               
    int end_row;                /* end row of the matrix decomposition */
    int sta_col;                /* start column of the matrix decomposition */               
    int end_col;                /* end column of the matrix decomposition */
};

/* initialize window  */
void initPrcWin(struct PrcWin* w)
{
    w->sta_row = 0;
    w->end_row = 0;
    w->sta_col = 0;
    w->end_row = 0;
}

/* neighbor processes of the current process needed for the shift */
struct NghPrc
{
    int left_rank;              /* rank of the neighbor left */
    int rght_rank;              /* rank of the neighbor right */
    int belw_rank;              /* rank of the neighbor below */
    int abov_rank;              /* rank of the neighbor above */
};

/* initialize window  */
void initNghPrc(struct NghPrc* n)
{
    n->left_rank = 0;
    n->rght_rank = 0;
    n->belw_rank = 0;
    n->abov_rank = 0;
}

/* initialize array with an integer */
void initarr(int arr[], int size, int init)
{
    int i;

    for (i = 0; i < size; i++)
    {
        arr[i] = init;
    }
}/* initarr */

/* find above, below, left and right neighbors of the current process */
/* implemented using two shifts: one shift per each dimension */
void findnghbrs(
        MPI_Comm comm2d,        /* 2D Cartesian communicator */
        struct NghPrc* n        /* ranks of the neighboring processes */
        )
{
    /* horizontal shift */
    MPI_Cart_shift(
        comm2d,                 /* Cartesian communicator */
        0,                      /* coordinate dimension of shift; 0 is x on the grid */
        1,                      /* displacement > 0 up, < 0 down */
        &n->left_rank,          /* rank of the source process */
        &n->rght_rank           /* rank of the destination process */
        );

    /* vertical shift */
    MPI_Cart_shift(
        comm2d,                 /* Cartesian communicator */
        1,                      /* coordinate dimension of shift; 1 is y on the grid */
        1,                      /* displacement > 0 up, < 0 down */
        &n->belw_rank,          /* rank of the source process */
        &n->abov_rank           /* rank of the destination process */
        );
}/* findnghbrs */

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

/* compute start and end index for both rows and columns in 2D decomposition for the current process */
void decomp2d(
        MPI_Comm comm2d,        /* 2D Cartesian communicator */
        int numb_intr,          /* size of the problem */
        struct PrcWin* w        /* process window */
        )
{
    int dime_prcs[NDIM];        /* array with number of processes in each dimension */
    int dime_peri[NDIM];        /* logical array, periodic or not in each dimension */
    int curr_crds[NDIM];        /* coordinates of the current process in Cartesian structure */

    initarr(dime_prcs, NDIM, 0);
    initarr(dime_peri, NDIM, 0);
    initarr(curr_crds, NDIM, 0);

    /* get Cartesian topology information */
    MPI_Cart_get( 
            comm2d,             /* communicator */
            NDIM,               /* length of the arrays */
            dime_prcs,          /* processes for each Cartesian dimension */
            dime_peri,          /* periodicity (true/false) for each Cartesian dimension */
            curr_crds           /* coordinates of calling process in Cartesian structure */
            );

    /* compute start and end index in rows */
    decomp1d(numb_intr, dime_prcs[0], curr_crds[0], &w->sta_col, &w->end_col);

    /* compute start and end index in columns */
    decomp1d(numb_intr, dime_prcs[1], curr_crds[1], &w->sta_row, &w->end_row);
}/* decomp2d */

/* initialize the right-hand side f, and initial solution guess a and b */
/* f is initialized all with zeroes */
/* a and b: 0 inside, 1 left boundary, 1 upper boundary, 0 right boundary, 0 down boundary */
/* note that upper left [0][0] and upper right corners [0][numb_intr + 1] are both zero */
void init2d(
        double a[][CPCT],       /* initial solution guess */ 
        double b[][CPCT],       /* solution estimation after iteration */
        double f[][CPCT],       /* right-hand side */ 
        int numb_intr,          /* number of intervals */ 
        struct PrcWin w         /* process window */
        )
{
    int i;
    int j;

    i = 0;
    j = 0;

    /* values inside the domain */
    for (i = w.sta_row - 1; i <= w.end_row + 1; i++)
    {
        for (j = w.sta_col - 1; j <= w.end_col + 1; j++)
        {
            a[i][j] = 0.0;
            b[i][j] = 0.0;
            f[i][j] = 0.0;
        }
    }

    /* boundary conditions */
    if (w.sta_col == 1)
    {
        for (i = w.sta_row; i <= w.end_row; i++)
        {
            a[i][0] = 1.0;
            b[i][0] = 1.0;
        }
    }
    if (w.end_col == numb_intr)
    {
        for (i = w.sta_row; i <= w.end_row; i++)
        {
            a[i][numb_intr + 1] = 0.0;
            b[i][numb_intr + 1] = 0.0;
        }
    }
    if (w.sta_row == 1)
    {
        for (i = w.sta_col; i <= w.end_col; i++)
        {
            a[0][i] = 1.0;
            b[0][i] = 1.0;
        }
    }
}/* init2d */

/* exchange ghost points */
void exchng2d(
        double a[][CPCT],       /* initial solution guess */
        struct PrcWin w,        /* process window */
        MPI_Comm comm2d,        /* 2D Cartesian communicator */
        MPI_Datatype svct_type, /* stride vector data type */
        struct NghPrc n         /* ranks of the neighboring processes */
        )
{
    int intr_col;               /* size of the window in columns */
    double* send_bfr;           /* send buffer address */
    double* recv_bfr;           /* receive buffer address */

    intr_col = w.end_col - w.sta_col + 1;

    /* ghost points above and below the current window, similar to 1D version */
    send_bfr = &(a[w.end_row][w.sta_col]);
    recv_bfr = &(a[w.sta_row - 1][w.sta_col]);
    MPI_Sendrecv(
            send_bfr,               /* initial address of send buffer */
            intr_col,               /* number of elements to send */
            MPI_DOUBLE_PRECISION,   /* data type of send elements */
            n.abov_rank,            /* destination rank */
            0,                      /* send tag */
            recv_bfr,               /* initial address of receive buffer */
            intr_col,               /* maximum number of elements to receive */
            MPI_DOUBLE_PRECISION,   /* data type of receive elements */
            n.belw_rank,            /* source rank */
            0,                      /* receive tag */
            comm2d,                 /* communicator */
            MPI_STATUS_IGNORE       /* status object: ignore it */            
            );
    send_bfr = &(a[w.sta_row][w.sta_col]);
    recv_bfr = &(a[w.end_row + 1][w.sta_col]);
    MPI_Sendrecv(
            send_bfr,               /* initial address of send buffer */
            intr_col,               /* number of elements to send */
            MPI_DOUBLE_PRECISION,   /* data type of send elements */
            n.belw_rank,            /* destination rank */
            1,                      /* send tag */
            recv_bfr,               /* initial address of receive buffer */
            intr_col,               /* maximum number of elements to receive */
            MPI_DOUBLE_PRECISION,   /* data type of receive elements */
            n.abov_rank,            /* source rank */
            1,                      /* receive tag */
            comm2d,                 /* communicator */
            MPI_STATUS_IGNORE       /* status object: ignore it */            
            );

    /* ghost points left and right from the current window, based on stride vector data type */
    send_bfr = &(a[w.sta_row][w.end_col]);
    recv_bfr = &(a[w.sta_row][w.sta_col - 1]);
    MPI_Sendrecv(
            send_bfr,               /* initial address of send buffer */
            1,                      /* number of elements to send */
            svct_type,              /* data type of send elements */
            n.rght_rank,            /* destination rank */
            0,                      /* send tag */
            recv_bfr,               /* initial address of receive buffer */
            1,                      /* maximum number of elements to receive */
            svct_type,              /* data type of receive elements */
            n.left_rank,            /* source rank */
            0,                      /* receive tag */
            comm2d,                 /* communicator */
            MPI_STATUS_IGNORE       /* status object: ignore it */            
            );
    send_bfr = &(a[w.sta_row][w.end_col]);
    recv_bfr = &(a[w.sta_row][w.sta_col + 1]);    
    MPI_Sendrecv(
            send_bfr,               /* initial address of send buffer */
            1,                      /* number of elements to send */
            svct_type,              /* data type of send elements */
            n.left_rank,            /* destination rank */
            1,                      /* send tag */
            recv_bfr,               /* initial address of receive buffer */
            1,                      /* maximum number of elements to receive */
            svct_type,              /* data type of receive elements */
            n.rght_rank,            /* source rank */
            1,                      /* receive tag */
            comm2d,                 /* communicator */
            MPI_STATUS_IGNORE       /* status object: ignore it */            
            );
}/* exchng2d */

/* Jacobi sweep for 2D decomposition */
/* a contains initial estimation, b contains the result */
void jsweep2d(
        double a[][CPCT],       /* initial solution guess */ 
        double b[][CPCT],       /* solution estimation after iteration */
        double f[][CPCT],       /* right-hand side */ 
        int numb_intr,          /* number of intervals */ 
        struct PrcWin w         /* process window */
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
    for (i = w.sta_row; i <= w.end_row; i++)
    {
        for (j = w.sta_col; j <= w.end_col; j++)
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
}/* jsweep2d */

/* difference between two arrays */
double diff2d(
        double a[][CPCT],       /* initial solution guess */ 
        double b[][CPCT],       /* solution estimation after iteration */ 
        struct PrcWin w         /* process window */
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

    for (i = w.sta_row; i <= w.end_row; i++)
    {
        for (j = w.sta_col; j <= w.end_col; j++)
        {
            local_diff = a[i][j] - b[i][j];
            result += local_diff * local_diff; 
        }
    }

    return result;
}/* diff */

int main(int argc, char* argv[])
{
    /* MPI variables */
    MPI_Comm comm2d;        /* 2D Cartesian communicator */
    MPI_Datatype svct_type; /* stride vector data type */
    int numb_prcs;          /* number of processes */
    int curr_rank;          /* rank of the current process */
    struct NghPrc ngh_prc;  /* ranks of the neighboring processes */
    int dime_prcs[NDIM];    /* array with number of processes in each dimension */
    int dime_peri[NDIM];    /* logical array, periodic or not in each dimension */
    double time_start;      /* starting time */
    double time_end;        /* ending time */

    /* algorithm variables */
    int i;                  /* loop counter */
    int more;               /* more iterations or the solution has converged */
    int numb_intr;          /* number of intervals to divide the domain / size of the problem */
    struct PrcWin prc_win;  /* window for each process in matrix coordinates */
    double diff_work;
    double diff_norm;
    double a[CPCT][CPCT];   /* initial solution guess */
    double b[CPCT][CPCT];   /* solution estimation after iteration */
    double f[CPCT][CPCT];   /* right-hand side */

    /* init */
    numb_prcs = 0;
    curr_rank = 0;
    initNghPrc(&ngh_prc);
    time_start = 0.0;
    time_end = 0.0;

    i = 0;
    more = 0;
    numb_intr = 0;
    initPrcWin(&prc_win);
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

    /* decompose domain in Cartesian domain */
    initarr(dime_prcs, NDIM, 0);/* init array with number of processes per dimension */
    MPI_Dims_create(            /* creates a division of processors in a Cartesian grid */
        numb_prcs,              /* number of processes available */
        NDIM,                   /* number of Cartesian dimensions */
        dime_prcs               /* number of processes per dimension */ 
        ); 
    initarr(dime_peri, NDIM, 0);/* dimensions are not periodic */ 
    MPI_Cart_create(            /* new Cartesian communicator */
        MPI_COMM_WORLD,         /* input communicator */
        NDIM,                   /* number of dimensions of Cartesian grid */
        dime_prcs,              /* number of processes in each dimension */
        dime_peri,              /* whether the grid is periodic */
        1,                      /* true to reorder processes */
        &comm2d                 /* the new Cartesian communicator */
        );

    /* get the position of the current process in the Cartesian communicator */
    MPI_Comm_rank(comm2d, &curr_rank);

    /* find above, below, left and right neighbors of the current process */
    findnghbrs(comm2d, &ngh_prc);

    /* find the window of the current process in the 2D mesh */
    decomp2d(comm2d, numb_intr, &prc_win);

    /* stride vector data type to exchange ghost points in horizontal direction */
    MPI_Type_vector(            /* create the data type */
        prc_win.end_row - prc_win.sta_row + 1,  /* number of blocks */
        1,                      /* number of elements in each block */
        prc_win.end_col - prc_win.sta_col + 3,  /* stride: number of elements between start of each block */        
        MPI_DOUBLE_PRECISION,   /* old data type */
        &svct_type              /* resulting stride vector data type */
        );
    MPI_Type_commit(&svct_type);/* commit the data type to the system */

    /* initialize the right-hand side f, and initial solution guess a and b */
    init2d(a, b, f, numb_intr, prc_win);

     /* synchronization between MPI processes */
    MPI_Barrier(MPI_COMM_WORLD);    /* completes after all group members have entered the barrier */

    time_start = MPI_Wtime();       /* starting time */

    /* iterative process to converge */
    more = 1;
    for (i = 0; more && i < MAXI; i++)
    {
        exchng2d(a, prc_win, comm2d, svct_type, ngh_prc);
        jsweep2d(a, b, f, numb_intr, prc_win);
        exchng2d(b, prc_win, comm2d, svct_type, ngh_prc);
        jsweep2d(b, a, f, numb_intr, prc_win);

        diff_work = diff2d(a, b, prc_win);

        MPI_Allreduce(
                &diff_work,             /* address of sending buffer */
                &diff_norm,             /* address of receive buffer */
                1,                      /* number of elements */
                MPI_DOUBLE_PRECISION,   /* data type of communicated data */
                MPI_SUM,                /* operation applied */
                comm2d                  /* communicator */
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
            printf("Failed to converge with difference %f.\n", diff_norm);
        }
        else
        {
            printf("Converged: %d iterations, %f sec.\n", 2 * i, time_end - time_start);
        }
    }

    /* clear the stride vector data type */
    MPI_Type_free(&svct_type);

    /* clear MPI resources */
    MPI_Finalize();
    
    return 0;
}/* main */
