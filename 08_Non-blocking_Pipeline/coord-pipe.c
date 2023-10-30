#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define ROOT 0              /* rank of the root process */
#define RQUST_SIZE 2        /* size of the array of communication handles (send and receive) */
#define MAX_PARTICLES 4000  /* total maximum number of particles */
#define MAX_PART_PRCS 128   /* maximum number of particles per process */
#define REPLIC 4            /* replication for contiguous data type */
#define MAX_ITER 10         /* number of iterations of the simulation */

/* point in the Cartesian plane */
typedef struct
{
    double x;               /* x-coordinate */
    double y;               /* y-coordinate */
} Coord;

/* print an array */
void printArr(
    int arr[],
    int size
    )
{
    int i;

    for (i = 0; i < size; i++)
    {
        printf("[%d]: %d\n", i, arr[i]);
    }
}/* print */

/* print array of particles coordinates */
void printPart(
    Coord part[],
    int size
    )
{
    int i;

    for (i = 0; i < size; i++)
    {
        printf("%d: (%f, %f)\n", i, part[i].x, part[i].y);
    }
}/* printPart */

/* read number of particles per process from the command line arguments */
int readNubPart(
    int argc,               /* number of command line arguments */
    char* argv[],           /* command line arguments */
    int numb_prcs           /* number of processes */
    )
{
    int numb_part;

    numb_part = 0;

    /* handle command line arguments */
    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s <number of particles>\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* number of particles per process */
    numb_part = atoi(argv[1]) / numb_prcs;

    /* number of particles per process must not be too big */
    if (numb_part * numb_prcs > MAX_PARTICLES)
    {
        fprintf(
            stderr, 
            "%d number of particles is more than the max %d\n", 
            numb_part * numb_prcs, 
            MAX_PARTICLES
            );
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return numb_part;
}/* readNubPart */

/* create pipe ring communicator */
void createRing(
        MPI_Comm* comm_ring,    /* 1D Cartesian communicator forms a ring */
        int* rank_left,         /* rank of the left neighbor */
        int* rank_rght,         /* rank of the right neighbor */
        int* numb_prcs          /* number of processes */
        
    )
{
    int periodic;               /* 0 not periodic, 1 periodic Cartesian communicator */
    
    /* create 1D periodic Cartesian communicator */
    periodic = 1;           /* the communicator is periodic, communicator is a ring */
    MPI_Cart_create(
        MPI_COMM_WORLD,     /* input communicator */
        1,                  /* number of dimensions */
        numb_prcs,          /* number of processes in the dimension */
        &periodic,          /* true if dimension is periodic */
        1,                  /* true to reorder processes */
        comm_ring           /* resulting Cartesian communicator that is ring */
    );

    /* get the ranks of the right and left neighbors in the ring */
    MPI_Cart_shift(
        *comm_ring,         /* ring communicator */
        0,                  /* coordinate dimension of shift */
        1,                  /* displacement */
        rank_left,          /* rank of the left neighbor */
        rank_rght           /* rank of the right neighbor */
    );
}/* createRing */

/* generate random particles as coordinates in the plane */
void randPart(
    Coord part[],           /* array of particles */
    int numb_part           /* number of particles */    
    )
{
    int i;

    for (i = 0; i < numb_part; i++)
    {
        /* pseudo-random coordinates in space */
        /* uniformly distributed in [0.0, 1.0) */
        part[i].x = drand48();
        part[i].y = drand48();
    }
}/* init */

double calcCurrDist(
    Coord curr[],           /* particles in the current process */
    Coord rest[],           /* the rest particles */
    int numb_part           /* number of particles */
    )
{
    int i, j;
    double max_dist;
    double loc_dist;

    i = 0;
    j = 0;
    max_dist = 0.0;
    loc_dist = 0.0;

    for (i = 0; i < numb_part; i++)
    {
        for (j = 0; j < numb_part; j++)
        {
            loc_dist = 
                sqrt( 
                    (curr[i].x - rest[j].x) * (curr[i].x - rest[j].x)
                    + (curr[i].y - rest[j].y) * (curr[i].y - rest[j].y)
                );
            max_dist = loc_dist > max_dist ? loc_dist : max_dist;
        }
    }

    return max_dist;
}/* calcCurrDist */

int main(int argc, char* argv[])
{
    /* MPI variables */
    MPI_Comm comm_ring;                 /* 1D Cartesian communicator forms a ring */
    MPI_Datatype part_type;             /* contiguous datatype for particles */
    MPI_Request rqust[RQUST_SIZE];      /* array of communication handles */
    MPI_Status stats[RQUST_SIZE];       /* array of statuses */
    int numb_prcs;                      /* number of processes */
    int curr_rank;                      /* rank of the current process */
    int periodic;                       /* 0 not periodic, 1 periodic Cartesian communicator */
    int rank_left;                      /* rank of the left neighbor */
    int rank_rght;                      /* rank of the right neighbor */
    double exec_time;                   /* computation time */

    /* algorithm variables */
    int indx_pipe;                      /* counter for the current process in the pipe */
    int numb_part;                      /* number of particles per process */
    double max_dist;                    /* maximum distance */
    double max_dist_curr;               /* maximum distance per segment */
    int part_per_prcs[MAX_PART_PRCS];   /* number of particles on each process */
    Coord part[MAX_PARTICLES];          /* array of particles for all nodes */
    Coord sbuf[MAX_PARTICLES];          /* send buffer of particles */
    Coord rbuf[MAX_PARTICLES];          /* receive buffer of particles */

    /* init */
    comm_ring = 0;
    part_type = 0;
    numb_prcs = 0;
    curr_rank = 0;
    periodic = 0;
    rank_left = 0;
    rank_rght = 0;
    exec_time = 0.0;
    indx_pipe = 0;
    numb_part = 0;
    max_dist = 0.0;
    max_dist_curr = 0.0;

    /* init MPI */
    MPI_Init(&argc, &argv);

    /* get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &numb_prcs);

    /* get the rank of the current process */
    MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);

    /* root process reads number of particles and broadcasts it */
    if (curr_rank == ROOT)
    {
        /* handle command line arguments */
        numb_part = readNubPart(argc, argv, numb_prcs);
    }
    MPI_Bcast(
        &numb_part,             /* data to be communicated */
        1,                      /* number of items */
        MPI_INT,                /* data type of communicated data */
        ROOT,                   /* rank of the sender process */
        MPI_COMM_WORLD          /* communicator */
    );

    /* 1D periodic Cartesian communicator: ring pipeline */
    createRing(&comm_ring, &rank_left, &rank_rght, &numb_prcs);

    /* each process gets the number of particles in each process */
    /* it is possible to rework that number of particles are not equally distributed among processes */
    MPI_Allgather(
        &numb_part,         /* send buffer */
        1,                  /* number of elements in send buffer */
        MPI_INT,            /* send data type */
        part_per_prcs,      /* number of particles on each processes */
        1,                  /* number of elements in receive buffer */
        MPI_INT,            /* receive data type */
        comm_ring           /* communicator */       
    );

    /* root process reads number of particles and broadcasts it */
    if (curr_rank == ROOT)
    {
        /* print the number of particles per process */
        printf("Number of particles per process:\n");
        printArr(part_per_prcs, numb_prcs);
    }

    /* create contiguous datatype for particles */
    MPI_Type_contiguous(
        REPLIC,             /* replication */
        MPI_DOUBLE,         /* source data type */
        &part_type          /* resulting contiguous datatype */
    );
    MPI_Type_commit(&part_type);

    /* initialize particles with random coordinates */
    srand48(time(0) * (curr_rank + 1));
    randPart(part, numb_part);

    printf("Particles in %d\n", curr_rank);
    printPart(part, numb_part);

    exec_time = MPI_Wtime();    /* starting time */

    /* implement the pipeline */
    /* copy particles into send buffer */
    memcpy(sbuf, part, numb_part * sizeof(Coord));
    max_dist = 0.0;
    for (indx_pipe = 0; indx_pipe <  numb_prcs; indx_pipe++)
    {
        if (indx_pipe != numb_prcs - 1)
        {
            /* non-blocking send */
            MPI_Isend(
                sbuf,           /* send buffer of particles */
                numb_part,      /* number of particles */
                part_type,      /* contiguous particle data type */
                rank_rght,      /* destination is rank of the right neighbor in the ring */
                indx_pipe,      /* message tag is the pipe index */           
                comm_ring,      /* ring communicator */
                &rqust[0]       /* handler to the request */
            );

            /* non-blocking receive */
            MPI_Irecv(
                rbuf,           /* receive buffer of particles */
                numb_part,      /* number of particles */
                part_type,      /* contiguous particle data type */
                rank_left,      /* source is rank of the left neighbor in the ring */
                indx_pipe,      /* message tag is the pipe index */           
                comm_ring,      /* ring communicator */
                &rqust[1]       /* handler to the request */
            );
        }

        /* calculate distances with the known points */
        max_dist_curr = calcCurrDist(part, sbuf, numb_part);
        
        /* current maximum in the pipe */
        max_dist = max_dist_curr > max_dist ? max_dist_curr : max_dist;

        /* push the pipe */
        if (indx_pipe != numb_prcs - 1)
        {
            /* block until send and receive communications to complete */
            MPI_Waitall(
                RQUST_SIZE,     /* size of arrays */ 
                rqust,          /* array of requests */
                stats           /* array of statuses */
            );
        }

        /* copy the received into the send buffers */
        memcpy(sbuf, rbuf, part_per_prcs[indx_pipe] * sizeof(Coord));
    }

    exec_time = MPI_Wtime() - exec_time;    /* ending time */

    /* output execution */
    if (curr_rank == ROOT)
    {
        printf("Max distance: %f. Execution time: %f sec.\n", max_dist, exec_time);
    }

    /* free contiguous data type */
    MPI_Type_free(&part_type);    

    /* clear MPI resources */
    MPI_Finalize();

    return 0;
}/* main */
