#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <mpi.h>

#define ROOT 0              /* rank of the root process */
#define RQUST_SIZE 2        /* size of the array of communication handles (send and receive) */
#define MAX_PARTICLES 4000  /* total maximum number of particles */
#define MAX_PART_PRCS 128   /* maximum number of particles per process */
#define REPLIC 4            /* replication for contiguous data type */
#define MAX_ITER 10         /* number of iterations of the simulation */
#define EPS 1.0e-14

/* two-dimensional vector */
typedef struct
{
    double x;               /* x-coordinate */
    double y;               /* y-coordinate */
} Vector2;

/* particle structure */
typedef struct
{
    Vector2 coord;          /* coordinates of the particle in space */
    double mass;            /* mass of the particle */
} Particle;

/* particle velocity structure */
typedef struct
{
    Vector2 coord;          /* old coordinates of the particle in space */
    Vector2 force;          /* three components of the force */
} Velocity;

/* generate initial positions and velocities */
void init(
    Particle part[],       /* array of particles */
    Velocity pvel[],       /* array of velocities */
    int numb_part          /* number of particles */    
    )
{
    int i;

    for (i = 0; i < numb_part; i++)
    {
        /* random coordinates in space */
        part[i].coord.x = drand48();
        part[i].coord.y = drand48();

        /* mass of all particles is 1 */
        part[i].mass = 1.0;

        /* the corresponding velocity of the particle */
        /* old positions */
        pvel[i].coord.x = part[i].coord.x;
        pvel[i].coord.y = part[i].coord.y;
        /* force components */
        pvel[i].force.x = 0.0;
        pvel[i].force.y = 0.0;
    }
}/* init */

/* calculate forces between particles */
double calcForce(
        Particle curr[],            /* particles in the current process */
        Particle rest[],            /* the rest particles */
        Velocity pvel[],            /* particles velocities */
        int numb_part               /* number of particles */
    )
{
    int i;
    int j;
    double max_force;
    double r;
    double r_min;
    double fx;
    double fy;

    max_force = 0.0;
    r = 0.0;
    r_min = 0.0;
    fx = 0.0;
    fy = 0.0;

    for (i = 0; i < numb_part; i++)
    {
        r_min = INT_MAX;
        fx = 0.0;
        fy = 0.0;
        for (j = 0; j < numb_part; j++)
        {
            r = (curr[i].coord.x - rest[j].coord.x) * (curr[i].coord.x - rest[j].coord.x)
                + (curr[i].coord.y - rest[j].coord.y) * (curr[i].coord.y - rest[j].coord.y);

            /* check if the particles do not coincide */
            if (fabs(r) >= EPS)
            {
                /* set minimum distance */
                if (r < r_min)
                {
                    r_min = r;
                }

                /* calculate forces */
                fx -=  rest[j].mass * (curr[i].coord.x - rest[j].coord.x) / r * sqrt(r);
                fy -=  rest[j].mass * (curr[i].coord.y - rest[j].coord.y) / r * sqrt(r);
            }
        }
        pvel[i].force.x += fx;
        pvel[i].force.y += fy;

        /* estimate of (1 / m) |df / dx| */
        fx = sqrt(fx * fx + fy * fy) / r_min;
        if (fx > max_force)
        {
            max_force = fx;
        }
    }
 
    return max_force;
}//force

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
    double time;                        /* computation time */

    /* algorithm variables */
    int i;                              /* counter */
    int pipe;                           /* counter for the current process in the pipe */
    int numb_part;                      /* number of particles per process */
    int totl_part;                      /* total number of particles */
    double max_force;                   /* maximum force */
    double max_force_seg;               /* maximum force per segment ???? */
    int counts[MAX_PART_PRCS];          /* number of particles on each process */
    int displc[MAX_PART_PRCS];          /* number of particles on each process */
    Particle part[MAX_PARTICLES];       /* array of particles for all nodes */
    Velocity pvel[MAX_PARTICLES];       /* array of velocities of all particles */
    Particle sbuf[MAX_PARTICLES];       /* send buffer of particles */
    Particle rbuf[MAX_PARTICLES];       /* receive buffer of particles */

    /* init */
    comm_ring = 0;
    part_type = 0;
    numb_prcs = 0;
    curr_rank = 0;
    periodic = 0;
    rank_left = 0;
    rank_rght = 0;
    time = 0.0;
    i = 0;
    numb_part = 0;
    max_force = 0.0;
    max_force_seg = 0.0;   
    
    /* init MPI */
    MPI_Init(&argc, &argv);

    /* get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &numb_prcs);

    /* get the rank of the current process */
    MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);

    /* handle command line arguments */
    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s <number of particles>\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    numb_part = atoi(argv[1]) / numb_prcs;

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

    /* create 1D periodic Cartesian communicator */
    periodic = 1;           /* the communicator is periodic, communicator is a ring */
    MPI_Cart_create(
        MPI_COMM_WORLD,     /* input communicator */
        1,                  /* number of dimensions */
        &numb_prcs,         /* number of processes in the dimension */
        &periodic,          /* true if dimension is periodic */
        1,                  /* true to reorder processes */
        &comm_ring          /* resulting Cartesian communicator that is ring */
    );

    /* get the ranks of the right and left neighbors in the ring */
    MPI_Cart_shift(
        comm_ring,          /* ring communicator */
        0,                  /* coordinate dimension of shift */
        1,                  /* displacement */
        &rank_left,         /* rank of the left neighbor */
        &rank_rght          /* rank of the right neighbor */
    );

    /* create contiguous datatype for particles */
    MPI_Type_contiguous(
        REPLIC,             /* replication */
        MPI_DOUBLE,         /* source data type */
        &part_type          /* resulting contiguous datatype */
    );
    MPI_Type_commit(&part_type);

    /* each process gets the number of particles in each process */
    MPI_Allgather(
        &numb_part,         /* send buffer */
        1,                  /* number of elements in send buffer */
        MPI_INT,            /* send data type */
        counts,             /* number of particles on each processes */
        1,                  /* number of elements in receive buffer */
        MPI_INT,            /* receive data type */
        comm_ring           /* communicator */       
    );

    /* calculate displacements */
    displc[0];
    for (i = 1; i < numb_prcs; i++)
    {
        displc[i] = displc[i - 1] + counts[i - 1];
    }
    /* and total number of particles; each process could have different number of particles */
    totl_part = displc[numb_prcs - 1] + counts[numb_prcs - 1];
    
    /* initial positions and velocities */
    init(part, pvel, numb_part);

    time = MPI_Wtime();                 /* starting time */

    /* iterate simulation */
    for (i = 0; i < MAX_ITER; i++)
    {
        /* copy particles into send buffer */
        memcpy(sbuf, part, numb_part * sizeof(Particle));

        /* total maximum force */
        max_force = 0.0;
        for (pipe = 0; pipe < numb_prcs; pipe++)
        {
            if (pipe != numb_prcs - 1)
            {
                /* non-blocking send */
                MPI_Isend(
                    sbuf,           /* send buffer of particles */
                    numb_part,      /* number of particles */
                    part_type,      /* contiguous particle data type   */
                    rank_rght,      /* destination is rank of the right neighbor in the ring */
                    pipe,           /* message tag is the pipe index */           
                    comm_ring,      /* ring communicator */
                    &rqust[0]       /* handler to the request */
                );

                /* non-blocking receive */
                MPI_Irecv(
                    rbuf,           /* receive buffer of particles */
                    numb_part,      /* number of particles */
                    part_type,      /* contiguous particle data type   */
                    rank_left,      /* source is rank of the left neighbor in the ring */
                    pipe,           /* message tag is the pipe index */           
                    comm_ring,      /* ring communicator */
                    &rqust[1]       /* handler to the request */
                );
            }

            /* compute forces */
            max_force_seg = calcForce(part, sbuf, pvel, numb_part);

            /* current maximum force */
            if (max_force_seg > max_force)
            {
                max_force = max_force_seg;
            }

            /* push pipe */
            if (pipe != numb_prcs - 1)
            {
                /* wait for the send and receive communications to complete */
                MPI_Waitall(
                    RQUST_SIZE,     /* size of arrays */ 
                    rqust,          /* array of requests */
                    stats           /* array of statuses */
                );
            }

            /* copy the received into the send buffers */
            memcpy(sbuf, rbuf, counts[pipe] * sizeof(Particle));
        }
        printf("%f\n", max_force);

        /* sim_t <- new positions */
        /* graphics display of particles on the screen */
    }

    time = MPI_Wtime() - time;          /* ending time */

    /* output execution */
    if (curr_rank == ROOT)
    {
        printf("Number of particles: %d. Execution time: %f sec.\n", totl_part, time);
    }

    /* free contiguous data type */
    MPI_Type_free(&part_type);

    /* clear MPI resources */
    MPI_Finalize();

    return 0;
}/* main */
