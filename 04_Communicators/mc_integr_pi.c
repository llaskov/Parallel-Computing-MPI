#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>

/* number of random numbers */
#define RAND_ARR_SIZE 1000

/* maximum number of random points */
#define MAX_POINTS 1000000000

/* message tags */
#define RQUST 1
#define REPLY 2

int main(int argc, char* argv[])
{
    /* MPI variables */
    int numb_prcs;          /* number of processes */
    int curr_rank;          /* rank of the current process */
    int srvr_rank;          /* rank of the server process */
    int arr_ranks[1];       /* array of ranks of processes to exclude from a group */

    MPI_Comm wrld;          /* communicator for world */
    MPI_Comm wrkr;          /* communicator for workers */
    MPI_Group wrld_group;   /* world group of processes */
    MPI_Group wrkr_group;   /* workers group of processes */
    MPI_Status status;      /* status structure; information about received message */

    /* algorithm variables */
    int i;                  /* loop counter */
    int iter;               /* count iterations */
    int in;                 /* points inside circle */
    int out;                /* points outside circle */
    int done;               /* when to stop approximation */
    int max;                /* for normalization */
    int request;            /* request */
    int arr_rand[RAND_ARR_SIZE]; /* array of random integers */
    int total_in;           /* points inside from all workers */
    int total_out;          /* points outside from all workers */    
    double eps;             /* accuracy */
    double x;               /* random x-coordinate */
    double y;               /* random y-coordinate */
    double pi;              /* approximated value of pi */
    double err;             /* error compared to M_PI */

    /* init */
    numb_prcs = 0;
    curr_rank = 0;
    srvr_rank = 0;
    arr_ranks[0] = 0;    

    i = 0;
    iter = 0;
    in = 0;
    out = 0;
    done = 0;
    max = INT_MAX;
    request = 0;
    eps = 0.0;  
    x = 0.0;
    y = 0.0;
    pi = 0.0;
    err = 0.0;
    
    /* init MPI */
    MPI_Init(&argc, &argv);

    /* set the world communicator */
    wrld = MPI_COMM_WORLD;

    /* get number of processes */
    MPI_Comm_size(wrld, &numb_prcs);

    /* get the rank of the current process */
    MPI_Comm_rank(wrld, &curr_rank);

    /* set the rank of the random number server process */
    srvr_rank = numb_prcs - 1; /* the last process */

    /* process with rank 0 reads accuracy and broadcasts it */
    if (curr_rank == 0)
    {
        if (argc < 2) 
        {
	        fprintf(stderr, "Usage: %s accuracy\n", argv[0]);
	        MPI_Abort(wrld, 1);
	    }

        /* read formatted data from a string / from command-line argument */
        sscanf(argv[1], "%lf", &eps);
    }

    MPI_Bcast(
        &eps,                   /* data to be communicated */
        1,                      /* number of items */
        MPI_DOUBLE,             /* data type of communicated data */
        0,                      /* rank of the sender process */
        MPI_COMM_WORLD          /* communicator */
        );

    /* extract the group from the MPI_COMM_WORLD communicator */
    MPI_Comm_group(
        wrld,                   /* communicator */
        &wrld_group             /* group */
        );

    /* exclude the random numbers server process to create workers group */
    arr_ranks[0] = srvr_rank;   /* exclude only the server */
    MPI_Group_excl(
        wrld_group,             /* initial group */
        1,                      /* number of processes to exclude */
        arr_ranks,              /* array of ranks of processes to exclude */ 
        &wrkr_group             /* the new group */
        );

    /* create the new communicator */
    /* collective operations that do not involve the random number server */
    MPI_Comm_create(
        wrld,                   /* initial communicator */
        wrkr_group,             /* group for the new communicator */
        &wrkr                   /* the new communicator */
        );

    /* release the worker group */
    /* free reference but does not destroy the group */
    MPI_Group_free(&wrkr_group);

    /* the random number server */
    if (curr_rank == srvr_rank)
    {
        srandom(time(0));
        do
        {
            /* receive a request from a worker */
            MPI_Recv(
                &request,       /* address of receive buffer */
                1,              /* number of elements */
                MPI_INT,        /* data type of communicated data */
                MPI_ANY_SOURCE, /* rank of the source */
                RQUST,          /* message tag is request */
                wrld,           /* communicator */
                &status         /* status object */
                );

            if (request)
            {
                /* generate array of random integers */
                while (i < RAND_ARR_SIZE)
                {
                    arr_rand[i] = random();
                    if (arr_rand[i] <= INT_MAX)
                    {
                        i++;
                    }
                }

                /* send the array of random integers to the worker */
                MPI_Send(
                    arr_rand,           /* address of send buffer */
                    RAND_ARR_SIZE,      /* number of elements */
                    MPI_INT,            /* data type of communicated data */
                    status.MPI_SOURCE,  /* rank of the destination process */
                    REPLY,              /* message tag is replay */
                    wrld                /* communicator */
                    ); 
                
            }
        }
        while (request > 0);
    }
    /* a worker process */
    else
    {
        request = 1;
        done = 0;
        in = 0;
        out = 0;
        max = INT_MAX;

        /* send request to the random numbers server */
        MPI_Send(
            &request,           /* address of send buffer */
            1,                  /* number of elements */
            MPI_INT,            /* data type of communicated data */
            srvr_rank,          /* rank of the destination process */
            RQUST,              /* message tag is request */
            wrld                /* communicator */
            );

        iter = 0;
        while (!done)
        {
            iter++;
            request = 1;

            MPI_Recv(
                &arr_rand,          /* address of receive buffer */
                RAND_ARR_SIZE,      /* number of elements */
                MPI_INT,            /* data type of communicated data */
                srvr_rank,          /* rank of the source */
                REPLY,              /* message tag is replay */
                wrld,               /* communicator */
                MPI_STATUS_IGNORE   /* ignore the status */
                );
            /* count points inside and outside the circle */
            i = 0;
            while (i < RAND_ARR_SIZE)
            {
                x = (((double) arr_rand[i++]) / max) * 2 - 1;
                y = (((double) arr_rand[i++]) / max) * 2 - 1;
                (x * x + y * y < 1.0) ? in++ : out++;
            }

            /* sum points inside and outside from all workers */
            MPI_Allreduce(
                &in,                    /* address of sending buffer */
                &total_in,              /* address of receive buffer */
                1,                      /* number of elements */
                MPI_INT,                /* data type of communicated data */
                MPI_SUM,                /* operation applied */
                wrkr                    /* communicator */
                );
            MPI_Allreduce(
                &out,                   /* address of sending buffer */
                &total_out,             /* address of receive buffer */
                1,                      /* number of elements */
                MPI_INT,                /* data type of communicated data */
                MPI_SUM,                /* operation applied */
                wrkr                    /* communicator */
                );
            
            /* calculate pi approximation */
            pi = (4.0 * total_in) / (total_in + total_out);
            err = fabs(pi - M_PI);
            done = (err < eps) || (total_in + total_out > MAX_POINTS);
            request = (done) ? 0 : 1;

            if (curr_rank == 0)
            {
                /* printf("Approximated: %.20f Error: %.20f\n", pi, err); */

                MPI_Send(
                    &request,           /* address of send buffer */
                    1,                  /* number of elements */
                    MPI_INT,            /* data type of communicated data */
                    srvr_rank,          /* rank of the destination process */
                    RQUST,              /* message tag is request */
                    wrld                /* communicator */
                    );
            }
            else
            {
                if (request)
                {
                    MPI_Send(
                        &request,           /* address of send buffer */
                        1,                  /* number of elements */
                        MPI_INT,            /* data type of communicated data */
                        srvr_rank,          /* rank of the destination process */
                        RQUST,              /* message tag is request */
                        wrld                /* communicator */
                        );
                }
            }
        }

        /* release the worker communicator */
        /* the group does not cease to exist until both references to the group have been freed */
        MPI_Comm_free(&wrkr);
    }

    if (curr_rank == 0)
    {
        printf("Approximated: %.20f Error: %.20f\n", pi, err);
    }

    /* clear MPI resources */
    MPI_Finalize();

    return 0;
}/* main */
