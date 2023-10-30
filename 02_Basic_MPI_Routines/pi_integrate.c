#include <mpi.h>
#include <stdio.h>
#include <math.h>

double func(double x)
{
    return (4.0 / (1.0 + x * x));
}/* func */

int main(int argc, char* argv[])
{
    /* variables */
    int world_size;         /* number of processes */
    int prc_rank;           /* rank of the process */
    int numb_inter;         /* number of intervals */
    int more;               /* more iterations or stop */
    int i;                  /* loop counter */
    double h;               /* size of the sub-interval */
    double sum;             /* accumulate the sum */
    double x;               /* argument of the function to integrate */
    double prc_pi;          /* the approximation in a process */
    double pi;              /* the general approximation */

    /* init */
    world_size = 0;
    prc_rank = 0;
    numb_inter = 0;
    more = 1;
    i = 0;
    h = 0.0;
    sum = 0.0;
    x = 0.0;
    prc_pi = 0.0;
    pi = 0.0;

    /* init MPI */
    MPI_Init(&argc, &argv);

    /* get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /* get the rank of the process */
    MPI_Comm_rank(MPI_COMM_WORLD, &prc_rank);

    while (more)
    {
        /* manager process asks uses for the number of intervals */
        if (prc_rank == 0)
        {
            printf("Number of intervals, or 0 to stop: ");
            fflush(stdout);
            scanf("%d", &numb_inter);
        }

        /* manager process send numb_inter to all other processes */
        MPI_Bcast(
                &numb_inter,    /* data to be communicated */
                1,              /* number of items */
                MPI_INT,        /* data type of communicated data */
                0,              /* rank of the sender process */
                MPI_COMM_WORLD  /* communicator */
                );

        if (numb_inter > 0)
        {
            /* perform integration */
            h = 1.0 / (double) numb_inter;
            sum = 0.0;
            for (i = prc_rank + 1; i <= numb_inter; i += world_size)
            {
                x = h * ((double)i - 0.5);
                sum += func(x);
            }
            prc_pi = h * sum;

            /* add up individual processes */
            MPI_Reduce(
                    &prc_pi,        /* source address */
                    &pi,            /* result address */
                    1,              /* number of items */
                    MPI_DOUBLE,     /* data type of communicated data */
                    MPI_SUM,        /* performed operation is sum */
                    0,              /* rank of the process to place the result */
                    MPI_COMM_WORLD  /* communicator */
                    );
            if (prc_rank == 0)
            {
                printf("PI approximation: %.16f, error: %.16f\n", pi, fabs(pi - M_PI));
            }
        }
        else
        {
            /* stop the main loop in the process */
            more = 0;
        }
    }

    /* clear MPI resources */
    MPI_Finalize();

    return 0;
}/* main */
