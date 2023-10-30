#include <stdio.h>

int main()
{
    /* variables */
    int i;
    int vect_size;
    int weight;

    /* init */
    unsigned int vect[] = 
        {0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1};
    i = 0;
    weight = 0;
    vect_size = sizeof(vect) / sizeof(vect[0]);

    for (i = 0; i < vect_size; i++)
    {
        if (vect[i])
        {
            weight++;
        }
    }

    printf("Weight of the vector: %d\n", weight);

    return 0;
}//main
