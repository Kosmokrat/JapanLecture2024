/**
 *    Written for my Japan Lectures at Yokosuka, October 2024
 *    Computational Materials Science with Applications in Molecular Dynamics
 *    -------------------------------------------------------------------------
 *    Author: Prof. Dr. rer. nat. habil. Martin O. Steinhauser
 *
 */

#include  <stdio.h>
#include  <stdlib.h>
#include  <time.h>

#define NUM_ARRAY_ENTRIES 10 

void insert_sort(int n, int z[])
{
    int  i, j, t, k;
    
    for (i=1; i<n; i++) { /* insert i.th  element */
       for (j=i;  j>0 && z[j] < z[j-1]; j--) {
          t      = z[j];  /*-- Exchange  z[j] and z[j-1] */
          z[j]   = z[j-1];
          z[j-1] = t;
       }
       /*------------ Output -------------------*/
       printf("%2d. Run: ", i);
       for (k=0; k < NUM_ARRAY_ENTRIES; k++)
          printf("%3d", z[k]);
       printf("\n");
       /*-------------------------------------------------------------------*/
    }
}

int  main(void)
{
   int   i, k, numbers[NUM_ARRAY_ENTRIES];
   
   clock_t start, stop;
   start = clock();

   srand(time(NULL));

   for (i=0 ; i<NUM_ARRAY_ENTRIES ; i++)
      numbers[i] = rand()%100;

   printf("---- Before Insertion Sort -----------------\n");
   for (k=0; k<NUM_ARRAY_ENTRIES; k++)
      printf("%3d", numbers[k]);
   printf("\n\n");

   insert_sort(NUM_ARRAY_ENTRIES, numbers);

   printf("---- After Insertion Sort ----------------\n");
   for (k=0; k<NUM_ARRAY_ENTRIES; k++)
      printf("%3d", numbers[k]);
   printf("\n\n");

   stop = clock();
   printf("\n===========================================\n");
   printf("Time spent with sorting: %.16f seconds.\n", ((float)stop-start)/CLOCKS_PER_SEC);
   printf("\n===========================================\n");
   return(0);
}
