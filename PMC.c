/* Japan Lecture, Computational Science with Applications in Molecular Dynamics;
   October 2024
   Written by Martin O. Steinhauser, Oct. 2024
*/


/* IMPORTANT: 
   For velocity Verlet, we need an additional array Fold[3] which has to be declared within the structure for the Particle !
*/


/* IMPORTANT: The output function, that writes the particle data into a file is most conveniently written
   as a function computeResults() within the function timeIntegration() - you can just do simple
   print-statements into a file after each timestep.
*/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>


/* ============================================================== */
#define DIM 3


/** 
    \def VERSION
    A global variable indicating the current version of the application.
*/
#define VERSION 1.0


/** \def sqr(x)
    This is the definition of the squared value of a number.
*/
#define sqr(x) ((x) * (x));


/**
   \struct Particle
*/
typedef struct {
  double m;           /* mass */
  double x[DIM];      /* position */
  double v[DIM];      /* velocity */
  double f[DIM];      /* force */
  double fOld[DIM];   /* Old force - needed for Verlet-velocity */
} Particle;

/* ============================================================== */



void computeStatistics(Particle *p, int N, double time){

  int i, d;
  double velocity;
  double energy = 0.;

  for (i = 0; i < N; i++){

    velocity = 0.;

    for (d = 0; d < DIM; d++)
      velocity += sqr(p[i].v[d]);
    
    energy += .5 * p[i].m * velocity;
  }
  /* If you want: insert a function here, that prints out the kinetic energy 
     at this timestep time */
  /* printf("Check: E_kin = %f\n",energy); */
}


void updateX(Particle *p, double deltaTime){
  int i, d;

  double a = deltaTime * .5 / p -> m;

  for ( d = 0; d < DIM; d++) {
    p -> x[d] += deltaTime * (p -> v[d] + a * p -> f[d]);  /* According to Eq. (A) */
    p -> fOld[d] = p -> f[d];
  }
}


void updateV(Particle *p, double timeDelta){
  int i, d;
  double a = timeDelta * .5 / p -> m;
  for ( d = 0; d < DIM; d++)
    p -> v[d] += a * (p -> f[d] + p -> fOld[d]);  /* According to Eq. (B) */
}


void computeX(Particle *p, int N, double deltaTime){
  int i;

  for (i = 0; i < N; i++)
    updateX(&p[i], deltaTime);
}



void computeV(Particle *p, int N, double deltaTime){
  int i;

  for (i = 0; i < N; i++)
    updateV(&p[i], deltaTime);
}




void forceCalculate(Particle *i, Particle *j){
  int d;
  double distance = 0.0;
  double force = 0.0;

  for (d = 0; d < DIM; d++)
    distance += sqr(j -> x[d] - i ->x[d]);

  force = i -> m * j -> m / (sqrt(distance) * distance);

  for (d = 0; d < DIM; d++)
    i -> f[d] += force * (j -> x[d] - i -> x[d]);
}




/*================================================================================ */
/**
   void computeForces(Particle *p, int N)

   \brief Computation of the new forces based after \em one timestep 
   using the Verlet-velocity method with the all particle approach

   \param *p A pointer to the particle information
   \param N The total number of particles

   In the basis version of the code the force is calculatied according to 
   \f[ 
   \vec{F}_i=-\nabla_{\vec{x}_i} V(\vec{x}_1,...,\vec{x}_N) = -\sum_{j=1,j\neq i}^{N}
   \nabla_{\vec{x}_i}U(r_{ij})=\sum_{j=1,j\neq i}^{N}\vec{F}_{ij} 
   \f], i.e. the two forces \f$\vec{F}_{ij}\f$ and  \f$\vec{F}_{ji}\f$ are calculated 
   separately.\n\n
   Alternatively, in an improved version, Newton's third law is taken into account, i.e
   the force  \f$\vec{F}_{ji}\f$ need not be calculated if \f$\vec{F}_{ij}\f$ already
   has been evaluated. It is \f$\vec{F}_{ji} = -\vec{F}_{ij}\f$. That way one can save 
   half of the computations of the forces. 

*/
void computeForces(Particle *p, int N)
{
  int d,i,j;

  for (i = 0; i < N; i++) // Comment 
    for (d = 0; d < DIM; d++)
      p[i].f[d] = 0.0;
  /* All forces are now set to zero */

  for (i = 0; i < N; i++)
    for (j = 0; j < N ; j++)
      if (i != j) forceCalculate(&p[i], &p[j]);
}
/*-------------------------------------------------------------------------------- */


void outputResults(Particle *p, int N, double time){

  static int once = 1;

  if (once){
    printf("# Sun:  time X  Y  Z\n");
    once = 0;
  }

  printf(" %f %f %f %f\n",time, p[4].x[0], p[4].x[1], p[4].x[2]);
}


void timeIntegration(double time, 
		     double deltaTime, 
		     double timeEnd, 
		     Particle *p, 
		     int N){
  computeForces(p, N);
   
  while (time < timeEnd) {
    time += deltaTime;
     
    computeX (p, N, deltaTime);
     
    computeForces (p, N);
      
    computeV (p, N, deltaTime);
     
    computeStatistics (p, N, time);

    outputResults(p, N, time);
  }
}
   

/**
   \fn void  initAllData(Particle *p)
   
   \brief This function initializes all particle attributes in the structure Particle. For 2D and 3D this has to be
   done explicitly here in the code.

   The case study of planetary motion. This study works in 2D AND 3D. The particle information is always allocated as 3D, but only in
   the 2D case the 3D information is ignored in the corresponding loops.
   
   \param *p The pointer to all particle properties.
   \param *params The pointer to all parameters.
*/


void initAllData(Particle *p)
{
  /* The case study of planetary motion. This study works in 2D AND 3D. 
     The particle information is always allocated as 3D, but only in
     the 2D case the 3D information is ignored in the corresponding loops. */
 
  /* Sun */
  p[0].m = 1.; p[0].x[0] = p[0].x[1] = p[0].x[2] = 0.;   p[0].v[0] = p[0].v[1] = p[0].v[2] = 0.;
      
  /* Earth */
  p[1].m = .000003; p[1].x[0] = 0; p[1].x[1] = 1.; p[1].x[2] = 0.;  p[1].v[0] = -1.0; p[1].v[1] = 0.; p[1].v[2] = 0.;
      
  /* Jupiter */
  p[2].m = .000955;  p[2].x[0] = 0.; p[2].x[1] = 5.36; p[2].x[2] = 0.; p[2].v[0] = -0.425; p[2].v[1] = p[2].v[2] = 0.;
      
  /* Comet */
  p[3].m = .00000000000001;   p[3].x[0] = 34.75; p[3].x[1] = p[3].x[2] = 0.;   p[3].v[0] = 0; p[3].v[1] = 0.0296; p[3].v[2] = 0.;

  /* Planet2 */
  p[4].m = 0.00005; p[4].x[0] = 25.; p[4].x[1] = 0.; p[4].x[2] = 0.;   p[4].v[0] = 0.; p[4].v[1] = 0.02; p[4].v[2] = 0.02;

  /* BigPlanet */
  p[5].m = 0.1; p[5].x[0] = 0.; p[5].x[1] = 17.; p[5].x[2] = 0.;   p[5].v[0] = -0.25; p[5].v[1] = 0.; p[5].v[2] = 0.;
}



/* And Finally, the main function */

int main (void){
  /* REMARK: Assign here the values provided in Problem 2 of Problem Set 5 */
  int    N         = 6;
  double timeDelta = 0.015;
  double timeEnd   = 468.5;
   
  Particle *p = (Particle*) malloc(N * sizeof (*p));

  /* THE FOLLOWING FUNCTION HAS TO BE WRITTEN BY YOU!! */
  /* Here, you can hard-code all the initial conditions
     provided in ProblemSet5, i.e. you fill here the members 
     of the 4 different particles which are of type 
     struct Particle, e.g. by using commands like
     p[0].m = 1.;
     p[0].x[0] = p[0].x[1] = p[0].x[2] = 0.;
     p[1].m = .000003; 
     and so on.
  */
  initAllData(p);
   
  timeIntegration(0, timeDelta, timeEnd, p, N);
  free(p);
  return (0);  
}
