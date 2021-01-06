/*FUNCTIONS FOR PACKAGE COLLPCM
 
 Author:	Jason Wyse,
 School of Computer Science and Statistics,
 Lloyd Institute,
 Trinity College,
 Dublin 2,
 Ireland.
 mailto: wyseja@tcd.ie
 
 Last modification of this code: Wed 16 Sep 2020 21:14:42 IST   */


#include "COLLPCM_mixmod.h"

struct mix_mod *COLLPCM_allocate_mixmod(int n, int d, int G_max, int G )
  /*this function allocates and returns a pointer to a mixmod structure...*/
{
  
  int i,j,k;

  struct mix_mod *mm = (struct mix_mod *)malloc(sizeof(struct mix_mod));
  
  mm->G_max = G_max;
  mm->G = G;
  mm->n = n;
  mm->d = d;
  
 // mm->z = (int *)calloc( n,sizeof(int) );
  mm->whereis = (int *)calloc( G_max, sizeof(int) );
  mm->components = (struct component **)malloc( sizeof(struct component *) * G_max );
  
  for( i=0; i<G_max; i++ )
    mm->components[i] = (struct component *)malloc(sizeof(struct component));
  
  
  // priors
  mm->prior_mu = (double *)calloc( d, sizeof(double) );
  
  // need to populate the hyperparameters
  
  return(mm);
  
}

struct mix_mod *COLLPCM_clone_mixmod( struct mix_mod *mm )
{
  
  int k;
  double *hparams = (double *)calloc( 2, sizeof(double) );
  hparams[0] = mm->alpha;

  struct mix_mod *mm_clone; 
  
  mm_clone = COLLPCM_allocate_mixmod( mm->n, mm->d, mm->G, mm->G );
  
  for( k=0; k<mm->G; k++ )
    COLLPCM_copy_component( mm->components[k], mm_clone->components[k] );
  
  free( hparams );
  
  return( mm_clone );
}


void COLLPCM_free_mixmod(struct mix_mod *mm)
  /*frees the memory used by mixmod object*/
{
  int n = mm->n,d = mm->d, G = mm->G, G_max = mm->G_max, k;
  
  //free up components
  for( k=0; k<G_max; k++ )
  {
    COLLPCM_free_component(mm->components[k]);
    //free(mm->components[k]);
  }
  free(mm->components);
  
  //free wehreis
  free(mm->whereis);

  return;
}

void COLLPCM_set_prior_hyperparameters(  struct mix_mod *mixmod , double gamma_0, double kappa_0)
{
  
  double c = 4., d = 4.;
  
  mixmod->shape_gamma = 2. * c * c ;
  mixmod->rate_gamma = mixmod->shape_gamma / gamma_0 ;
  
  mixmod->shape_kappa = 2. * d * d ;
  mixmod->rate_kappa = mixmod->shape_kappa / kappa_0 ;
  
}


struct results *COLLPCM_create_results()
{
  struct results *r= (struct results *)malloc(sizeof(struct results));
  r->proposed_m1 = 0;
  r->accepted_m1 = 0;
  r->proposed_m2 = 0;
  r->accepted_m2 = 0;
  r->proposed_m3 = 0;
  r->accepted_m3 = 0;
  r->proposed_eject = 0;
  r->accepted_eject = 0;
  r->proposed_absorb = 0;
  r->accepted_absorb = 0;
  return(r);
}



