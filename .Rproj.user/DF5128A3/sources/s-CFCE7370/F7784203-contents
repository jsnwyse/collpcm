/*FUNCTIONS FOR PACKAGE COLLPCM
 
 Author:	Jason Wyse,
 School of Computer Science and Statistics,
 Lloyd Institute,
 Trinity College,
 Dublin 2,
 Ireland.
 mailto: wyseja@tcd.ie
 
 Last modification of this code: Wed 16 Sep 2020 21:14:42 IST   */



#include "COLLPCM_component.h"

struct component * COLLPCM_create_component( struct mix_mod *mixmod )
{
  struct component *comp = (struct component *)malloc(sizeof(struct component));
  COLLPCM_allocate_component( comp, mixmod );
  return( comp );
}

void COLLPCM_allocate_component(struct component *component,struct mix_mod *mixmod)
{
  //allocate memory for a component
  int j, d=mixmod->d;
  
  component->d = d;
  component->n_g = 0;
  component->sum = (double *)calloc(d,sizeof(double));
  component->sum_squared_norm = 0.0;
  component->log_prob = -DBL_MAX;

  return;
}

void COLLPCM_free_component( struct component *component )
{
  free( component->sum );
  free( component );
  return;
}

void COLLPCM_copy_component(struct component *component_original,struct component *component_target )
{
  // copy a component
  
  int i, j, d = component_original->d;
  double *sum_t = component_target->sum, *sum_o = component_original->sum ;
  
  component_target->n_g = component_original->n_g;
  component_target->sum_squared_norm = component_original->sum_squared_norm;
  for( i=0; i<d; i++ ) sum_t[i] = sum_o[i] ;
  
  component_target->log_prob = component_original->log_prob;
  
  return;
}


void COLLPCM_add_to_component( struct component *component, double *x, int sgn )
{
  int j, d = component->d;
  
  component->n_g += sgn ;
  
  for( j=0; j<d; j++ ) component->sum[j] += sgn * x[j] ;
  
  component->sum_squared_norm += sgn * ( x[j] * x[j] );
  
  return;
}

void COLLPCM_add_components( struct component *component_original, struct component *component_target )
{
  // add two components together- result is returned  in target
  int j, d = component_original->d ; 
  double *sum_t = component_target->sum, *sum_o = component_original->sum;
  
  component_target->n_g += component_original->n_g;
  
  for( j=0; j<d; j++ ) sum_t[j] += sum_o[j];
  
  component_target->sum_squared_norm += component_original->sum_squared_norm;
  
  return;
}


// to do after

/*double COLLPCM_compute_log_data_probability_with_inclusion_in_component(int *x,struct component *component,struct mix_mod *mixmod)
{
  
  int i, j, c, k=mixmod->component_compute, d=mixmod->d, *vind=mixmod->varindicator, *ncat=mixmod->ncat, *N;
  double log_prob = 0., s, ind, *beta;
  
  if( mixmod->prior_type > 0 )
  {
    for( j=0; j<d; j++ )
    {
      if( vind[j] )
      {
        s =  0.; 
        beta = mixmod->beta_prior[k][j] ;
        N = component->N[j];
        for( c=0; c<ncat[j]; c++ )
        { 
          ind = (x[j] == c) ? 1.0 : 0.0 ;
          log_prob += lgamma( N[c] + ind + beta[c] ); // - lgamma( beta[c] )
          s += beta[c] ;
        } 
        log_prob += mixmod->lg_sum_beta[k][j] - mixmod->lg_beta_sum[k][j] - lgamma( component->n_g + 1. + s ) ; 
      }
    }
  }
  
  return(log_prob);
}*/


double COLLPCM_compute_log_marginal_likelihood_with_inclusion_in_component( double  *x, struct component *component, struct mix_mod *mixmod )
{
  
  COLLPCM_add_to_component( component, x, 1 );
  
  double lp = COLLPCM_return_log_marginal_likelihood_component( component, mixmod );
  
  COLLPCM_add_to_component( component, x, -1 );
  
  return( lp );
}


/*double COLLPCM_compute_log_data_probability_component(struct component *component,struct mix_mod *mixmod)
{
  
  int i, j, c, k=mixmod->component_compute, d=mixmod->d, *ncat=mixmod->ncat, *vind=mixmod->varindicator, *N;
  double log_prob = 0., s, *beta;
  
  if( mixmod->prior_type ==  0 )
  {
    
    for( j=0; j<d; j++ )
    {
      
      if( vind[j] )
      {
        
        log_prob += lgamma( ncat[j] * mixmod->beta ) - ncat[j]*lgamma(mixmod->beta) - lgamma(component->n_g + ncat[j]*mixmod->beta );
        N = component->N[j];
        for( i=0; i<ncat[j]; i++ ) log_prob += lgamma(N[i] + mixmod->beta);
      }
      
    }
    
  }
  
  if( mixmod->prior_type > 0 )
  {
    for( j=0; j<d; j++ )
    {
      if( vind[j] )
      { 
        s = 0.;
        beta = mixmod->beta_prior[k][j];
        N = component->N[j];
        for( c=0; c<ncat[j]; c++ )
        { 
          log_prob += lgamma( N[c] + beta[c] ) - lgamma( beta[c] ) ;
          s += beta[c] ;
        } 
        log_prob += lgamma(s) - lgamma( component->n_g + 1.0 + s ) ; 
      } 
    }
  }
  
  return(log_prob);	
}*/

double COLLPCM_return_log_marginal_likelihood_component( struct component *component, struct mix_mod *mixmod )
{
  
  int i, n = component->n_g, d = mixmod->d ;
  double lp, a, sq_norm=0.0, delta=mixmod->delta, kappa=mixmod->kappa, 
    *sum = component->sum,
    *mu = mixmod->prior_mu ;
    
    lp = lgamma( n +  mixmod->alpha ) + lgamma( 0.5*( n * d + delta ) ) - 0.5* d * log( n + kappa ) ;
    
    for( i=0; i<d; i++ )
    {
      a = sum[i] + kappa * mu[i] ;
      sq_norm += a * a ;
    }
    
    lp -= .5 * ( n * d + delta ) * log( component->sum_squared_norm - sq_norm/( n + kappa ) + kappa * mixmod->xi2 + mixmod->gamma  );
    
    return( lp ) ; 
    
}

void COLLPCM_recompute_log_marginal_likelihood_component(struct component *component,struct mix_mod *mixmod)
{
  
  double log_prob = COLLPCM_return_log_marginal_likelihood_component( component, mixmod );
  
  component->log_prob = log_prob;
  
  return;
}


void COLLPCM_initialize_simple( struct mix_mod *mixmod ,int G )
{
  // do a simple initialisation
  
  int i, o, j, k, d=mixmod->d, n=mixmod->n;
  double *y;
  
  int *order = (int *)calloc( n , sizeof(int) );
  for( i=0; i<n; i++ ) order[i] = i;
  COLLPCM_random_ranshuffle( order, n );
  
  struct component *comp;
  
  int m = n/G; //gives the number of segments
  for( k=0; k<G-1; k++ )
  {
    comp = mixmod->components[k] ;
    comp->in_use = TRUE;
    
    //cycle through appropriate data
    for(i=k*m;i<(k+1)*m;i++)
    {
      o = order[i];
      mixmod->z[ o ] = k ;
      y = mixmod->y + o * d ;
      COLLPCM_add_to_component( comp, y, 1 );
    }
  }
  
  //special case for last group
  comp = mixmod->components[G-1];
  comp->in_use = TRUE;

  for( i=(G-1)*m; i<n; i++ )
  {
    o = order[i];
    mixmod->z[o] = G-1;
    y = mixmod->y + o * d ;
    COLLPCM_add_to_component( comp, y, 1 );
  }
  
  //compute the log_prob for each of the components
  for(k=0;k<G;k++)
  {
    COLLPCM_recompute_log_marginal_likelihood_component( mixmod->components[k], mixmod );
    mixmod->whereis[k] = k;
  }
  
  free( order );
  
  return;
}

