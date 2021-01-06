/*FUNCTIONS FOR PACKAGE COLLPCM
 
 Author:	Jason Wyse,
 School of Computer Science and Statistics,
 Lloyd Institute,
 Trinity College,
 Dublin 2,
 Ireland.
 mailto: wyseja@tcd.ie
 
 Last modification of this code: Wed 16 Sep 2020 21:14:42 IST   */

#ifndef __COLLPCM_MIXMOD_H__
#define __COLLPCM_MIXMOD_H__

#include "COLLPCM_required_libs.h"
#include "COLLPCM_utils.h"

struct component
{
  int in_use; //indicator saying whether the component is in use or not
  int d; //dimension
  int n_g; //number of members of component
  double *sum; //gives the sum of the items
  double sum_squared_norm; 
  double log_prob; //quantity that can be updated to save computations
};

struct mix_mod
{
  //structure to hold mixture essentials
  
  int G; //number of groups
  int n; //number of items
  int d; //item dimension
  
  int G_max; //max number of groups
  
  double *y; //pointer to data
  int *z; //class membership
  int *whereis; //indexes of components
  struct component **components; //pointer to array of components
  
  // priors
  double *prior_mu; //prior mean for component means
  double xi2; //squared norm of prior_mu (precompute)
  double *log_prior_G; //this is the log prior on the number of groups
  
  // hyper-parameters
  double alpha; //dirichlet prior on weights (symmetric)
  double delta; //twice the prior shape on the component precisions
  double gamma; //twice the prior rate on the component precisions
  double kappa; //fold change component error precisions are by for the prior on the within component means*/
  double lambda; //prior mean of a Poisson on the no. of components
  
  // oft-used pre-computables
  
  // hyper-prior-parameters
  double shape_lambda; //twice shape parameter for the gamma hyperprior on lambda
  double rate_lambda; //twice rate parameter for the gamma hyperprior on lambda
  double shape_kappa; //twice shape parameter for the gamma hyperprior on kappa
  double rate_kappa; //twice rate parameter for the gamma hyperprior on kappa
  double shape_gamma; //twice shape parameter for the gamma hyperprior on gamma
  double rate_gamma; //twice parameter for the gamma hyperprior on gamma
  double precision_prior_mu; // precision on spherical prior for mu
    
  // mcmc control parameters
  int update_lambda; //logical whether lambda should be updated.
  int update_kappa; //logical whether kappa should be updated or not
  int update_gamma; //logical whether gamma should be updated after each sweep  
  int update_prior_mu; //logical whether to draw prior mean of component means again or not
  
  double *table_a; //this is a lookup table for vals of a (eject/absorb moves)
};

struct results
{
  // structure to hold runtime statistics
  // label moves
  int proposed_m1;
  int accepted_m1;
  int proposed_m2;
  int accepted_m2;
  int proposed_m3;
  int accepted_m3;
  // eject / absorb
  int proposed_eject;
  int accepted_eject;
  int proposed_absorb;
  int accepted_absorb;
};


#include "COLLPCM_component.h"
#include "COLLPCM_label_updates.h"
#include "COLLPCM_eject_absorb.h"
#include "COLLPCM_analysis.h"


struct mix_mod *COLLPCM_allocate_mixmod(int n, int d, int G_max, int G );

struct mix_mod *COLLPCM_clone_mixmod( struct mix_mod *mm );

void COLLPCM_free_mixmod(struct mix_mod *mm);

void COLLPCM_set_prior_hyperparameters(  struct mix_mod *mixmod , double gamma_0, double kappa_0);

struct results *COLLPCM_create_results();

#endif
