/*FUNCTIONS FOR PACKAGE COLLPCM
 
 Author:	Jason Wyse,
 School of Computer Science and Statistics,
 Lloyd Institute,
 Trinity College,
 Dublin 2,
 Ireland.
 mailto: wyseja@tcd.ie
 
 Last modification of this code: Wed 16 Sep 2020 21:14:42 IST   */

#ifndef __COLLPCM_COMPONENT_H__
#define __COLLPCM_COMPONENT_H__

#include "COLLPCM_required_libs.h"
#include "COLLPCM_mixmod.h"

struct component * COLLPCM_create_component( struct mix_mod *mixmod ) ;

void COLLPCM_allocate_component(struct component *component,struct mix_mod *mixmod) ;

void COLLPCM_free_component(struct component *component ) ;

void COLLPCM_copy_component(struct component *component_original,struct component *component_target ) ;

void COLLPCM_add_to_component( struct component *component, double *x, int sgn ) ;

void COLLPCM_add_components( struct component *component_original, struct component *component_target ) ;

// still to do these

//double COLLPCM_compute_log_data_probability_with_inclusion_in_component(int *x,struct component *component,struct mix_mod *mixmod) ;

double COLLPCM_compute_log_marginal_likelihood_with_inclusion_in_component(int *x, struct component *component, struct mix_mod *mixmod ) ; // done

//double COLLPCM_compute_log_data_probability_component(struct component *component,struct mix_mod *mixmod) ;

double COLLPCM_return_log_marginal_likelihood_component( struct component *component, struct mix_mod *mixmod ); // done

void COLLPCM_recompute_log_marginal_likelihood_component(struct component *component,struct mix_mod *mixmod) ; // done

void COLLPCM_initialize_simple( struct mix_mod *mixmod ,int G );

#endif
