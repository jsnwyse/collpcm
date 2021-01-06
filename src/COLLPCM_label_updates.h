/*FUNCTIONS FOR PACKAGE COLLPCM
 
 Author:	Jason Wyse,
 School of Computer Science and Statistics,
 Lloyd Institute,
 Trinity College,
 Dublin 2,
 Ireland.
 mailto: wyseja@tcd.ie
 
 Last modification of this code: Wed 16 Sep 2020 21:14:42 IST   */

#ifndef __COLLPCM_LABEL_UPDATES_H__
#define __COLLPCM_LABEL_UPDATES_H__

#include "COLLPCM_mixmod.h"

int COLLPCM_update_allocations_with_gibbs(struct mix_mod *mixmod) ;

/*double BLCA_get_log_relative_probability_for_gibbs(int *x,struct component *component_k,struct component *component_g,struct mix_mod *mixmod) ;

int BLCA_update_allocations_with_metropolis_move_1(struct mix_mod *mixmod,int *accepted,int *proposed) ;

int BLCA_update_allocations_with_metropolis_move_2(struct mix_mod *mixmod,int *accepted,int *proposed) ;

int BLCA_update_allocations_with_metropolis_move_3(struct mix_mod *mixmod,int *accepted,int *proposed) ;*/

#endif
