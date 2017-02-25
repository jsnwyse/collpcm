/*C functions to do inference for a collapsed Gaussian mixture model

	Author:	Jason Wyse, 
			School of Computer Science and Statistics,
			Trinity College Dublin,
			Dublin 2, Ireland.
			email: wyseja@tcd.ie
			
Last modified: Fri 14 Mar 2014 13:00:01 GMT  */ 



#ifndef _GAUSSIAN_MIXTURE_MODEL_H_
#define _GAUSSIAN_MIXTURE_MODEL_H_

#include "required_libs.h"

//#include "efficient_gamma.h"

#define TRUE 1
#define FALSE 0
#define log_2_pi 1.83787706640934533908
#define log_pi 1.1447298858494001639


/*definition for different prior types*/
#define RICHARDSON_AND_GREEN 0
#define NOBILE_AND_FEARNSIDE 1
#define CUSTOM 2
#define TRIONA 3
#define HYPERPRIOR_LAMBDA 4 /*probably won't use that much, but to include the */
#define TRIONABASIC 5

/*definition for initialization types*/
#define INITIALIZE_EM 0
#define INITIALIZE_SIMPLE 1


struct component
{
	int in_use; /*indicator saying whether the component is in use or not*/
	
	int n_g; /*number of members of component*/
	
	int d; //dimension of the vector
	
	//int *indexes_g; /*indexes of members of component in Y (row-wise)*/
	
	double *sum; /*sum of data in component*/
	
	double sum_squared_norm; /*sum of squared norm*/
	
	double log_prob; /*quantity that can be updated to save computations*/

};


struct mix_mod
/*structure to hold mixture essentials*/
{
	int G; /*total number of groups*/
	
	int n; /*total number of data*/
	
	int d; /*dimension of the data i.e. R^d*/
	
	int maxgroups;
	
	double **Y; /*the raw data stacked row-wise i.e. y_i^t*/
	
	double *y_uni; //univariate y
	
	int *z; /*group memberhsip indicator*/
	
	int *whereis; /*stores the index in components of a component*/

	struct component **components; /*pointer to array of pointers to components*/
	
	int vcov_type; /*form of vcov matrix to use for the model SPHERICAL (only)*/
	
	double *prior_mu; /*vector giving the prior mean of the component means*/
	
	double xi2; /*the squared norm of the prior mean on the component means*/
	
	/*other hyperparameters*/
	
	double alpha; /*alpha: dirichlet prior on weights symmetric*/
	
	double delta; /*delta: TWICE the prior shape on the component precisions*/
	
	double gamma; /*gamma: TWICE the prior rate on the component precisions*/
	
	double kappa; /*kappa: amount by which within component error precisions are 
							by for the prior on the within component means*/
							
	double lambda; /*this is the prior mean of a Poisson on the no. of components*/
	
	double shape_lambda;
	
	double rate_lambda;
	
	int update_lambda; /*logical indicating whether lambda should be updated.*/
	
	double *log_prior_G; /*this is the prior for the number of groups... poisson(1)*/
							
	double *table_a; /*this is a lookup table for the values of a when ejecting or
								absorbing components*/
								
	int update_kappa; /*logical indicating whether kappa should be updated or not*/
	
	double shape_kappa; /*TWICE shape parameter for the gamma hyperprior on kappa*/
	
	double rate_kappa; /*TWICE rate parameter for the gamma hyperprior on kappa*/
	
	int update_gamma; /*logical indicating whether gamma should be updated after each sweep*/
	
	double shape_gamma; /*TWICE shape parameter for the gamma hyperprior on gamma*/
	
	double rate_gamma; /*TWICE parameter for the gamma hyperprior on gamma*/
	
	int update_prior_mu; /*logical value indicating whether to draw prior mean of component means again or not*/
	
	double precision_prior_mu;
	
	FILE *fp_log; /*file pointer for the debugger log*/

};


struct results
/*a structure to store all of the results of the analysis*/
{
	
	int niterations; /*number of iterations*/
	
	int nburnin; /*number of burnin iterations*/
		
	int *ngroups; /*gives number of groups at each iteration*/
	
	int **memberships; /*gives membership at each iteration*/
	
	int *MAP_memberships; /*gives membership with maximum prob for each i=1,\dots,n*/
	
	int proposed_m1;
	
	int accepted_m1;
	
	int proposed_m2;
	
	int accepted_m2;
	
	int proposed_m3;
	
	int accepted_m3;
	
	int proposed_eject;
	
	int accepted_eject;
	
	int proposed_absorb;
	
	int accepted_absorb;
	
};

/*functions*/

void mixmod_warning(int warning);

struct mix_mod *mixmod_create( int datasize, int datadimension, int maxgroups,  int initgroups );

void mixmod_destroy( struct mix_mod *mm );

struct component *component_create( int d );

void component_destroy( struct component *c );

void component_refresh( struct component *c );

void copy_component( struct component *component_original, struct component *component_target );

void component_add_to_component( struct component *comp, double *x, int sgn ) ;

void component_add_to_component_uni( struct component *comp, double x, int sgn ) ;

void allocate_results(struct results *results,int iterations,int burn_in,int len);

void free_results(struct results *results,int iterations,int burn_in);

int initialize_simple(struct mix_mod *mixmod,int numgroups);

double get_max(double *x,int len);

double get_min(double *x,int len);

int get_imax(int *x,int len);

int sample_discrete(double *weights, int len );

int get_ind_max(double *x,int len);

void set_prior_hyperparameters_x(  struct mix_mod *mixmod , double gamma_0, double kappa_0);

void get_hyperparameters_x( double *h, double gamma_0, double kappa_0 );

void set_prior_hyperparameters(struct mix_mod *mixmod,int type);

void set_prior_on_number_of_components(struct mix_mod *mixmod,int type);

int update_allocations_with_gibbs(struct mix_mod *mixmod);

void update_allocations_with_metropolis_move_1(struct mix_mod *mixmod,int *accepted,int *proposed);

void update_allocations_with_metropolis_move_2(struct mix_mod *mixmod,int *accepted,int *proposed);

void update_allocations_with_metropolis_move_3(struct mix_mod *mixmod,int *accepted,int *proposed);

void update_allocations_with_ejection_move(struct mix_mod *mixmod,int *accepted,int *proposed,double pr_ej_G,double pr_ej_Gp1);

void update_allocations_with_absorb_move(struct mix_mod *mixmod,int *accepted,int *proposed,double pr_ej_G,double pr_ej_Gm1);

//int update_allocations_with_split_move(struct mix_mod *mixmod,int *accepted,int *proposed,double pr_ej_G,double pr_ej_Gp1);

//int update_allocations_with_combine_move(struct mix_mod *mixmod,int *accepted,int *proposed,double pr_ej_G,double pr_ej_Gm1);

void update_hyperparameters(struct mix_mod *mixmod);

void GMM_recompute_marginal_likelihood_component( int idx, struct mix_mod *mm );

void GMM_recompute_marginal_likelihood_component_0( struct component *cmp, struct mix_mod *mm );

double GMM_return_marginal_likelihood_component( struct component *cmp, struct mix_mod *mm );

double  GMM_compute_marginal_likelihood_with_inclusion_in_component( double  *x, struct component *comp, struct mix_mod *mm );

double  GMM_compute_marginal_likelihood_with_inclusion_in_component_uni( double  x, struct component *comp, struct mix_mod *mm );

//double integrated_likelihood_component(struct component *component, double scalem, double shape,double rate ,int d);

//ouble integrated_gradient_calc(struct component *component, double scalem, double shape,double rate ,int d);

double log_normalizing_constant_model( int G, struct mix_mod *mm );

void do_mixmod_analysis_one_sweep(struct results *pres,struct mix_mod *mixmod,int fix_G, int iter);

//int write_out_results(struct results *results,int N,int datasize,FILE *fp1,FILE *fp2);

void random_ranshuffle( int *a, int n );

#endif

