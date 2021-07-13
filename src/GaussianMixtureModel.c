/*C functions to do inference for a collapsed Gaussian mixture model

	Author:	Jason Wyse, 
			School of Computer Science and Statistics,
			Trinity College Dublin,
			Dublin 2, Ireland.
			email: wyseja@tcd.ie
			
Last modified: Fri 14 Mar 2014 13:00:01 GMT  */ 


#include "GaussianMixtureModel.h"

#define debug TRUE

void mixmod_warning(int warning)
/*warning handler*/
{
	#define NumWarningMsgs 4

  /*static const char *warningmsg[NumWarningMsgs] = {"The EM algorithm reached the maximum number of steps when initializing. Only a lower tolerance may be obtainable.",
  	"Computation of the marginal likelihood for a cluster returned a nan",
  	"There is a counting bug in metropolis move 1","Triona's initialization type cannot use the range as in conventional as only distances between points stored"};*/

 // printf("\n***Warning*** \t Reason: %s\n\n",warningmsg[warning]);
  return;

}


struct mix_mod *mixmod_create( int datasize, int datadimension, int maxgroups,  int initgroups )
{

	int i;

	struct mix_mod *mm = (struct mix_mod *)malloc( sizeof(struct mix_mod) );
	
	mm->G = initgroups;
	mm->n = datasize;
	mm->d = datadimension;
	mm->maxgroups = maxgroups;
	
	if( mm->d == 1 )
	{
		mm->y_uni = calloc( datasize, sizeof(double) );
	}
	else
	{
		mm->Y = calloc(datasize,sizeof(double *));
		for(i=0;i<datasize;i++){
			mm->Y[i] = calloc(datadimension,sizeof(double));
		}
	}

	mm->z = calloc(datasize,sizeof(int));
	
	/*allocate this memory-- only initialize what is needed for initial conditions*/
	mm->components = (struct component **)malloc(sizeof(struct component *)*maxgroups);

	for(i=0;i<maxgroups;i++){
		mm->components[i] = component_create( datadimension );
	}


	/*allocate whereis*/
	
	mm->whereis = calloc(maxgroups,sizeof(int));
	for(i=0;i<maxgroups;i++)
		mm->whereis[i] = -1;

	if( mm->d > 1 )
		mm->prior_mu = calloc( mm->d ,sizeof(double));
	else
		mm->prior_mu = calloc( 2, sizeof(double) );
	
	mm->log_prior_G = calloc(maxgroups+1,sizeof(double));

	return( mm );

}


void mixmod_destroy( struct mix_mod *mm )
{
	int n = mm->n,i,k;
	
	//free up components
	for(k=0;k<mm->maxgroups;k++){
		component_destroy( mm->components[k] );
	}
	free(mm->components);
	
	//free wehreis
	free(mm->whereis);
	
	//free data
	if( mm->d > 1 )
	{
		for(i=0;i<n;i++){
			free(mm->Y[i]);
		}
		free(mm->Y);
	}else{
		free( mm->y_uni );
	}
	
	
	//free others
	free(mm->z);
	free(mm->prior_mu);	
	
	free(mm->log_prior_G);
	
	free( mm );
	
	return;

}


struct component *component_create( int d )
{
	/*allocate memory for a component*/
	
	struct component *c = (struct component *)malloc( sizeof( struct component) );
	c->in_use = FALSE ;
	c->n_g = 0;
	c->d = d ;
	c->sum_squared_norm = 0.;
	if( d == 1 ) 
		c->sum = calloc( 2, sizeof(double) );
	else
		c->sum = calloc( d, sizeof(double) );
		
	c->log_prob = -DBL_MAX;
	return( c );
}

void component_destroy( struct component *c )
{
	free( c->sum );
	free( c );
}

void component_refresh( struct component *c )
{
	int d = c->d;
	struct component *n = component_create( d );
	copy_component( n, c );
	component_destroy( n );
	return;
}


void component_add_to_component( struct component *comp, double *x, int sgn )
{
	int k;
	
	comp->n_g += sgn;
	
	for( k=0; k<comp->d; k++ )
	{
	 	comp->sum[k] += sgn * x[k] ;
		comp->sum_squared_norm += sgn * x[k] * x[k] ;
	}
	
}

void component_add_to_component_uni( struct component *comp, double x, int sgn )
{
	int k;
	
	comp->n_g += sgn;
	
	for( k=0; k<comp->d; k++ )
	{
	 	comp->sum[0] += sgn * x ;
		comp->sum_squared_norm += sgn * x * x ;
	}
	
}

void copy_component( struct component *component_original, struct component *component_target )
/*copy the contents of the first argument into the second component argument*/
{

	int i, d = component_original->d;
	
	component_target->in_use = component_original->in_use;
	component_target->n_g = component_original->n_g;
	for(i=0;i<d;i++){
		component_target->sum[i] = component_original->sum[i];
	}
	component_target->sum_squared_norm = component_original->sum_squared_norm;
	component_target->log_prob = component_original->log_prob;
	
	return;

}


void allocate_results(struct results *results,int iterations,int burn_in,int len)
/*allocates space to store post burn-in iterations*/
{
	
	int N = iterations-burn_in;
	
	results->ngroups = (int *)calloc(N,sizeof(int));
	
	results->niterations = iterations;
	results->nburnin = burn_in;
	
	/*results->memberships = calloc(N,sizeof(int *));
	for(i=0;i<N;i++){
		results->memberships[i] = calloc(len,sizeof(int));
	}*/
	
	results->MAP_memberships = (int *)calloc(len,sizeof(int));
	
	results->proposed_m1 = 0;
	results->accepted_m1 = 0;
	results->proposed_m2 = 0;
	results->accepted_m2 = 0;
	results->proposed_m3 = 0;
	results->accepted_m3 = 0;
	results->proposed_eject = 0;
	results->accepted_eject = 0;
	results->proposed_absorb = 0;
	results->accepted_absorb = 0;
	
	return;
}

void free_results(struct results *results,int iterations,int burn_in)
{
	
	free(results->ngroups);
	
	/*for(i=0;i<N;i++){
		free(results->memberships[i]);
	}
	free(results->memberships);*/
	
	free(results->MAP_memberships);
	
	return; 
}


int initialize_simple(struct mix_mod *mixmod,int numgroups)
//gives a very simple initialization of the model by just dumping the
//	first n/G observations in the first group, second n/G in second group...
{
	int i, o, j,k,G=numgroups,d=mixmod->d,n=mixmod->n;
	
	int *order = (int *)calloc( n , sizeof(int) );
	for( i=0; i<n; i++ ) order[i] = i;
	random_ranshuffle( order, n );
	
	struct component *comp;

	int m = n/G; //gives the number of segments
	for(k=0;k<G-1;k++)
	{
		comp = mixmod->components[k] ;
		
		comp->in_use = TRUE;
	
		for(j=0;j<d;j++) comp->sum[j] = 0.;
		comp->sum_squared_norm = 0.;
	
		//cycle through appropriate data
		for(i=k*m;i<(k+1)*m;i++)
		{
			o = order[i];
			mixmod->z[ o ] = k ;
			if( mixmod->d   == 1  )
				component_add_to_component_uni( comp, mixmod->y_uni[o], 1);
			else
				component_add_to_component( comp, mixmod->Y[o], 1 );
		}
	}
	
	//special case for last group
	comp = mixmod->components[G-1];
	
	comp->in_use = TRUE;

	for(j=0;j<d;j++) comp->sum[j] = 0.;
	comp->sum_squared_norm = 0.;
	
	
	for(i=(G-1)*m;i<n;i++){
		o = order[i];
		mixmod->z[o] = G-1;
		if( mixmod->d == 1 )
			component_add_to_component_uni( comp, mixmod->y_uni[o], 1 );
		else
			component_add_to_component( comp, mixmod->Y[o], 1 );
	}
	
	//compute the log_prob for each of the components
	for(k=0;k<G;k++){
		GMM_recompute_marginal_likelihood_component( k, mixmod );
		//recompute_marginal_likelihood_component(mixmod->components[k],mixmod->alpha,mixmod->delta,mixmod->kappa,mixmod->gamma,mixmod->prior_mu,mixmod->d);
	}	

	for(k=0;k<G;k++){
		mixmod->whereis[k] = k;
	}

	free( order );

	return(TRUE);
}


double get_max(double *x,int len)
/*returns maximum of a vector x*/
{
	int i;
	double max=x[0];
	
	if(len > 1){
		for(i=1;i<len;i++){
			if(x[i]>max)
				max = x[i];
		}
	
	}
	return(max);
}

double get_min(double *x,int len)
/*returns maximum of a vector x*/
{
	int i;
	double min=x[0];
	
	if(len > 1){
		for(i=1;i<len;i++){
			if(x[i]<min)
				min = x[i];
		}
	
	}
	return(min);
}

int get_imax(int *x,int len)
/*returns maximum of a vector x*/
{
	int i;
	int max=x[0];
	
	if(len > 1){
		for(i=1;i<len;i++){
			if(x[i]>max)
				max = x[i];
		}
	
	}
	return(max);
}


int sample_discrete( double *weights, int len )
{
	/*sample once from a multinomial distribution with weights*/
	
	int i=0;
	double w , u;
	
	u = runif(0.0,1.0) ;
	
	w = weights[0];
	
	while( w < u && i < len )
	{
		i++ ;
		w += weights[i];		
	}
	
	return(i);	
}


int get_ind_max(double *x,int len)
/*returns index of the maximum of vector x
  cases with ties: lowest index from 0 returned*/
{
	double max;
	int i=0;
	
	max = get_max(x,len);
	
	while(1){	
		if(x[i]==max) break;
		i++;
	}
	
	return i;
}

void set_prior_hyperparameters_x(  struct mix_mod *mixmod , double gamma_0, double kappa_0)
/*This function initializes the prior hyperparmaters *rougly* analagously to Richardson and Green
... althoguh their ideas do not carrry over directley to d-dimensions, an attempt is made to counteract this*/
{
	
	double c = 4., d = 4.;
	
	mixmod->shape_gamma = 2. * c * c ;
	mixmod->rate_gamma = mixmod->shape_gamma / gamma_0 ;
	
	mixmod->shape_kappa = 2. * d * d ;
	mixmod->rate_kappa = mixmod->shape_kappa / kappa_0 ;

}


void get_hyperparameters_x( double *h, double gamma_0, double kappa_0 )
{
	double c = 4., d = 4.;
	
	h[0] = 2. * c * c ;
	h[1] = h[0] / gamma_0 ;
	
	h[2] = 2.* d * d ;
	h[3] = h[2] / kappa_0 ;
	
}



void set_prior_hyperparameters(struct mix_mod *mixmod,int type)
/*This function initializes the prior hyperparmaters *rougly* analagously to Richardson and Green
... althoguh their ideas do not carrry over directley to d-dimensions, an attempt is made to counteract this*/
{
	
	int i,j;
	double max,min,sumR2=0.,*x,*range,*low_range,h;
	
	switch(type)
	{
	
	
		case RICHARDSON_AND_GREEN:
	
			/*priors are assigned the same way as R&G allowing for the fact that some of the
			parameterizations here are slightly different*/
	
			/*need to determine range in each direction*/
			x = (double *)calloc(mixmod->n,sizeof(double));
			range = (double *)calloc(mixmod->d,sizeof(double));
			low_range = (double *)calloc(mixmod->d,sizeof(double));
	
			/*determine range and low_range for each dimension*/
			for(j=0;j<mixmod->d;j++){
	
				for(i=0;i<mixmod->n;i++)
					x[i] = mixmod->Y[i][j];
		
				max = get_max(x,mixmod->n);
				min = get_min(x,mixmod->n);
			
				range[j] = max-min;
				low_range[j] = min;
		
				sumR2 += pow(range[j],2.);
		
			}
	
			mixmod->gamma = 0.02*(sumR2/mixmod->d);
			
			/*updating of gamma would be default in R&G's model, so here I overwrite mixmod->update_gamma*/
			
			mixmod->update_gamma = TRUE;
	
			if(mixmod->update_gamma){
			
				mixmod->shape_gamma = 0.4;
	
				h = 10.*mixmod->d/sumR2;
	
				mixmod->rate_gamma = 4.*h;
	
			}
			
			mixmod->update_gamma = FALSE;
	
			/*for(j=0;j<mixmod->d;j++){
				mixmod->prior_mu[j] = low_range[j] + 0.5*range[j];
			}*/
	
			mixmod->delta = 4.;
		
			mixmod->kappa = mixmod->d/sumR2;

			///TRIONA HAS INCLUDED THIS TO TEST VALGRIND
			mixmod->update_kappa = TRUE;

			
			/*include option here for update of kappa... slightly different to
				R&G's implementation*/
			if(mixmod->update_kappa){
			
				mixmod->shape_kappa = 2.;
				
				mixmod->rate_kappa = 0.02;	
				
			}
			
			mixmod->update_kappa = FALSE;
	
			/*finally parameter for symmetric Dirichlet on weights*/
			mixmod->alpha = 1.;
		

			free(x);
			free(range);
			free(low_range);
	
		break;
		
		case NOBILE_AND_FEARNSIDE:
		
		break;
		
		case CUSTOM:
		
		break;
		
	}

}

void set_prior_on_number_of_components(struct mix_mod *mixmod,int type)
/*set prior to either that by Nobile and Fearnside or Richardson and Green*/
{

	int i;
	
	switch(type)
	{
		case RICHARDSON_AND_GREEN:
		
			//take all models equally likely apriori
			for(i=1;i<mixmod->maxgroups+1;i++){
				mixmod->log_prior_G[i] = -log((double)mixmod->maxgroups);
			}
			//do not include hprameter updates on apoisson rate
			mixmod->update_lambda = FALSE;
			
		break;
		
		case NOBILE_AND_FEARNSIDE:
		
			/*take a poisson one prior on number of components*/
			for(i=1;i<mixmod->maxgroups+1;i++){
				mixmod->log_prior_G[i] = -lgamma((double)i+1.);
			}
			/*do not include hparameter updates on a poisson rate*/
			mixmod->update_lambda = FALSE;
		
		break;
		
		case CUSTOM:

			/*keep this option in here with an aim to putting in later...*/

		break;
		
		case HYPERPRIOR_LAMBDA:
		
			/*in this case include a hyperprior on lambda which gives Poisson 
				rate of no. of compoonentns*/
				
			mixmod->lambda = 1.;
			for(i=1;i<mixmod->maxgroups+1;i++){
				mixmod->log_prior_G[i] = i*log(mixmod->lambda) - mixmod->lambda - lgamma(i+1.);
			}
			
			mixmod->update_lambda = TRUE;
		
		break;
	
	}


}




int update_allocations_with_gibbs( struct mix_mod *mixmod )
{

	int i, ind, k, g, d=mixmod->d, G=mixmod->G, *order, *wis;
	double *probs, x_, *x, *lpp_store , max, z, lc, lp, lcm, lpp ;
	
	probs = (double *)calloc( G, sizeof(double)  );
	lpp_store = (double *)calloc( G, sizeof(double) );
	order = (int *)calloc( mixmod->n, sizeof(int) );
	
	struct component *comp_c, *comp_p ;
	
	//rewrite this using the new function for returning the marginal likelihood of a component
	
	//randomize the order
	for(i=0;i<mixmod->n;i++) order[i] = i;
		
	random_ranshuffle( order, mixmod->n );
	
	//for( i=0; i< (int)(runif(0.,1.)* mixmod->n); i++ ){
	for( i=0; i<mixmod->n; i++ ){
	
		ind = order[i];
		
		//point to the entry
		if( d > 1 )	x = mixmod->Y[ind]; else x_ = mixmod->y_uni[ind] ;
	
	
		//current group
		g = mixmod->z[ind];
		
		//component
		comp_c = mixmod->components[ mixmod->whereis[g] ];
		
		lc = comp_c->log_prob;
		
		if( d > 1 ) 
			component_add_to_component( comp_c, x, -1 ) ;
		else
			component_add_to_component_uni( comp_c, x_, -1);
		
		lcm = GMM_return_marginal_likelihood_component( comp_c, mixmod ) ;
		
		max = 0.;
		
		wis =  mixmod->whereis ; 
		
		//cycle through remaining groups
		for( k=0; k<G; k++ )
		{
			
			if( k != g ){
			
				comp_p = mixmod->components[ wis[k] ] ;
				
				lp = comp_p->log_prob ; 
				
				if( d > 1 )
					component_add_to_component( comp_p, x, 1 ) ;
				else
					component_add_to_component_uni( comp_p, x_, 1 );
				
				lpp = GMM_return_marginal_likelihood_component( comp_p, mixmod ) ;
				
				lpp_store[k] = lpp;
			
				//get log relative probability for Gibbs
				probs[k] = lpp + lcm - ( lc + lp ) ;
				
				if( d > 1 )
					component_add_to_component( comp_p, x, -1 ) ;
				else
					component_add_to_component_uni( comp_p, x_, -1 );
					
			}else{
				probs[k] = 0.0;
			}
			
			if( probs[k] > max ) max = probs[k] ;
		
		}
		
		//max = get_max( probs, G );
		
		z = 0.;
		
		for( k=0; k<G; k++ ){
			probs[k] -= max;
			probs[k] = exp( probs[k] );
			z += probs[k];
		}
			
		for( k=0; k<G; k++ ) probs[k] /= z;

		k = sample_discrete( probs, G );

		if( k != g )
		{
			
			comp_p = mixmod->components[ wis[k] ] ;
			
			mixmod->z[ind] = k;
			
			if( d > 1 )
				component_add_to_component( comp_p, x, 1 );
			else
				component_add_to_component_uni( comp_p, x_, 1 );
			
			comp_p->log_prob = lpp_store[k] ;
			
			comp_c->log_prob = lcm ;
			
		}
		else
		{
			if( d > 1 )
				component_add_to_component( comp_c, x, 1 );
			else
				component_add_to_component_uni( comp_c, x_, 1 );
		}
	
	}
	
	
	free( probs );
	free( order );
	free( lpp_store );
	
	return(TRUE);
}


void update_allocations_with_metropolis_move_1(struct mix_mod *mixmod,int *accepted,int *proposed)
	/*this performs the move M1 taken from Nobile and Fearnside (2007) Stats and Computing 17: p147-162*/
{

	//printf("\nwithin Move 1...\n");
	
	/*DO NOT PERFORM THIS MOVE UNLESS THERE ARE AT LEAST TWO COMPONENTS*/
	if(mixmod->G < 2){
		
		return;
	
	}

	*proposed += 1;

	int i, k , d = mixmod->d, g1, g2 , ntot,*indexes,*proposed_alloc, *z;
	double p,log_acceptance,ag1=1.,ag2=1., y_, *y;
	struct component *component_g1,*component_g2, *curr_component_g1, *curr_component_g2;
	
	/*the integers g1 and g2 give the components, and ig1 and ig2 their whereis value*/
	/*the doubles ag1 and ag2 give the parameters to the Beta distribution for generating p_1
		take default uniform*/
	
	/*sample the two components*/
	g1 = (int)( runif(0.0,1.0) * mixmod->G ) ;  
	g2 = g1;
	while(g2 == g1){
		g2 = (int)( runif(0.0,1.0) * mixmod->G ) ;
	}
	
	/*find where in mixmod->components g1 and g2 are*/
	
	curr_component_g1 = mixmod->components[ mixmod->whereis[g1] ] ;
	curr_component_g2 = mixmod->components[ mixmod->whereis[g2] ] ;  

	
	ntot = curr_component_g1->n_g + curr_component_g2->n_g;
	
	/*DO NOT PERFORM THIS MOVE UNLESS ntot > 0*/
	if(ntot == 0){
	
		return;
		
	}
	
	/*allocate space for proposal stuff...*/
	
	component_g1 = component_create( mixmod->d );
	component_g2 = component_create( mixmod->d );
	
	/*allocate a vector to keep track of indexes*/
	indexes = (int *)calloc(ntot,sizeof(int));
	proposed_alloc = (int *)calloc(ntot,sizeof(int));
	
	k=0;
	for(i=0;i<mixmod->n;i++){
		if(mixmod->z[i] == g1 || mixmod->z[i] == g2){
			indexes[k] = i;
			k+=1;
		}
	}

	/*generate p and begin reallocation*/
	p = rbeta( ag1, ag2 ) ; 
	
	/*now propose reallocation*/
	
	for(i=0;i<ntot;i++){
	
		if( d > 1 )
			y = mixmod->Y[ indexes[i] ];
		else
			y_ = mixmod->y_uni[ indexes[i] ];
	
		if( runif(0.0,1.0) < p){
		
			/*reallocate to g1*/
			proposed_alloc[i] = g1;
			if( d > 1 )
				component_add_to_component( component_g1, y, 1 ) ;
			else
				component_add_to_component_uni( component_g1, y_, 1 );
			
		}else{
		
			/*reallocate to g2*/
			proposed_alloc[i] = g2;
			if( d > 1 )
				component_add_to_component( component_g2, y, 1 ) ;
			else
				component_add_to_component_uni( component_g2, y_, 1 );
		
		}		
		
	}
	
	/*evaluate log of acceptance probability*/
	GMM_recompute_marginal_likelihood_component_0( component_g1, mixmod );
	GMM_recompute_marginal_likelihood_component_0( component_g2, mixmod );
		
   log_acceptance = component_g1->log_prob + component_g2->log_prob - curr_component_g1->log_prob - curr_component_g2->log_prob
                       + lgamma(mixmod->alpha + curr_component_g1->n_g) + lgamma(mixmod->alpha + curr_component_g2->n_g) 
                       - lgamma(mixmod->alpha + component_g1->n_g) - lgamma(mixmod->alpha + component_g2->n_g);
	//printf("\nThe value of log acceptance is %.10f,",log_acceptance);
	
	if( log( runif(0.0,1.0) ) < log_acceptance ){
	
		/*do the swap!!!*/
		
		*accepted += 1;
		
		z = mixmod->z;
		
		for(i=0;i<ntot;i++){
			z[ indexes[i] ] = proposed_alloc[i];
		}
		
		/*copy over the accepted components*/
		component_g1->in_use = component_g2->in_use = TRUE;
		
		copy_component( component_g1, curr_component_g1 );
		copy_component( component_g2, curr_component_g2 );
	
	}
	

	component_destroy( component_g1 );
	component_destroy( component_g2 );
	
	free(indexes);
	free(proposed_alloc);
	
	return;
}

void update_allocations_with_metropolis_move_2(struct mix_mod *mixmod,int *accepted,int *proposed)
/*this performs the move M2 taken from Nobile and Fearnside (2007) Stats and Computing 17: p147-162*/
{

	/*DO NOT PERFORM THIS MOVE UNLESS THERE ARE AT LEAST TWO COMPONENTS*/
	if(mixmod->G < 2){
		
		return;
	
	}

	int i,ii,k,d=mixmod->d,g1,g2,curr_n_g1,m,c=0;
	int *indexes,*order;
	double log_acceptance;
	struct component *component_g1,*component_g2, *curr_component_g1, *curr_component_g2;
	
	/*sample the two components*/
	g1 = (int)( runif(0.0,1.0) * mixmod->G ) ; 
	g2 = g1;
	while(g2 == g1){
		g2 = (int)( runif(0.0,1.0) * mixmod->G ) ;
	}	

	/*find where in mixmod->components g1 and g2 are*/
	curr_component_g1 = mixmod->components[ mixmod->whereis[g1] ] ;
	curr_component_g2 = mixmod->components[ mixmod->whereis[g2] ] ;

	if( curr_component_g1->n_g == 0 ){
		/*cannot perform move in this case... return*/
		return;
	}	
	
	/*current component sizes*/
	curr_n_g1 = curr_component_g1->n_g;
	
	/*allocate the candidate components*/
	
	component_g1 = component_create( mixmod->d );//(struct component *)malloc(sizeof(struct component));
	component_g2 = component_create( mixmod->d );//(struct component *)malloc(sizeof(struct component));
	
	*proposed += 1;
	
	order = (int *)calloc(curr_n_g1,sizeof(int));
	for(i=0;i<curr_n_g1;i++){
		order[i] = i;
	}
	
	/*shuffle the order*/
	random_ranshuffle( order, curr_n_g1 );
	
	indexes = (int *)calloc(curr_n_g1,sizeof(int));
	for(i=0;i<mixmod->n;i++){
		if(mixmod->z[i] == g1){
			indexes[c] = i;
			c += 1;
		}
	}
	
	m = (int)( runif(0.0,1.0) * curr_n_g1 ) ;     
	
	/*try to reallocate the first i=1,...,m in indexes[order[i]]*/
	
	copy_component( curr_component_g1, component_g1 );
	copy_component( curr_component_g2, component_g2 );
	
	for(i=0;i<m;i++){
	
		k = order[i];
		ii = indexes[k];
		
		if( d > 1 )
		{
			component_add_to_component( component_g1, mixmod->Y[ii], -1 );
			component_add_to_component( component_g2, mixmod->Y[ii],  1 );
		}
		else
		{
			component_add_to_component_uni( component_g1, mixmod->y_uni[ii], -1 );
			component_add_to_component_uni( component_g2, mixmod->y_uni[ii],  1 );
		}
	}
	
	/*compute the log acceptance probability*/
	
	GMM_recompute_marginal_likelihood_component_0( component_g1, mixmod );
	GMM_recompute_marginal_likelihood_component_0( component_g2, mixmod );
	
	log_acceptance = component_g1->log_prob + component_g2->log_prob - curr_component_g1->log_prob - curr_component_g2->log_prob
							+ log(component_g1->n_g + m) - log(component_g2->n_g) + lgamma(component_g1->n_g+m+1.) + lgamma(component_g2->n_g-m+1.)
								- lgamma(component_g1->n_g + 1.) - lgamma(component_g2->n_g+1.);
	
	if( log( runif(0.0,1.0) ) < log_acceptance ){
		
		/*do the swap*/
		
		*accepted += 1;
		
		for(i=0;i<m;i++){
			k = order[i];
			ii = indexes[k];
			mixmod->z[ii] = g2;		
		}
		
		/*copy over the accepted components*/
		copy_component( component_g1, curr_component_g1 );
		copy_component( component_g2, curr_component_g2 );
	
	}
	

	component_destroy(component_g1);
	component_destroy(component_g2);

	free(order);
	free(indexes);
	
	return;
}



void update_allocations_with_metropolis_move_3(struct mix_mod *mixmod ,int *accepted,int *proposed)
/*this performs the move M3 taken from Nobile and Fearnside (2007) Stats and Computing 17: p147-162*/
{

	/*DO NOT PERFORM THIS MOVE UNLESS THERE ARE AT LEAST TWO COMPONENTS*/
	if(mixmod->G < 2){
		
		return;
	
	}

	int i,ii,k,g1,g2,ntot,c=0,d = mixmod->d;
	int *indexes,*order,*proposed_allocation;
	double w,log_acceptance, log_transition_z_to_zprime=0.,log_transition_zprime_to_z=0.,l1,l2,p1, lm1, lm2;
	struct component *component_g1, *component_g2, *curr_component_g1, *curr_component_g2;
	
	*proposed += 1;
	
	/*sample the two components*/
	g1 = (int)( runif(0.0,1.0) * mixmod->G ) ; 
	g2 = g1;
	while(g2 == g1){
		g2 = (int)( runif(0.0,1.0) * mixmod->G ) ;
	}	

	/*find where in mixmod->components g1 and g2 are*/
	curr_component_g1 = mixmod->components[ mixmod->whereis[g1] ];
	curr_component_g2 = mixmod->components[ mixmod->whereis[g2] ];
	
	ntot = curr_component_g1->n_g + curr_component_g2->n_g;
	
	if(ntot == 0){
		/*cannot perform if both components empty*/
		return;
	}
	
	indexes = (int *)calloc(ntot,sizeof(int));
	order = (int *)calloc(ntot,sizeof(int));
	proposed_allocation = (int *)calloc(ntot,sizeof(int));
	
	/*this move can still be done if either component empty*/
	
	for(i=0;i<mixmod->n;i++){	
		if(mixmod->z[i] == g1 || mixmod->z[i] == g2){
			indexes[c] = i;
			c += 1;
		}
	}
	
	for(i=0;i<ntot;i++){
		order[i] = i;
	}
	
	/*shuffle the order*/
	random_ranshuffle( order, ntot ) ;
	
	component_g1 = component_create( mixmod->d );
	component_g2 = component_create( mixmod->d );
	
	/*randomizing this part should give a 50% acceptance*/

	k = order[0];
	ii = indexes[k];
	
	if( mixmod->z[ii] == g1 ){
	
		if( d > 1 )
			component_add_to_component( component_g1, mixmod->Y[ii], 1 ) ;
		else
			component_add_to_component_uni( component_g1, mixmod->y_uni[ii], 1 ) ;

		proposed_allocation[0] = g1;
		
		log_transition_z_to_zprime = log(0.5) - log(1.);
	
	}else{
	
		/*put it into component_g2*/
		
		if( d > 1 )
			component_add_to_component( component_g2, mixmod->Y[ii], 1 );
		else
			component_add_to_component_uni( component_g2, mixmod->y_uni[ii], 1 );
			
		proposed_allocation[0] = g2;

		log_transition_z_to_zprime = log(0.5) - log(1.);

	}
	
	/*identify this as being in the same component as original
		...this allows the computing of the backwards probability*/
	
	log_transition_zprime_to_z += log(0.5) - log(1.);
	
	GMM_recompute_marginal_likelihood_component_0( component_g1, mixmod );
	GMM_recompute_marginal_likelihood_component_0( component_g2, mixmod );
	
	for(i=1;i<ntot;i++){
	
		k = order[i];
		ii = indexes[k];
		
		/*do the proposed component values, then sample and make changes*/
		
		/*compute probability generated from g1*/
		if( d > 1 )
			lm1 = GMM_compute_marginal_likelihood_with_inclusion_in_component( mixmod->Y[ii], component_g1, mixmod );
		else
			lm1 = GMM_compute_marginal_likelihood_with_inclusion_in_component_uni( mixmod->y_uni[ii], component_g1, mixmod ) ;
		
		l1 = lm1 + component_g2->log_prob ; 
		
		//l1 = compute_log_data_probability_with_inclusion_in_component(mixmod->Y[ii],component_g1,delta,kappa,gamma,xi,d)
			//	+ compute_log_data_probability_component(component_g2,delta,kappa,gamma,xi,d);
		
		if( d > 1 )
			lm2 = GMM_compute_marginal_likelihood_with_inclusion_in_component( mixmod->Y[ii], component_g2, mixmod ) ;
		else
			lm2 = GMM_compute_marginal_likelihood_with_inclusion_in_component_uni( mixmod->y_uni[ii], component_g2, mixmod ) ;
		
		l2 = component_g1->log_prob + lm2 ;
		
		//l2 = compute_log_data_probability_component(component_g1,delta,kappa,gamma,xi,d)
			//	+ compute_log_data_probability_with_inclusion_in_component(mixmod->Y[ii],component_g2,delta,kappa,gamma,xi,d);
		
		w = exp( l1 - l2 );
		
		p1 = w/(1.+w);
		
		/*make a draw*/
		
		if( runif(0.0,1.0) < p1 ){
		
			/*put it in g1*/
			
			if( d > 1 )
				component_add_to_component( component_g1, mixmod->Y[ii], 1 );
			else
				component_add_to_component_uni( component_g1, mixmod->y_uni[ii], 1 );

			if(!(mixmod->z[ii] == g1)){
				log_transition_z_to_zprime += log(p1);
				log_transition_zprime_to_z += log(1.-p1);
			}
			
			proposed_allocation[i] = g1;
			
			component_g1->log_prob = lm1;
			
		}else{

			/*put it in g2*/
			
			if( d > 1 )
				component_add_to_component( component_g2, mixmod->Y[ii], 1 );
			else
				component_add_to_component_uni( component_g2, mixmod->y_uni[ii], 1 );
			
			if(!(mixmod->z[ii] == g2)){
				log_transition_z_to_zprime += log(1.-p1);
				log_transition_zprime_to_z += log(p1);
			}
		
			proposed_allocation[i] = g2;
			
			component_g2->log_prob = lm2;
		}
		
	}
	

	/*compute the acceptance probability*/
	
	//GMM_recompute_marginal_likelihood_component_0( component_g1, mixmod );
	//GMM_recompute_marginal_likelihood_component_0( component_g2, mixmod );
	
	//recompute_marginal_likelihood_component(component_g1,alpha,delta,kappa,gamma,xi,d);
	//recompute_marginal_likelihood_component(component_g2,alpha,delta,kappa,gamma,xi,d);
	
	log_acceptance = component_g1->log_prob + component_g2->log_prob 
						- curr_component_g1->log_prob - curr_component_g2->log_prob
						+ log_transition_zprime_to_z - log_transition_z_to_zprime;
	
	//printf("\n%.4f \t %.4f \t %.4f ",log_acceptance,log_transition_zprime_to_z,log_transition_z_to_zprime );
		//printf("\n move M3 lg acc %.10f",log_acceptance);				
	if( log( runif(0.0,1.0) ) < log_acceptance ){
	
		*accepted += 1;
	
		/*accept the move and update all quantities*/
		
		/*allocations first*/
		
		component_g1->in_use = TRUE;
		component_g2->in_use = TRUE;

		copy_component( component_g1, curr_component_g1 );
		copy_component( component_g2, curr_component_g2 );
		
		
		for(i=0;i<ntot;i++){
				
			k = order[i];
			ii = indexes[k];
			mixmod->z[ii] = proposed_allocation[i];
			
		}
	
	}
	
	/*free up all memory*/
	
	component_destroy(component_g1);
	component_destroy(component_g2);

	free(indexes);
	free(order);
	free(proposed_allocation);
	
	return;

}


/*eject and absorb moves*/

void update_allocations_with_ejection_move(struct mix_mod *mixmod ,int *accepted,int *proposed,double pr_ej_G,double pr_ej_Gp1)
/*this is the ejection move for one comonent ejecting another*/
{

	int i,ii,g1,g2,ig1,ig2,ntot,c=0,d = mixmod->d;
	int *indexes,*proposed_allocation,G = mixmod->G;
	double w,a,prob_put_in_g2,log_acceptance, log_transition_z_to_zprime=0.,log_transition_zprime_to_z=0.;
	struct component *component_g1, *component_g2, *curr_component_g1;
	
	*proposed += 1;
	
	/*sample the ejecting component*/
	g1 = (int)( runif(0.0,1.0) * mixmod->G ) ;
	g2 = G;

	/*find where in mixmod->components g1 is*/
	curr_component_g1 = mixmod->components[ mixmod->whereis[g1] ] ;
	
	/*if this component is empty we need a special case*/
	
	component_g1 = component_create( mixmod->d );
	component_g2 = component_create( mixmod->d );
	
	ntot = curr_component_g1->n_g;	
	
	if(ntot > 0){
		//printf("\nin there 1");
	
		/*this is the case for ejecting from a non-empty component*/
		
		
		indexes = (int *)calloc(ntot,sizeof(int));
		proposed_allocation = (int *)calloc(ntot,sizeof(int));
		
		c = 0;
		
		for(i=0;i<mixmod->n;i++){
			if(mixmod->z[i] == g1){
				indexes[c] = i;
				c += 1;
			}
		}	
		
		/*copy the contents of mixmod->components[ig1] into*/
		copy_component( curr_component_g1, component_g1 );	
		
		/*generate the probability of assignment to the new component*/
		
		if(ntot < 4){
			/*just set a = 100*/
			a = 100.;
			prob_put_in_g2 = rbeta(a,a) ; 
		}else{
			a = a_table[ntot-1];
			prob_put_in_g2 = rbeta(a,a) ;
		}
	
		/*now reassign or not*/
		
		for(i=0;i<ntot;i++){
		
			ii = indexes[i];
			
			if(runif(0.0,1.0)<prob_put_in_g2){
			
				/*then move this point to g2*/
				
				/*first take out of component_g1*/
				
				if( d > 1 )
				{
					component_add_to_component( component_g1, mixmod->Y[ii], -1 );
					component_add_to_component( component_g2, mixmod->Y[ii], 1 );
				}
				else
				{
					component_add_to_component_uni( component_g1, mixmod->y_uni[ii], -1 );
					component_add_to_component_uni( component_g2, mixmod->y_uni[ii], 1 );
				}
				
				proposed_allocation[i] = g2;
			
			}else{
			
				proposed_allocation[i] = g1;
				
			}
		
		
		}
	
	}
	
	GMM_recompute_marginal_likelihood_component_0( component_g1, mixmod );
	GMM_recompute_marginal_likelihood_component_0( component_g2, mixmod );
	
	/*compute the acceptance probability, remembering to add all necessary normalizing constants*/
	
	w = log_normalizing_constant_model(G+1,mixmod) - log_normalizing_constant_model(G,mixmod);
	
	//printf("\nBirth: w = %.10f ",w);

// 	log_transition_z_to_zprime = log(pr_ej_G);
// 		printf("log_tr... is %lf",log_transition_z_to_zprime);
	if(ntot > 0){
		log_transition_z_to_zprime += log(pr_ej_G) + lgamma(2.*a) - 2.*lgamma(a) + lgamma(a + component_g1->n_g) 
									+ lgamma(a + component_g2->n_g) - lgamma(2.*a + ntot);
	}
	
	log_transition_zprime_to_z = log(1.-pr_ej_Gp1);
	
	log_acceptance = w + component_g1->log_prob + component_g2->log_prob - curr_component_g1->log_prob 
							- log_transition_z_to_zprime + log_transition_zprime_to_z + mixmod->log_prior_G[G+1] - mixmod->log_prior_G[G];
							
	
	//printf("\nlog_acceptance = %.10f",log_acceptance);
	
	if(log(runif(0.0,1.0)) < log_acceptance){
	
		*accepted += 1;
	
		/*update the model structure*/
		
		mixmod->G += 1;
		
		/*relabel the appropriate indexes*/
		if(ntot>0){
			for(i=0;i<ntot;i++){
				ii = indexes[i];
				mixmod->z[ii] = proposed_allocation[i];
			}
		}
		
		int new_whereis=0;
		
		/*begin to encode an unused component by -1 in mixmod->whereis*/
		while(mixmod->components[new_whereis]->in_use == TRUE){
			new_whereis += 1;
		}
		
		
		/*put new componenet G in new_whereis*/
		/*this will now be indexed as component G*/
		
		mixmod->whereis[G] = new_whereis;
		//mixmod->components[new_whereis]->in_use = TRUE;
		
		component_g1->in_use = TRUE;
		component_g2->in_use = TRUE;
		
		copy_component( component_g1, curr_component_g1 );
		copy_component( component_g2, mixmod->components[new_whereis] );
		
		/*now do a swap between component G and one of the others...*/
		
		/*generate component randomly*/
		
		g1 = (int)( runif(0.0,1.0) * (G + 1) ) ;
		
		if( g1 != G ){
		
			/*do a swap!*/
			
			ig1 = mixmod->whereis[g1];
			ig2 = mixmod->whereis[G];
			
			mixmod->whereis[g1] = ig2;
			mixmod->whereis[G] = ig1;
			
			
			/*copy_component(mixmod->components[ig1],component_g1,d);
			copy_component(mixmod->components[ig2],mixmod->components[ig1],d);
			copy_component(component_g1,mixmod->components[ig2],d);*/

			/*relabel the appropriate indexes*/
			for(i=0;i<mixmod->n;i++){
				//ii = indexes[i];
				if(mixmod->z[i] == g1){
					mixmod->z[i] = G;
				}else if(mixmod->z[i] == G){
					mixmod->z[i] = g1;
				}
			}
			
		}
		
		/*for(i=0;i<mixmod->G;i++){
			print_component(i,mixmod);
		}*/
	
	}
	
	
	/*for(k=0;k<mixmod->maxgroups;k++)
		printf("\nComponent %d is in %d",k,mixmod->whereis[k]);*/
	
	component_destroy(component_g1);
	component_destroy(component_g2);
	
	//free(component_g1);
	//free(component_g2);

	if(ntot>0){
		free(indexes);
		free(proposed_allocation);
	}
	
	return;
}


void update_allocations_with_absorb_move(struct mix_mod *mixmod,int *accepted,int *proposed,double pr_ej_G,double pr_ej_Gm1)
{
	int i,ii,j,k,g1,g2,ntot,c=0,d = mixmod->d;
	int *indexes,*proposed_allocation,G = mixmod->G,n_g2;
	double w,a,log_acceptance,log_transition_z_to_zprime=0.,log_transition_zprime_to_z=0.;
	struct component *component_g1, *curr_component_g1, *curr_component_g2;
	
	*proposed += 1;
	
	/*choose component to absorb into and to absorb*/
	g1 = (int)( runif(0.0,1.0) * mixmod->G ) ;
	g2 = g1;
	while(g2 == g1){
		g2 = (int)( runif(0.0,1.0) * mixmod->G ) ;
	}
	
	/*find where in mixmod->components g1 and g2 are*/
	curr_component_g1 = mixmod->components[ mixmod->whereis[g1] ] ;
	curr_component_g2 = mixmod->components[ mixmod->whereis[g2] ] ;

	/*use a component to store the proposed*/
	
	component_g1 = component_create( mixmod->d );
	
	/*form the proposed component by combining the two other components...*/
	
	ntot = curr_component_g1->n_g + curr_component_g2->n_g;

	n_g2 = curr_component_g2->n_g;

	copy_component( curr_component_g1, component_g1 );
	
	if(n_g2 > 0){
	
		/*this is the case for non-empty components*/
		indexes = (int *)calloc(n_g2,sizeof(int));
		proposed_allocation = (int *)calloc(n_g2,sizeof(int));
		
		c = 0;
		
		for(i=0;i<mixmod->n;i++){
			if(mixmod->z[i] == g2){
				indexes[c] = i;
				c += 1;
			}
		}
		
		/*now just place everything in component_g1*/
		
		for(i=0;i<n_g2;i++){
		
			ii = indexes[i];
			
			if( d > 1 )
				component_add_to_component( component_g1, mixmod->Y[ii], 1 );	
			else
				component_add_to_component_uni( component_g1, mixmod->y_uni[ii], 1 );		
		}	
	
	}
	
	GMM_recompute_marginal_likelihood_component_0( component_g1, mixmod );
	//recompute_marginal_likelihood_component(component_g1,alpha,delta,kappa,gamma,xi,d);	
	
	/*compute the acceptance probability, remembering to add all necessary normalizing constants*/
	
	/* w = log of difference in normalizing constants*/
	
	w = log_normalizing_constant_model(G-1,mixmod) - log_normalizing_constant_model(G,mixmod);
	
	//printf("\nDeath: w = %.10f ",w);
	
	log_transition_zprime_to_z = log(pr_ej_Gm1);
		
	if(ntot > 0){
	
		if(ntot < 4){
			a = 100.;
		}else{
			a = a_table[ntot-1];
		}
	
		log_transition_zprime_to_z +=  lgamma(2.*a) - 2.*lgamma(a) + lgamma(a + curr_component_g1->n_g) + lgamma(a + curr_component_g2->n_g) - lgamma(2.*a + ntot);
		
	}
	
	log_transition_z_to_zprime = log(1.-pr_ej_G);	
	
	log_acceptance = w + component_g1->log_prob - curr_component_g1->log_prob - curr_component_g2->log_prob
							- log_transition_z_to_zprime + log_transition_zprime_to_z + dpois( (double)(G-1), 1., 1 ) - dpois( (double)G, 1., 1 ) ;
	

	//printf("\nlog_acceptance absorb move is %.10f\n",log_acceptance);

	if(log(runif(0.0,1.0)) < log_acceptance){
	
		//printf("\nlog_acceptance absorb move is %.10f\n",log_acceptance);
	
		*accepted += 1;
	
		/*update the model structure*/
		
		mixmod->G -= 1;
		
		/*relabel the appropriate indexes*/
		if(n_g2 > 0){
			for(i=0;i< n_g2;i++){
				ii = indexes[i];
				mixmod->z[ii] = g1;
			}
		}
		
		copy_component( component_g1, curr_component_g1 ) ;
		
		component_refresh( curr_component_g2 );
		
		
		//mixmod->whereis[g2] = -1;
		
		/*should relabel everything from component g2 upwards*/
		
		for(k=g2+1;k<G;k++){
		
			for(i=0;i<mixmod->n;i++){
				if(mixmod->z[i] == k){
					mixmod->z[i] = k-1;
				}
			}
			
			j = mixmod->whereis[k];
			//printf("\nvalue of j = %d",j);
			mixmod->whereis[k-1] = j;
			//mixmod->whereis[k] = -1;
		
		}
		
		mixmod->whereis[G-1] = -1;
		
		/*for(k=0;k<mixmod->maxgroups;k++)
			printf("\nComponent %d is in %d",k,mixmod->whereis[k]);*/

		
	}
	
	component_destroy( component_g1 );
	
	if(n_g2 > 0){
		free(indexes);
		free(proposed_allocation);
	}	

	return ;

}

/****************experimental split/combine moves*****************/


/*int update_allocations_with_split_move(struct mix_mod *mixmod,int *accepted,int *proposed,double pr_ej_G,double pr_ej_Gp1)
{

	int i,ii,j,k,g1,g2,ig1,ig2,ntot,c=0,d=mixmod->d;
	int *indexes,*order,*proposed_allocation;
	double w,q,log_acceptance,alpha=mixmod->alpha,delta=mixmod->delta,kappa=mixmod->kappa,gamma=mixmod->gamma,*xi=mixmod->prior_mu
			,log_transition_z_to_zprime=0.,log_transition_zprime_to_z=0.,l1,l2,p1;
	struct component *component_g1,*component_g2;

	*proposed += 1;
	
	
	//ample a component to split
	
	g1 = (int)( runif(0.0,1.0) * mixmod->G ) ;
	ig1 = mixmod->whereis[g1];
	
	g2 = mixmod->G;
	
	int G = mixmod->G;
	
	ntot = mixmod->components[ig1]->n_g;
	
	//if(ntot==0){
		//return(TRUE); //abort move if component empty
	//}	
	
	indexes = calloc(ntot,sizeof(int));
	order = calloc(ntot,sizeof(int));
	proposed_allocation = calloc(ntot,sizeof(int));
		
	for(i=0;i<mixmod->n;i++){
		if(mixmod->z[i] == g1){
			indexes[c] = i;
			c += 1;
		}
	}
	
	for(i=0;i<ntot;i++){
		order[i] = i;
	}
		
	if(  ntot > 0 ) random_ranshuffle( order, ntot ) ;
	
	component_g1 =(struct component *)malloc(sizeof(struct component));
	component_g2 =(struct component *)malloc(sizeof(struct component));
	
	allocate_component(component_g1,mixmod);
	allocate_component(component_g2,mixmod);
	
	component_g1->sum_squared_norm = 0.;
	component_g2->sum_squared_norm = 0.;
	component_g1->n_g = 0;
	component_g2->n_g = 0;
	
	//k = order[0];
	//ii = indexes[k];
	
	if(ntot>0){
	for(i=0;i<ntot;i++)
	{
	
		k = order[i];
		ii = indexes[k];
		
		l1 = compute_log_data_probability_with_inclusion_in_component(mixmod->Y[ii],component_g1,delta,kappa,gamma,xi,d) 
		+ compute_log_data_probability_component(component_g2,delta,kappa,gamma,xi,d);
		
		l2 = compute_log_data_probability_component(component_g1,delta,kappa,gamma,xi,d)
				+ compute_log_data_probability_with_inclusion_in_component(mixmod->Y[ii],component_g2,delta,kappa,gamma,xi,d);
		
		w = ((alpha + component_g1->n_g)/(alpha + component_g2->n_g))*exp(l1-l2);
		
		//printf("\n w = %.10f,l1-l2 = %.10f",w,l1);
		
		p1 = w/(1.+w);
		
		if(runif(0.0,1.0)<p1){
		
			component_g1->n_g += 1;
			
			for(j=0;j<d;j++){
				component_g1->sum[j] += mixmod->Y[ii][j];
				component_g1->sum_squared_norm += mixmod->Y[ii][j] * mixmod->Y[ii][j] ;
			}
			
			log_transition_z_to_zprime += log(p1);

			proposed_allocation[i] = g1;	
		
		}else{
		
			component_g2->n_g += 1;
			
			for(j=0;j<d;j++){
				component_g2->sum[j] += mixmod->Y[ii][j];
				component_g2->sum_squared_norm += mixmod->Y[ii][j] * mixmod->Y[ii][j] ;
			}
			
			log_transition_z_to_zprime += log(1.-p1);
		
			proposed_allocation[i] = g2;
		
		}
	
	}
	
	}
	
	//compute the acceptance probability
	
	GMM_recompute_marginal_likelihood_component_0( component_g1, mixmod );
	GMM_recompute_marginal_likelihood_component_0( component_g2, mixmod );
	
	//recompute_marginal_likelihood_component(component_g1,alpha,delta,kappa,gamma,xi,d);
	//recompute_marginal_likelihood_component(component_g2,alpha,delta,kappa,gamma,xi,d);
	
	q = log_normalizing_constant_model( G+1, mixmod ) - log_normalizing_constant_model( G, mixmod );
	
	log_transition_z_to_zprime += log(pr_ej_G);
	
	log_transition_zprime_to_z = log(1-pr_ej_Gp1);
	
	log_acceptance = q + component_g1->log_prob + component_g2->log_prob - mixmod->components[ig1]->log_prob - log_transition_z_to_zprime + log_transition_zprime_to_z + mixmod->log_prior_G[G+1] - mixmod->log_prior_G[G];
	
	
	//printf("\nValue of log acc for split %.10f",log_acceptance);	
	
	if( log(runif(0.0,1.0))  <  log_acceptance  ){
	
		*accepted += 1;
		
		mixmod->G += 1;
		
		for(i=0;i<ntot;i++){
			ii = indexes[i];
			mixmod->z[ii] = proposed_allocation[i];
		}
		
		int new_whereis = 0;
		
		while(mixmod->components[new_whereis]->in_use){
			new_whereis += 1;
		}
		
		mixmod->whereis[G] = new_whereis;
		mixmod->components[new_whereis]->in_use = TRUE;
		
		copy_component(component_g1,mixmod->components[ig1]);
		copy_component(component_g2,mixmod->components[new_whereis]);
	
		//do a swap of G with one other component
		
		g1 = (int)( runif(0.0,1.0) * (G + 1) ) ;  
		
		if(g1 != G){
		
			ig1 = mixmod->whereis[g1];
			ig2 = mixmod->whereis[G];
			
			mixmod->whereis[g1] = ig2;
			mixmod->whereis[G] = ig1;
			
			for(i=0;i<mixmod->n;i++){
				if(mixmod->z[i] == g1){
					mixmod->z[i] = G;
				}else if(mixmod->z[i] == G){
					mixmod->z[i] = g1;
				}
			
			}
		
		
		}
	
	}
	
	free_component(component_g1);
	free_component(component_g2);
	
	free(component_g1);
	free(component_g2);
	
	free(indexes);
	free(proposed_allocation);
	free(order);
	
	return(TRUE);

}*/

/*int update_allocations_with_combine_move(struct mix_mod *mixmod,int *accepted,int *proposed,double pr_ej_G,double pr_ej_Gm1)
{

	
	

	int i,ii,j,k,g1,g2,ig1,ig2,ntot,c=0,d = mixmod->d;
	int *indexes,*proposed_allocation,*order,G = mixmod->G,n_g2;
	double w,l1,l2,p1,q,a,log_acceptance,alpha=mixmod->alpha,delta=mixmod->delta,kappa=mixmod->kappa,gamma=mixmod->gamma,*xi=mixmod->prior_mu
			,log_transition_z_to_zprime=0.,log_transition_zprime_to_z=0.;
	struct component *component_g1,*component_g2,*component_new;
	
	
	
	*proposed += 1;

	//choose components to combine
	g1 = (int)( runif(0.0,1.0) * G ) ;  
	g2 = g1;
	while(g2 == g1){
		g2 = (int)( runif(0.0,1.0) * G ) ;
	}
	
	//find where in mixmod->components g1 and g2 are
	ig1 = mixmod->whereis[g1];
	ig2 = mixmod->whereis[g2];	
	
	//form the proposed component by combining the two other components...
	
	ntot = mixmod->components[ig1]->n_g + mixmod->components[ig2]->n_g;

	//if(ntot == 0){
	//	return(TRUE); //abort moves
	//}

	//use a component to store the proposed
	
	component_g1 =(struct component *)malloc(sizeof(struct component));
	component_g2 =(struct component *)malloc(sizeof(struct component));
	component_new = (struct component *)malloc(sizeof(struct component));
	
	allocate_component(component_g1,mixmod);
	allocate_component(component_g2,mixmod);
	allocate_component(component_new,mixmod);
	
	component_g1->sum_squared_norm = 0.;
	component_g1->n_g = 0;
	
	component_g2->sum_squared_norm = 0.;
	component_g2->n_g = 0;
	
	component_new->sum_squared_norm = 0.;
	component_new->n_g = 0;
	
	
	indexes = calloc(ntot,sizeof(int));
	order = calloc(ntot,sizeof(int));
	proposed_allocation = calloc(ntot,sizeof(int));

	for(i=0;i<mixmod->n;i++){	
		if(mixmod->z[i] == g1 || mixmod->z[i] == g2){
			indexes[c] = i;
			c += 1;
		}
	}
	
	for(i=0;i<ntot;i++){
		order[i] = i;
	}
	
	if( ntot > 0 ) random_ranshuffle( order, ntot ) ; 
	
	if(FALSE)
	{
	
	k = order[0];
	ii = indexes[k];
	
	if(mixmod->z[ii] == g1){
	
		component_g1->n_g += 1;
	
		for(j=0;j<d;j++){
			component_g1->sum[j] += mixmod->Y[ii][j];
			component_g1->sum_squared_norm += pow(mixmod->Y[ii][j],2.);
		}
		
	
	}else{
	
		component_g2->n_g += 1;
		
		for(j=0;j<d;j++){
			component_g2->sum[j] += mixmod->Y[ii][j];
			component_g2->sum_squared_norm += pow(mixmod->Y[ii][j],2.);
		}
	
	}
	
	component_new->n_g += 1;
		
	for(j=0;j<d;j++){
		component_new->sum[j] += mixmod->Y[ii][j];
		component_new->sum_squared_norm += pow(mixmod->Y[ii][j],2.);
	}
	
	log_transition_zprime_to_z = log(.5);
	
	log_transition_z_to_zprime = 0.;
	
	}
	
	if( ntot > 0 )
	{
	
	for(i=0;i<ntot;i++){
	
		k = order[i];
		ii = indexes[k];
	
		//compute probability generated from g1
		
		l1 = compute_log_data_probability_with_inclusion_in_component(mixmod->Y[ii],component_g1,delta,kappa,gamma,xi,d)
				+ compute_log_data_probability_component(component_g2,delta,kappa,gamma,xi,d);
		
		l2 = compute_log_data_probability_component(component_g1,delta,kappa,gamma,xi,d)
				+ compute_log_data_probability_with_inclusion_in_component(mixmod->Y[ii],component_g2,delta,kappa,gamma,xi,d);		
	
		w = ((alpha + component_g1->n_g)/(alpha + component_g2->n_g))*exp(l1-l2);
		
		p1 = w/(1.+w);
		
		if(mixmod->z[ii] == g1){
		
			component_g1->n_g += 1;
			
			for(j=0;j<d;j++){
				component_g1->sum[j] += mixmod->Y[ii][j];
				component_g1->sum_squared_norm += mixmod->Y[ii][j] * mixmod->Y[ii][j];
			}
			
			log_transition_zprime_to_z += log(p1);
		
		
		}else{
		
			component_g2->n_g += 1;
			
			for(j=0;j<d;j++){
				component_g2->sum[j] += mixmod->Y[ii][j];
				component_g2->sum_squared_norm += mixmod->Y[ii][j] * mixmod->Y[ii][j];		
			}
			
			log_transition_zprime_to_z += log(1.-p1);
		
		
		}
	
		component_new->n_g += 1;
		
		for(j=0;j<d;j++){
			component_new->sum[j] += mixmod->Y[ii][j];
			component_new->sum_squared_norm += mixmod->Y[ii][j] * mixmod->Y[ii][j] ;
		}
	
	}
	
	}
	
	//recompute_marginal_likelihood_component(component_new,alpha,delta,kappa,gamma,xi,d);
	
	GMM_recompute_marginal_likelihood_component_0( component_new, mixmod );
	
	q = log_normalizing_constant_model( G-1, mixmod ) - log_normalizing_constant_model( G, mixmod );
	
	log_transition_zprime_to_z += log(pr_ej_Gm1);
	
	log_transition_z_to_zprime = log(1.-pr_ej_G);
	
	log_acceptance = q + component_new->log_prob - mixmod->components[ig1]->log_prob - mixmod->components[ig2]->log_prob - log_transition_z_to_zprime + log_transition_zprime_to_z + mixmod->log_prior_G[G-1] - mixmod->log_prior_G[G];
							
    //printf("\nValue of log acc for combine %.10f",log_acceptance);	
	
	if(log(runif(0.0,1.0)) < log_acceptance){
	
		*accepted += 1;
		
		mixmod->G -= 1;
		
		for(i=0;i<mixmod->n;i++){
			if(mixmod->z[i] == g1 || mixmod->z[i] == g2){
				mixmod->z[i] = g1;
			}
		}
	
		copy_component(component_new,mixmod->components[ig1]);
	
		mixmod->components[ig2]->in_use = FALSE;
		
		for(k=g2+1;k<G;k++){
			
			for(i=0;i<mixmod->n;i++){
				if(mixmod->z[i] == k){
					mixmod->z[i] = k-1;
				}
			}
			
			j = mixmod->whereis[k];
			
			mixmod->whereis[k-1] = j;
		
		}
		
		mixmod->whereis[G-1] = -1;
	}
	
	
	free_component(component_g1);
	free(component_g1);
	free_component(component_g2);
	free(component_g2);
	free_component(component_new);
	free(component_new);
	
	free(indexes);
	free(order);
	free(proposed_allocation);
	

	
	return(TRUE);
}*/


/****************************************************************/

void update_hyperparameters(struct mix_mod *mixmod)
/*here try to only update the gamma hyperparameter and the precisions using Gibbs...*/
{
	
	int i, j, k, G = mixmod->G, d = mixmod->d, n_g, iter;
	double *tau_g, xi2, sq_norm, a, b, m, sd, sum, sum2, l, c2;
	
	tau_g = (double *)calloc(G,sizeof(double));
	
	/*draw tau_g*/
	
	xi2 = 0.;
	for(j=0;j<d;j++) xi2 += mixmod->prior_mu[j] * mixmod->prior_mu[j] ;
	
	/*gamma update*/
	
	if( mixmod->update_gamma )
	{
	
		for( iter=0; iter<1; iter++ )
		{

			for(i=0;i<G;i++)
			{
			
				k = mixmod->whereis[i];
		
				n_g = mixmod->components[k]->n_g;
		
				sq_norm = 0.;
				for(j=0;j<d;j++){
					l = mixmod->components[k]->sum[j] + mixmod->kappa * mixmod->prior_mu[j] ;
					sq_norm += l * l ;
				}
		
				/*precisions*/
				a = 0.5 * ( n_g * d + mixmod->delta );
				b = 0.5 * ( mixmod->components[k]->sum_squared_norm + mixmod->kappa*xi2 + mixmod->gamma - sq_norm/(n_g+mixmod->kappa) );
			
				tau_g[i] =  rgamma( a, 1.0/b ) ;
			
			}
		
			b = mixmod->rate_gamma;
		
			for(i=0;i<G;i++) b += tau_g[i];
		
			b *= 0.5;
		
			a = 0.5*( G * mixmod->delta + mixmod->shape_gamma );
		
			mixmod->gamma = rgamma( a , 1.0/b ) ;
		
		}
	
	}	
	
	
	//update both kappa and gamma

	if( mixmod->update_kappa && mixmod->update_gamma )
	{
	
		/*having generated the tau's imagine uncollapsing the mu's in order to simulate kappa and prior_mu/xi*/
	
		/*simulate the elements of mu_g conditioning on tau_g*/

	
		double **mu_g;
		mu_g = (double **)calloc(G,sizeof(double *));
		for(i=0;i<G;i++) mu_g[i] = (double *)calloc(d,sizeof(double));
		
		for( iter=0; iter<100; iter++) 
		{
		
			for(i=0;i<G;i++)
			{
	
				k = mixmod->whereis[i];
	
				for(j=0;j<d;j++)
				{
		
					m = ( mixmod->components[k]->sum[j] + mixmod->kappa*mixmod->prior_mu[j] )/( mixmod->components[k]->n_g + mixmod->kappa );
					sd = 1. / sqrt( tau_g[i] * ( mixmod->components[k]->n_g + mixmod->kappa ) ) ;
			
					mu_g[i][j] = m + rnorm( 0.0, sd );
		
				}
		
			}
	
	
	
			/*now draw a value for mixmod->kappa*/
	
		
			sum2=0.;
	
			for(i=0;i<G;i++)
			{
	
				sum = 0.;
				for(j=0;j<d;j++){
					sum += ( mu_g[i][j] - mixmod->prior_mu[j] ) * (  mu_g[i][j] - mixmod->prior_mu[j] ) ;
				}
		
				sum2 += tau_g[i] * sum;
			}
	
			a = 0.5*( G*d + mixmod->shape_kappa );
			b = 0.5*( mixmod->rate_kappa + sum2 );
	
			mixmod->kappa = rgamma( a, 1.0/b ) ; 

		
			for(i=0;i<G;i++)
			{
			
				k = mixmod->whereis[i];
		
				n_g = mixmod->components[k]->n_g;
		
				sq_norm = 0.;
				c2 = 0.;
				for(j=0;j<d;j++){
					l = mixmod->components[k]->sum[j] + mixmod->kappa * mixmod->prior_mu[j] ;
					c2 += l * mu_g[i][j] ;
					sq_norm += mu_g[i][j] * mu_g[i][j] ;
				}
		
				/*precisions*/
				a = 0.5 * ( n_g * d + mixmod->delta + 1. );
			
				b = 0.5 * ( mixmod->components[k]->sum_squared_norm  - 2. * c2 + ( n_g + mixmod->kappa ) * sq_norm + mixmod->kappa * xi2 + mixmod->gamma );
			
				tau_g[i] =  rgamma( a, 1.0/b ) ;
			
			}	
		
			b = mixmod->rate_gamma;
		
			for(i=0;i<G;i++) b += tau_g[i];
		
			b *= 0.5;
		
			a = 0.5*( G * mixmod->delta + mixmod->shape_gamma );

			//printf("\nUpdating gamma: Shape = %.10f, rate = %.10f",a,b);
		
			mixmod->gamma = rgamma( a , 1.0/b ) ; 
		
		}
	
		for(i=0;i<G;i++)
			free(mu_g[i]);
		free(mu_g);
	
	
	}
	
	if( FALSE ){
	
		/*update lambda and the prior on the number of components*/
		a = mixmod->shape_lambda+G;
		b = 1.+mixmod->rate_lambda;
		
		mixmod->lambda = rgamma( a, 1.0/b ) ;
		
		for(i=1;i<mixmod->maxgroups+1;i++){
			mixmod->log_prior_G[i] = i*log(mixmod->lambda) - mixmod->lambda - lgamma(i+1.);
		}
	
	}
	

	for(i=0;i<G;i++){
		k = mixmod->whereis[i];
		GMM_recompute_marginal_likelihood_component( k, mixmod );
		//recompute_marginal_likelihood_component(mixmod->components[k],mixmod->alpha,mixmod->delta,mixmod->kappa,mixmod->gamma,mixmod->prior_mu,d);
		//printf("\nComponent %d has marg like = %.10f",i,mixmod->components[k]->log_prob);
	}
	
	
	
	
	free(tau_g);

}


void GMM_recompute_marginal_likelihood_component( int idx, struct mix_mod *mm )
{

	GMM_recompute_marginal_likelihood_component_0( mm->components[idx], mm );
	
	return;
	
}

void GMM_recompute_marginal_likelihood_component_0( struct component *cmp, struct mix_mod *mm )
{
	
	cmp->log_prob = GMM_return_marginal_likelihood_component( cmp, mm );
	
	return;
}

double  GMM_compute_marginal_likelihood_with_inclusion_in_component( double  *x, struct component *comp, struct mix_mod *mm )
{
	component_add_to_component( comp, x, 1 );
	
	double lp = GMM_return_marginal_likelihood_component( comp, mm );
	
	component_add_to_component( comp, x, -1 );
	
	return( lp );
}

double  GMM_compute_marginal_likelihood_with_inclusion_in_component_uni( double  x, struct component *comp, struct mix_mod *mm )
{
	component_add_to_component_uni( comp, x, 1 );
	
	double lp = GMM_return_marginal_likelihood_component( comp, mm );
	
	component_add_to_component_uni( comp, x, -1 );
	
	return( lp );
}


double GMM_return_marginal_likelihood_component( struct component *cmp, struct mix_mod *mm )
{

	int i, n = cmp->n_g, d = mm->d ;
	double lp, a, sq_norm=0., 
			*sum = cmp->sum,
			*mu = mm->prior_mu ;
	
	lp = lgamma( n +  mm->alpha ) + lgamma( 0.5*( n * d + mm->delta ) ) - 0.5* d * log( n + mm->kappa ) ;
	
	for( i=0; i<mm->d; i++ )
	{
		 a = sum[i] + mm->kappa * mu[i] ;
		 sq_norm += a * a ;
	}
	
	lp -= .5 * ( n * d + mm->delta ) * log( cmp->sum_squared_norm - sq_norm/( n + mm->kappa ) + mm->kappa * mm->xi2 + mm->gamma  );
	
	return( lp ) ; 

}




double log_normalizing_constant_model( int G, struct mix_mod *mm )
/*returns the log of the normalizing constant for a model with G components*/
{
	double z;
	
	z = 0.5 * G * mm->d * log( mm->kappa ) + lgamma( G * mm->alpha ) - G * lgamma(mm->alpha) - lgamma(mm->n+G*mm->alpha) + 0.5*G*mm->delta*log(mm->gamma) - G*lgamma(0.5*mm->delta);

	return(z); 	
}


void do_mixmod_analysis_one_sweep(struct results *pres,struct mix_mod *mixmod,int fix_G, int iter)
/*fixedG takes the value either TRUE or FALSE as defined in the macros*/
{


	int ej_case,maxgroups = mixmod->maxgroups;

	double pr_ej_G,pr_ej_Gm1,pr_ej_Gp1;
	
	update_allocations_with_gibbs(mixmod);
 
 	//Move M1 always has a low acceptance rate... so don't bother
  update_allocations_with_metropolis_move_1(mixmod,&(pres->accepted_m1),&(pres->proposed_m1));

	update_allocations_with_metropolis_move_2(mixmod,&(pres->accepted_m2),&(pres->proposed_m2));

 	update_allocations_with_metropolis_move_3(mixmod,&(pres->accepted_m3),&(pres->proposed_m3));


	//update no. groups every 10th iteration to improve mixing
	if(!fix_G /*&& iter%10 == 0*/  ){

		if(mixmod->G == 1){
			ej_case = 0;
		}else if(mixmod->G == maxgroups){
			ej_case = 1;
		}else if(mixmod->G == 2){
			ej_case = 2;
		}else if(mixmod->G == maxgroups-1){
			ej_case = 3;
		}else{
			ej_case = 4;
		}
				
		switch(ej_case)
		{
			case 0:
				pr_ej_G = 1.;
				pr_ej_Gp1 = .5;
				pr_ej_Gm1 = 0.;
			break;
			
			case 1:
				pr_ej_G = 0.;
				pr_ej_Gp1 = 0.;
				pr_ej_Gm1 = .5;
			break;
		
			case 2:
				pr_ej_G = .5;
				pr_ej_Gp1 = .5;
				pr_ej_Gm1 = 1.;
			break;
			
			case 3:
				pr_ej_G = 0.5;
				pr_ej_Gp1 = 0.;
				pr_ej_Gm1 = 0.5;
			break;
			
			case 4:
				pr_ej_G = 0.5;
				pr_ej_Gp1 = 0.5;
				pr_ej_Gm1 = 0.5;
			break;
		}
		
		
		if(runif(0.0,1.0) < pr_ej_G){
		
 			
 			if( 0 ){
 			// update_allocations_with_split_move(mixmod,&(pres->accepted_eject),&(pres->proposed_eject),pr_ej_G,pr_ej_Gp1);
 			}else{
 			update_allocations_with_ejection_move(mixmod,&(pres->accepted_eject),&(pres->proposed_eject),pr_ej_G,pr_ej_Gp1);
 			}
			
		}
		else
		{
			if( 0 ){
			// update_allocations_with_combine_move(mixmod,&(pres->accepted_absorb),&(pres->proposed_absorb),pr_ej_G,pr_ej_Gm1);
			}else{
 			update_allocations_with_absorb_move(mixmod,&(pres->accepted_absorb),&(pres->proposed_absorb),pr_ej_G,pr_ej_Gm1);
			}

		}
			
	}
			
		
	//only bother doing this if some of the hyperparameters to be updated
	if(mixmod->update_gamma || mixmod->update_kappa || mixmod->update_prior_mu)
	{	
		update_hyperparameters(mixmod);
	}
				
		
	//store in results
	/*pres->ngroups[0] = mixmod->G;
	int i;
	for(i=0;i<mixmod->n;i++)
	{
		pres->memberships[0][i] = mixmod->z[i];
	}*/

}


/* a util to permute a vector  */
void random_ranshuffle( int *a, int n )
{
	//randomly permute the first n elements of a using Fisher-Yates algorithm
	//  using RNG from R
	
	int i, j, x;
	
	for( i=n-1; i>0; i-- )
	{
		j = (int)( i * runif(0.0,1.0) ) ;
		x = a[i] ;
		a[i] = a[j] ;
		a[j] = x;
	}
	
	return;
}

