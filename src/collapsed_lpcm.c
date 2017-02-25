#include "required_libs.h"
#include "GaussianMixtureModel.h"
#include "NetworkLib.h"
#include "cat.h"


//function decls

void collapsed_get_expected_Y( int *n, int *d, int *sample, int *dir, double *ExY, double *beta, double *latentpos );

void Relabel(int *n_obs,int *n_sample,int *n_groups,int *labels_in,int *labels_out) ;

void collapsed_lpcm( int *Y, int *nnodes, int *dimlatent, int *ncovariates, double *covariates, int *directed, int *maxgroups, int *initgroups, int *sample, int *burn, int *interval, int *modelsearch, double *hyparams, double *log_prior_groups, double *prparams, double *initparams, double *initialpositions, int *return_allocations, double *return_latent_positions, double *return_abundance, int *return_ngroups, double *return_llike, double *return_theta, double *return_delta, double *return_kappa, double *accrt_latent_positions,double *accrt_abundance, double *accrt_metmoves, double *accrt_ejectabsorb, double *accrt_theta, int *nthread, int *update_gamma, int *update_kappa, int *npilot, int *store_sparse, int *adapt, int *adapt_int, int *bradley_terry, int *verbose );

//void collapsed_lpcm_MLE_distances( int *Y, int *nnodes, int *directed, double *dist, double *beta );


void collapsed_lpcm( int *Y, int *nnodes, int *dimlatent, int *ncovariates, double *covariates, int *directed, int *maxgroups, int *initgroups, int *sample, int *burn, int *interval, int *modelsearch, double *hyparams, double *log_prior_groups, double *prparams, double *initparams, double *initialpositions, int *return_allocations, double *return_latent_positions, double *return_abundance, int *return_ngroups, double *return_llike, double *return_theta, double *return_delta, double *return_kappa, double *accrt_latent_positions,double *accrt_abundance, double *accrt_metmoves, double *accrt_ejectabsorb, double *accrt_theta, int *nthread, int *update_gamma, int *update_kappa, int *npilot, int *store_sparse, int *adapt, int *adapt_int, int *bradley_terry, int *verbose )
{
	
	int i,j,it,actori,s,fixedG=1-*modelsearch;
	double tinit, lps;
	
	int nburn = *burn, niterations = (*sample) * (*interval) + (*burn) ;
	
	struct results *results;

	struct resy *yres;
	
	results = (struct results *)malloc(sizeof(struct results));
	
	struct network *nw = network_create( *nnodes, *dimlatent, *ncovariates, *directed, *maxgroups, *initgroups, *bradley_terry  );	
	
	yres = (struct resy *)malloc(sizeof(struct resy));
	
	double *prior_hparams,*tc;
	
	prior_hparams = calloc(4,sizeof(double));
	tc = calloc(2,sizeof(double));
	
	//setstartime(tc);
	
	//cross reference with other code
	double sigmab = prparams[0], sigmaz = prparams[1], 
		   bcur = initparams[0], thetainit = initparams[1];
	
	double *thetainitpar = calloc( 2, sizeof(double) );
	
	if(*ncovariates == 0){
		tinit = 0.;
	}else{
		tinit = initparams[1];
	}
	
	prior_hparams[0] = hyparams[0]; prior_hparams[1] = hyparams[1];
	prior_hparams[2] = 0.; prior_hparams[3] = 1.; //the last two are for covariates.. ignore for now. 
	
	network_initialize( nw, Y, bcur, thetainitpar, prior_hparams,  sigmab, sigmaz, thetainitpar, initialpositions, log_prior_groups );
	
	allocate_results( results, niterations, nburn, *nnodes );
	
	initresy( yres, *ncovariates );
	
	if(ncovariates > 0){
	//	put_covariates(covariates,network);
	}

	int *Gcount;
	Gcount = calloc(*maxgroups,sizeof(int));
	
	nw->pmix->gamma = hyparams[2];
	nw->pmix->delta = hyparams[3];
	nw->pmix->alpha = hyparams[4];
	nw->pmix->kappa = 1./hyparams[5];


	if( *ncovariates > 0 ){

	}
	
	for(j=0;j<nw->pmix->d;j++){
		nw->pmix->prior_mu[j] = 0.;
	}
	nw->pmix->xi2 = 0.;
	
	nw->pmix->update_gamma = FALSE;
	
	initialize_simple( nw->pmix, *initgroups );

	nw->llike = llike_full( nw );

	nw->pmix->update_kappa = *update_kappa;
	nw->pmix->update_gamma = *update_gamma;
	nw->pmix->update_prior_mu = FALSE;
	
	set_prior_hyperparameters_x( nw->pmix, hyparams[2], 1./hyparams[5] ) ;
	
	//enable random sampling from R
	GetRNGstate();

	int jp; 
	
	double accrt, del, Del;
	
	int *order = calloc( *nnodes, sizeof(int) );
	
	if( *verbose ) Rprintf("\n\t Starting MCMC iterations... ");
	
	for( it=0; it<niterations + *npilot; it++){
	
		R_CheckUserInterrupt(); 
	
		if( !(*bradley_terry) )	betaupdate( nw, yres, it, nburn, 1.);
		
		if(*ncovariates > 0){
			//do update of the covariates
			//for(s=0;s<*ncovariates;s++)
			//	thetaupdate(network,yres,it,nburn,1.,s);
		}
		
		//update of actor positions 
		
		for( i=0; i<*nnodes; i++ ) order[i] = i;
		random_ranshuffle( order, nw->n );
		
		for( actori=0; actori<nw->n; actori++ )
		{
			zupdatemh( nw, yres, order[actori], it, nburn, 1. );
		}
		
		//adaptive proposals
		
		if( *adapt && it < nburn && ( it % (*adapt_int) == 0 ) )
		{
		
			//adapt the proposal variances here
			
			del =  1./sqrt( it ) ;
			Del = del < .01 ? del : .01 ;
			
			//beta
			accrt = (double) yres->accepted_beta / yres->proposed_beta ;
			
			if( accrt > .234 ) nw->sigmab *= exp( Del ) ; else nw->sigmab *= exp( -Del ) ;
			
			//z update
			
			accrt = (double) yres->accepted_z / yres->proposed_z ;
			
			if( accrt > .234 ) nw->sigmaz *= exp( Del ) ; else nw->sigmaz *= exp( -Del ) ;
		
		}
		
		if( it > (int)(.5*nburn) ) fixedG = 1-*modelsearch;  else fixedG = TRUE ; 
		
		do_mixmod_analysis_one_sweep( results, nw->pmix, fixedG, it );
		
		if( *verbose && it == *npilot ) Rprintf("\n\t Starting burn in iterations...");
		if( *verbose && it == nburn  + *npilot - 1 ) Rprintf("\n\t Finished burn in iterations...");
		
		if( it > nburn + *npilot - 1 ){
			
			jp = it-nburn-*npilot+1;
		
			if(jp%(*interval) == 0){
			
				if( *verbose && (jp/(*interval)-1)%1000 == 0 && jp/(*interval)-1 > 0 ) 
					Rprintf("\n\t Finished sample %d... ", jp/(*interval)-1 );
			
				//record beta and latent positions
				
				return_ngroups[jp/(*interval)-1] = nw->pmix->G;
				
				if( !(*store_sparse) )
				{
				
					return_abundance[jp/(*interval)-1] = nw->beta;
					return_delta[ jp/(*interval)-1] = nw->pmix->gamma;
					return_kappa[ jp/(*interval)-1] = nw->pmix->kappa;
					
					for(i=0;i<nw->n;i++){
						return_allocations[(jp/(*interval)-1)*nw->n + i] = nw->pmix->z[i];
						//printf(" %d ",network->pmix->z[i]);
					}
					//printf("\n");
					for(i=0;i<nw->n;i++){
						for(j=0;j<nw->pmix->d;j++){
							if( nw->pmix->d > 1 ) lps = nw->pmix->Y[i][j]; else lps = nw->pmix->y_uni[i]; 
							return_latent_positions[(jp/(*interval)-1)*( (*dimlatent)*nw->n) + (*dimlatent)*i + j] = lps;
						}
					}
				
					if(*ncovariates>0){
						//store the values of the covariates
						for(i=0;i<*ncovariates;i++){
							return_theta[(jp/(*interval)-1)*(*ncovariates) + i] = nw->theta[i];
						}
					}
					
					return_llike[jp/(*interval)-1] = nw->llike;
					//Rprintf("\n Likcur %lf ", network->likcur );
				
				}
				
			
			}
		
		}
		
		if( it == nburn + *npilot - 1 )
		{
		 	//reset all of the acceptance counters after pilot and burnin
		 	yres->accepted_z = 0;
		 	yres->proposed_z = 0;
			results->accepted_eject = 0;
			results->proposed_eject = 0;
			results->accepted_absorb = 0;
			results->proposed_absorb = 0;
			results->accepted_m1 = 0;
			results->proposed_m1 = 0;
			results->accepted_m2 = 0;
			results->proposed_m2 = 0;
			results->accepted_m3 = 0;
			results->proposed_m3 = 0;
		}
	
	}
	
	if( *verbose ) Rprintf("\n\t Finished MCMC...") ; 
	
	PutRNGstate();
	
	//store the acc rates
	*accrt_latent_positions =100.*(double)yres->accepted_z/((double)yres->proposed_z);
	*accrt_abundance = 100.*(double)yres->accepted_beta/((double)yres->proposed_beta);
	if(*ncovariates>0){
		for(i=0;i<*ncovariates;i++){
			//accrt_theta[i] = 100.*(double)yres->accepted_theta[i]/(*niterations);
		}
	}
	accrt_ejectabsorb[0] = 100.*(double)results->accepted_eject/results->proposed_eject;
	accrt_ejectabsorb[1] = 100.*(double)results->accepted_absorb/results->proposed_absorb;
	accrt_metmoves[0] = 100.*(double)results->accepted_m1/results->proposed_m1;
	accrt_metmoves[1] = 100.*(double)results->accepted_m2/results->proposed_m2;
	accrt_metmoves[2] = 100.*(double)results->accepted_m3/results->proposed_m3;

	prparams[0] = nw->sigmab * nw->sigmab;
	prparams[1] = nw->sigmaz * nw->sigmaz;
	
	free_results(results,niterations,nburn);
	free(results);
	free(yres);
	network_destroy( nw );
	free(prior_hparams);
	free( thetainitpar );
	free( order );
	free(tc);
	
	return;
}

//label switching algorithm

void Relabel(int *n_obs,int *n_sample,int *n_groups,int *labels_in,int *labels_out)
{

	int n,g,**raw,**relab,**summary,
		i,j,k,t,N,T,**cost,*lab;
	
	/*n and no. of groups*/
	n=n_obs[0];
	g=n_groups[0];
	/*N = number of samples to undo label switching for*/
	N=n_sample[0];

	T=0;

	/*allocate memory and initialize*/
	raw = imatrix(1,N,1,n);
	relab = imatrix(1,N,1,n);
	summary = imatrix(1,g,1,n);
	cost = imatrix(1,g,1,g+1);
	lab = ivector(1,g);


	/*copy from labels_in to raw_r*/
	for(i=0;i<N;i++){
		for(j=0;j<n;j++){
			raw[i+1][j+1] = labels_in[ i + j*N ];
		}	
	}


	for(i=1;i<g+1;i++){
		for(j=1;j<n+1;j++){
			summary[i][j]=0;
		}
	}

	/*use the first allocation in raw as the first labelling*/

	for(i=1;i<n+1;i++){
		summary[raw[1][i]][i]+=1;
	}

	for(i=1;i<n+1;i++){
		relab[1][i] = raw[1][i];
	}
	

	for(t=2;t<N+1;t++){

	/*row t*/
	/*compute the cost matrix*/
	for(i=1;i<g+1;i++){
	
		for(j=1;j<g+1;j++){
	
			cost[i][j]=0;
		
			for(k=1;k<n+1;k++){
		
				if(raw[t][k]==j){
					cost[i][j]+=summary[i][k];
				}
			
			} 
		
			cost[i][j]=n*(t-1)-cost[i][j];
		
		}
		cost[i][g+1]=0;
	
		}

		T=0;
		
		assct(g,cost,lab,&T);

		/*relabel based on output from assct*/
		/*update the summary matrix*/
		for(i=1;i<n+1;i++){
			relab[t][i]=lab[raw[t][i]];
			summary[lab[raw[t][i]]][i]+=1;
		}

	}


	for(i=0;i<N;i++){
		for(j=0;j<n;j++){
			labels_out[ i + j*N ] = relab[i+1][j+1];
		}	
	}	
	
	free_imatrix(raw,1,N,1,n);
	free_imatrix(relab,1,N,1,n);
	free_imatrix(summary,1,g,1,n);
	free_imatrix(cost,1,g,1,g+1);
	free_ivector(lab,1,g);

	return;
}

void collapsed_get_expected_Y( int *n, int *d, int *sample, int *dir, double *ExY, double *beta, double *latentpos )
{

	int iter, i, j, k;
	double *lpos = calloc( (*n)*(*d), sizeof(double) ), eta;

	//don't worry about no. of groups as not needed for this part
	struct network *nw = network_create( *n, *d, 0, *dir, 2, 2, 0 );		

	//cycle through samples and compute the expected prob
	for( iter=0; iter<*sample; iter++ )
	{
		
		k = iter * (*n) * (*d );
		for( i=0; i< (*n)*(*d); i++ ) lpos[i] = latentpos[ k + i ];
		
		put_latentpositions( lpos, nw );
		
		nw->beta = beta[ iter ];
		
		for( i=0; i<*n; i++ )
		{
			for( j=0; j<*n; j++ ) 
			{
				eta = nw->beta - nw->dist[i][j] ;
				ExY[ i + (*n) * j ] += 1./( 1. + exp( - 1. * eta ) );
			}
		}
		
	}	
	
	for( i=0; i<*n; i++ ) 
	{
		ExY[ i + (*n) * i ] = 0.;
		for( j=0; j<*n; j++ )
		{
			ExY[ i + (*n) * j ] /= (*sample) ; 
		}
	}

	free( lpos );

	return;
}
