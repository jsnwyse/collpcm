/*C functions to do evaluations for the latent position cluster model

	Authors: Caitriona Ryan (University of Limerick, Ireland)  
			 some slight modifications to Catriona's  functions made by Jason Wyse
	
	Corresponding: Jason Wyse,
			School of Computer Science and Statistics,
			Trinity College Dublin,
			Dublin 2, Ireland.
			email: wyseja@tcd.ie
			
Last modified: Fri 14 Mar 2014 13:00:01 GMT  */ 

#include "NetworkLib.h"
#define debug FALSE


/*new functions added by Jason for the resnet structure*/


void put_network(int *Y,struct network *nw)
{

	int i,j;
	
	for(i=0;i<nw->n;i++)
	{
		for(j=0;j<nw->n;j++)
		{
			nw->y[i][j] = Y[i + j*(nw->n)]; //remember column major format from R
			nw->y_transpose[j][i] = nw->y[i][j];
		}
	}

	return;
}

void put_latentpositions( double *z, struct network *nw )
{

	struct mix_mod *mixmod = nw->pmix ;

	int i, j, k, n = mixmod->n, d = mixmod->d;
	
	double **X, *x, a;
	
	if( nw->d > 1 ) X= mixmod->Y; else x = mixmod->y_uni;

	for(i=0;i<n;i++)
	{
		if( nw->d > 1 )
		{
			for(k=0;k<d;k++) X[i][k] = z[i + n*k];
			
		}
		else
		{
			x[i] = z[i] ; 
		}
	}
	
	//for( i=0; i<n; i++ ) dist_update( nw, i) ;
	
	//compute the distances
	if( nw->modty == 0 )
	{
		for( i=0; i<n-1; i++ )
		{
			for( j=i+1; j<n; j++ )
			{
				a = 0.;
				if( d > 1 )
				{
					for( k=0; k<d; k++ ) a += ( X[i][k] - X[j][k] ) * ( X[i][k] - X[j][k] ) ;
				}
				else
				{
					a += ( x[i] - x[j] ) * ( x[i] - x[j] ) ;
				}
				nw->dist[j][i] = nw->dist[i][j] = sqrt( a ); 
			}
		}
	}else{
		
		for( i=0; i<n; i++ )
		{
			for( j=0; j<n; j++ ) nw->dist[i][j] = x[i] - x[j] ;
		}
	
	}
	
	return;
}


void put_covariates(double *x,struct network *nw)
{

	int ncov = nw->p,n = nw->n,i,j;
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<ncov;j++)
			nw->xcovs[i][j] = x[i + j*n];
	}

	return;

}

/*this function allocates and returns a pointer to a ynet structure...*/

struct network *network_create( int n , int d, int p, int dir, int maxG, int initG, int modty )
{
	int i;
	
	struct network *nw = (struct network *)malloc( sizeof(struct network) ) ;
	
	nw->n = n;
	nw->p = p;
	nw->d = d;
	nw->dir = dir;
	nw->modty = modty;
	
	//network and covariates
	nw->y  = (int **)calloc( n, sizeof(int*)  );
	nw->y_transpose = (int **)calloc( n, sizeof(int*) );
	nw->dist =  (double **)calloc( n, sizeof(double*) );
	if( p > 0 )
	{
		nw->theta = (double *)calloc( p, sizeof(double) );
		nw->sigmatheta = (double *)calloc( p,  sizeof(double) );
		nw->xcovs = (double **)calloc( n, sizeof(double*) );
	}
	
	for( i=0; i<n; i++ )
	{
		nw->y[i] = (int *)calloc( n, sizeof(int) );
		nw->y_transpose[i] = (int *)calloc( n, sizeof(int) );
		nw->dist[i] = (double *)calloc( n, sizeof(double ) );
		if( p > 0 ) nw->xcovs[i] = (double *)calloc( n, sizeof(double) );
	}
	
	nw->llike = -DBL_MAX;
	
	nw->pmix = mixmod_create( n, d, maxG, initG );
	
	return( nw );
}

void network_destroy( struct network *nw )
{
	int i, n = nw->n;
	
	for( i=0; i<n; i++ )
	{
		free( nw->y[i] );
		free( nw->y_transpose[i] );
		free( nw->dist[i] );
		if( nw->p > 0 ) free( nw->xcovs[i] );
	}
	
	free( nw->y );
	free( nw->y_transpose );
	free( nw->dist );
	if( nw->p > 0 ) 
	{
		free( nw->xcovs );
		free( nw->theta );
		free( nw->sigmatheta );
	}
	
	mixmod_destroy( nw->pmix );
	
	free( nw );
	
	return;
}

void network_initialize( struct network *nw, int *Y, double beta, double *theta, double *hyper_params, double sigmab, double sigmaz, double *sigmatheta, double *initialpositions, double *log_prior_groups )
{
	int k, p = nw->p ;
	
	put_network( Y, nw );
	
	nw->beta = beta ;
	
	if( p > 0 ) 
	{
		for( k=0; k<p; k++ ) 
		{
			nw->theta[k] = theta[k] ;
			nw->sigmatheta[k] = sigmatheta[k] ;
		}
	} 
	
	nw->sigmab = sigmab;
	nw->sigmaz = sigmaz;
	
	nw->xi = hyper_params[0] ;
	nw->psi = hyper_params[1] ;
	nw->rho = hyper_params[2] ;
	nw->zeta = hyper_params[3] ;
	
	put_latentpositions( initialpositions, nw );
	
	for( k=0 ; k<nw->pmix->maxgroups+1; k++ ) nw->pmix->log_prior_G[k] = log_prior_groups[k];

	return;
	
}


void initresy	(struct resy *presy, int ncovs){
  presy->accepted_beta=0; 	/*counter for betaupdate acceptance rate*/	
  presy->proposed_beta=0;
  presy-> accepted_z=0; 	/*counter for zupdate acceptance rate*/	
  presy->proposed_z = 0;
  int i;

	presy->accepted_theta=(int *)calloc(ncovs,sizeof(int));
	for (i=0;i<ncovs;i++){
	presy->accepted_theta[i]=0;
	 }
  
}

void dist_update( struct network *nw, int i )
{
	int j, k, n = nw->n, d = nw->pmix->d ;
	double **X, *x, **dist = nw->dist, a ;
	
	if( d > 1 )
		X = nw->pmix->Y ;
	else
		x = nw->pmix->y_uni ;

	//compute the distances
	if( nw->modty == 0 )
	{
		for( j=0; j<n; j++ )
		{
			a = 0.;
		
			if( d > 1 )
			{
				for( k=0; k<d; k++ ) a += ( X[i][k] - X[j][k] ) * ( X[i][k] - X[j][k] ) ;
			}
			else
			{
				a += ( x[i] - x[j] ) * ( x[i] - x[j] ) ;
			}
		
			dist[j][i] = dist[i][j] = sqrt( a );
		}
	}else{
		
		for( j=0; j<n; j++ )
		{
			dist[i][j] = x[i] - x[j] ;
		}
	
	}
	
	return;
}


void zupdatemh( struct network * nw, struct resy * presy, int i, int itnum, int burnin, double c )
{
	int k;
	
	presy->proposed_z += 1;
	
	struct mix_mod *mix = nw->pmix ;
	
	int g = nw->pmix->z[i], d = mix->d ;

	double llik_curr = llike_node( nw, i ), 
			 *xdelta,
			 xdelta_,
			 *x ;

	//if( d > 1 ) xdelta = (double *)calloc( d, sizeof(double) ) ;
	
	//Rprintf("\n llikcurr = %.5f", llik_curr ); 

	if( d > 1 ) 
	{
		x = mix->Y[i]; 
		xdelta = (double *)calloc( d, sizeof(double) );
	}

	struct component *comp = mix->components[ mix->whereis[g] ] ;
	
	if( d > 1 )
		component_add_to_component( comp, x, -1 );
	else
		component_add_to_component_uni( comp, mix->y_uni[i], -1 );
	
	if( d > 1 )
	{
		for( k=0; k<d; k++ ) 
		{
			xdelta[k] = rnorm( 0., nw->sigmaz );
			x[k] += xdelta[k] ;
		}
	}
	else
	{
		xdelta_ = rnorm( 0, nw->sigmaz );
		mix->y_uni[i] += xdelta_ ;
	}
	
	dist_update( nw, i );
	
	if( d > 1 )
		component_add_to_component( comp, x, 1 );
	else
		component_add_to_component_uni( comp, mix->y_uni[i], 1 );
	
	double llik_prop = llike_node( nw, i ),
			 lmarg_prop = GMM_return_marginal_likelihood_component( comp, mix ) ;
			 
	/*Rprintf("\n lmargcurr = %.5f", comp->log_prob ); 
	Rprintf("\n llikprop = %.5f", llik_prop );
	Rprintf("\n xdelta_ = %.5f", xdelta_ );
	Rprintf("\n x = %.5f", mix->y_uni[i] ); 
	Rprintf("\n lmargprop = %.5f", lmarg_prop );*/
	
	double lratio = llik_prop + lmarg_prop - ( llik_curr + comp->log_prob ) ;
	
	if( lratio > log( runif(0.,1.) ) )
	{
		presy->accepted_z += 1;
		comp->log_prob = lmarg_prop;
		nw->llike += (llik_prop - llik_curr);
	}
	else
	{
		if( d > 1 )
		{	
			component_add_to_component( comp, x, -1 );
			for( k=0; k<d; k++ ) x[k] -= xdelta[k] ;
			component_add_to_component( comp, x, 1 );
		}
		else
		{
			component_add_to_component_uni( comp, mix->y_uni[i], -1 );
			mix->y_uni[i] -= xdelta_ ;
			component_add_to_component_uni( comp, mix->y_uni[i], 1 );
		}
		dist_update( nw, i );
	}

	if( d > 1 ) free( xdelta );
	
	return;
}




void betaupdate(struct network * nw, struct resy * presy, int itnum, int burnin, double c )
{
 	 // update for the intercept parameter

	presy->proposed_beta += 1;

	double llike_prop, llike_curr = nw->llike, beta_curr = nw->beta ;
  
  //llike_curr = llike_full( py );
  
  	nw->beta += rnorm( 0.0, nw->sigmab );
  
  // recalculate likelihood with this new value of beta
  llike_prop = llike_full( nw ); 
  
  double lratio = llike_prop + dnorm( nw->beta - nw->xi, 0.0, sqrt(nw->psi), 1 ) 
  					   - ( llike_curr + dnorm( beta_curr - nw->xi, 0.0, sqrt( nw->psi ) , 1 ) )  ; 

  if( lratio > log( runif( 0.0, 1.0 ) ) )
  { 
    //accept proposed value. 
    nw->llike = llike_prop;
    presy->accepted_beta += 1; 
  }
  else 
  { 
  	//change back to original value of beta
    nw->beta = beta_curr;
  }
  
  return;
}


double get_eta( double b, int d, double *x_1, double *x_2 )
{
	int k;
	double u = 0., a;
	for( k=0; k<d; k++ )
	{
		 a = x_1[k] - x_2[k];
		 u += a * a;
	}
	return( b - sqrt(u) );
}


//it appears that using the distance matrix is a little bit faster actually

double llike_node(struct network * nw, int i)
{
  int    j, *y_i = nw->y[i], *yt_i = nw->y_transpose[i] ;
  
  double loglike = 0., eta, b = nw->beta, /***X = py->pmix->Y ,*/*d_i = nw->dist[i], prob;
  
  if( nw->modty == 0 )
  {
  
	  if( nw->dir )
	  {
	  		for( j=0; j<nw->n; j++ )
	  		{
	  			eta = /*get_eta( b, d, X[i], X[j] );*/ b - d_i[j] ; 
	  			loglike += ( y_i[j] + yt_i[j] ) * eta - 2.0 * log( 1.0 + exp( eta ) ) ;
	  		}
	  		loglike += 2.0 * log( 1.0 + exp( b ) ) ;
	  		return( loglike ); 
	  }
	  else
	  {
	  		for( j=0; j<nw->n; j++ )
	  		{
	  			eta = /*get_eta( b, d, X[i], X[j] );*/ b - d_i[j] ; 
	  			loglike += y_i[j] * eta - log( 1.0 + exp( eta ) ) ;
	  		}
	  		loglike += log( 1.0 + exp( b ) ) ;
	  		return( loglike );
	  }
  
  }else{
  	
  		for( j=0; j<nw->n; j++ )
  		{
  			prob = 1./( 1. + exp( -d_i[j] ) ) ;
  			if( j != i ) loglike += dbinom( (double) y_i[j], (double)( y_i[j] + yt_i[j] ), prob, 1 ) ; 
  		}
  		return( loglike );
  }

}


double llike_full(struct network * nw )
{

  int i, j , n = nw->n, *y_i, *yt_i;
  			
  double loglike = 0.0, eta, b = nw->beta, /***X = py->pmix->Y,*/ *d_i, prob ;
  
  if( nw->modty == 0 )
  {
  
	  if( nw->dir )
	  {
	  		for( i=0; i<n-1; i++ )
	  		{
	  			y_i = nw->y[i];
	  			yt_i = nw->y_transpose[i];
	  			d_i = nw->dist[i] ;
	  			for( j=i+1; j<n; j++ )
	  			{
	  				eta = /*get_eta( b, d, X[i], X[j] );*/  b - d_i[j];//X[i][j]; //				
	  				loglike += ( y_i[j] + yt_i[j] ) * eta - 2.0 * log( 1.0 + exp( eta ) );
	  			}
	  		}
	  		return( loglike ); 
	  }
	  else
	  {
	  		for( i=0; i<n-1; i++ )
	  		{
	  			y_i = nw->y[i];
	  			d_i = nw->dist[i];
	  			for( j=i+1; j<n; j++ )
	  			{
	  				eta = /*get_eta( b, d, X[i], X[j] );*/ b - d_i[j];//X[i][j]; 
	  				loglike += y_i[j] * eta - log( 1.0 + exp( eta ) );
	  			}
	  		}
	  		return( loglike );
	  }
  
  }else{
  
  		for( i=0; i<n-1; i++ )
	  	{
	  		y_i = nw->y[i];
	  		yt_i = nw->y_transpose[i];
	  		d_i = nw->dist[i] ;
	  		for( j=i+1; j<n; j++ )
	  		{
	  			//eta = /*get_eta( b, d, X[i], X[j] );*/  b - d_i[j];//X[i][j]; //	
	  			prob = 1. / ( 1. + exp( -d_i[j] ) ) ;			
	  			loglike += dbinom( (double) y_i[j], (double)( y_i[j] + yt_i[j] ), prob, 1 ) ; 
	  		}
	  	}
	  	return( loglike );
  }
  
}

/********************** Computing the initial distances ***************************/





