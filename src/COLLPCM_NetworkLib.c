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



/*void put_network(int *Y,struct network *nw)
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
 }*/

/*void put_latentpositions( double *z, struct network *nw )
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
 }*/


/*void put_covariates(double *x,struct network *nw)
 {
 
 int ncov = nw->p,n = nw->n,i,j;
 
 for(i=0;i<n;i++)
 {
 for(j=0;j<ncov;j++)
 nw->xcovs[i][j] = x[i + j*n];
 }
 
 return;
 
 }*/

//this function allocates and returns a pointer to a ynet structure...

struct network *COLLPCM_network_create( int n , int d, int p, int dir, int maxG, int initG )
{
  int i;
  
  struct network *nw = (struct network *)malloc( sizeof(struct network) ) ;
  
  nw->n = n;
  nw->p = p;
  nw->d = d;
  nw->dir = dir;
  
  nw->dist =  (double *)calloc( n * n, sizeof(double) );
  if( p > 0 )
  {
    nw->theta = (double *)calloc( p, sizeof(double) );
    nw->sigmatheta = (double *)calloc( p,  sizeof(double) );
  }
  
  nw->llike = -DBL_MAX;
  
  nw->pmix = COLLPCM_allocate_mixmod( n, d, maxG, initG );
  
  return( nw );
}

void COLLPCM_network_destroy( struct network *nw )
{
  int i, n = nw->n;
  
  free( nw->dist );
  
  if( nw->p > 0 ) 
  {
    free( nw->theta );
    free( nw->sigmatheta );
  }
  
  //COLLPCM_free_mixmod( nw->pmix );
  
  return;
}

void COLLPCM_network_initialize( struct network *nw, int *Y, double beta, double *theta, double *hyper_params, double sigmab, double sigmaz, double *sigmatheta, double *initialpositions, double *log_prior_groups )
{
  int k, p = nw->p ;
  
  nw->y = Y; // point directly to the network data
  
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
  
  //nw->pmix->y = initialpositions;
  
  nw->pmix->log_prior_G = log_prior_groups;
  
  return;
  
}


void COLLPCM_initresy	(struct resy *presy, int ncovs)
{
  presy->accepted_beta = 0; 	/*counter for betaupdate acceptance rate*/	
  presy->proposed_beta = 0;
  presy->accepted_z = 0; 	/*counter for zupdate acceptance rate*/	
  presy->proposed_z = 0;
  int i;

  //presy->accepted_theta=(int *)calloc(ncovs,sizeof(int));
  //for (i=0;i<ncovs;i++)
    //presy->accepted_theta[i] = 0;

return;
}

void COLLPCM_dist_update( struct network *nw, int i )
{
  int j, k, n = nw->n, d = nw->pmix->d ;
  double *x, *xi, *disti, *dist0, a, r ;
  
  x = nw->pmix->y ;
  xi = nw->pmix->y + i*d ;
  dist0 = nw->dist + i ;
  disti = nw->dist + i*n ;
   
  //compute the distances against xi
  for( j=0; j<n; j++ )
  {
    a = 0.0;
    for( k=0; k<d; k++ )
    {
      r = xi[k] - x[k];
      a += r*r; 
    }
    disti[j] = sqrt(a);
    dist0[0] = disti[j];
    dist0 += n ;
    x += d ;
  }
  
  return;
}


void COLLPCM_zupdatemh( struct network * nw, struct resy * presy, int i, int itnum, int burnin, double c )
{
  int k;
  
  presy->proposed_z += 1;
  
  struct mix_mod *mix = nw->pmix ;
  
  int g = nw->pmix->z[i], d = mix->d, n = nw->n ;
  
  double llik_curr = COLLPCM_llike_node( nw, i ), 
    *xdelta,
    *x ;
  
  x = mix->y + i*n; 
  xdelta = (double *)calloc( d, sizeof(double) );
  
  struct component *comp = mix->components[ mix->whereis[g] ] ;
  
  COLLPCM_add_to_component( comp, x, -1 );
  
  for( k=0; k<d; k++ ) 
  {
    xdelta[k] = rnorm( 0., nw->sigmaz );
    x[k] += xdelta[k] ;
  }
  
  COLLPCM_dist_update( nw, i );
  
  COLLPCM_add_to_component( comp, x, 1 );
  
  double llik_prop = COLLPCM_llike_node( nw, i ),
    lmarg_prop = COLLPCM_return_log_marginal_likelihood_component( comp, mix ) ;
  
  double lratio = llik_prop + lmarg_prop - ( llik_curr + comp->log_prob ) ;
  
  if( lratio > log( runif(0.0,1.0) ) )
  {
    presy->accepted_z += 1;
    comp->log_prob = lmarg_prop;
    nw->llike += llik_prop - llik_curr;
  }
  else
  {
    COLLPCM_add_to_component( comp, x, -1 );
    for( k=0; k<d; k++ ) x[k] -= xdelta[k] ;
    COLLPCM_add_to_component( comp, x, 1 );
    COLLPCM_dist_update( nw, i );
  }
  
  free( xdelta );
  
  return;
}




void COLLPCM_betaupdate(struct network * nw, struct resy * presy, int itnum, int burnin, double c )
{
  // update for the intercept parameter
  
  presy->proposed_beta += 1;
  
  double llike_prop, llike_curr = nw->llike, beta_curr = nw->beta ;
  
  //llike_curr = llike_full( py );
  
  nw->beta += rnorm( 0.0, nw->sigmab );
  
  // recalculate likelihood with this new value of beta
  llike_prop = COLLPCM_llike_full( nw ); 
  
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


/*double get_eta( double b, int d, double *x_1, double *x_2 )
{
  int k;
  double u = 0., a;
  for( k=0; k<d; k++ )
  {
    a = x_1[k] - x_2[k];
    u += a * a;
  }
  return( b - sqrt(u) );
}*/


// compute log-likelihood contribution for one node
double COLLPCM_llike_node(struct network * nw, int i)
{
  int j, n = nw->n, *y_i = nw->y + i*n, *y_j = nw->y + i ;
  
  double loglike = 0.0, eta, b = nw->beta, *d_i = nw->dist + i*n, prob;
  
  if( nw->dir )
  {
    for( j=0; j<n; j++ )
    {
      eta = b - d_i[j] ; 
      loglike += ( y_i[j] + y_j[0] ) * eta - 2.0 * log( 1.0 + exp( eta ) ) ;
      y_j += n ;
    }
    loglike += 2.0 * log( 1.0 + exp( b ) ) ;
    return( loglike ); 
  }
  else
  {
    for( j=0; j<n; j++ )
    {
      eta = b - d_i[j] ; 
      loglike += y_i[j] * eta - log( 1.0 + exp( eta ) ) ;
    }
    loglike += log( 1.0 + exp( b ) ) ;
    return( loglike );
  }
  
}


double COLLPCM_llike_full(struct network * nw )
{
  
  int i, j , n = nw->n, *y_i = nw->y, *y_j;
  
  double loglike = 0.0, eta, b = nw->beta, *d_i = nw->dist, prob ;
  
  if( nw->dir )
  {
    for( i=0; i<n-1; i++ )
    {
      //y_i = nw->y[i];
      //yt_i = nw->y_transpose[i];
      //d_i = nw->dist[i] ;
      y_j = nw->y + i*n + i; // point to the i-th column
      
      for( j=1; j<n-i; j++)
      {
        y_j += n;
        eta = b - d_i[j];
        loglike += ( y_i[j] + y_j[0] ) * eta - 2.0 * log( 1.0 + exp( eta ) );
      }
      
      y_i += n + 1; 
      d_i += n + 1;
      
      /*for( j=i+1; j<n; j++ )
       {
       y_j += n ;
       eta = b - d_i[j];
       loglike += ( y_i[j] + y_j[0] ) * eta - 2.0 * log( 1.0 + exp( eta ) );
       }
       y_i += n ;
       d_i += n ;*/
    }
    return( loglike ); 
  }
  else
  {
    for( i=0; i<n-1; i++ )
    {
      //y_i = nw->y[i];
      //d_i = nw->dist[i];
      for( j=1; j<n-i; j++ )
      {
        eta = b - d_i[j];
        loglike += y_i[j] * eta - log( 1.0 + exp( eta ) );
      }
      y_i += n + 1 ;
      d_i += n + 1 ;
    }
    return( loglike );
  }
  
}





