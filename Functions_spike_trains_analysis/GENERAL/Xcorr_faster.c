  /*  Program in .C  to calculate xcorr faster !!!   (normalized so that xcorr(x,x)(tau= 0) = 1 ) 

   X_corr_faster (x,y,maxlags)   
   computes E(x(n)y(n+m)), m = 0 ... max_lags
 
   [ccg_x_y, time] = Xcorr_faster(x,y,max_lags)
    
*/     

 /*   % Ici on veut calculer E(x(n), y(n+m)) , Hey Jude !   xcorrfast(x,y,maxlags) =   xcorr(y,x,maxlags,'coeffs') */
 /*                             Be aware to put x and y as columns  !!! */
#include <math.h>
#include <matrix.h>
#include <mex.h>   

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

     mxArray *x_in_m, *y_in_m, *time_out_m, *cnorm_cccg_E_1_E2_main_out_m;
    const mwSize *dims;
    double   *x,      *y,      *time,       *cnorm_cccg_E_1_E2_main;
    int       N,       i,       M_L ,        n_lag,  n_x  ;
    double  xsqr_sum, ysqr_sum , temp       ;
     xsqr_sum = ysqr_sum = temp =   0; 
    
    x_in_m   = mxDuplicateArray(prhs[0]);
    y_in_m   = mxDuplicateArray(prhs[1]);
    M_L      = mxGetScalar(prhs[2]);
 
    dims = mxGetDimensions(prhs[0]);
    N = (int)dims[0];  
 
    x        = mxGetPr(x_in_m);
    y        = mxGetPr(y_in_m);
    
      
    time_out_m = plhs[1]         = mxCreateDoubleMatrix(2*M_L+1,1,mxREAL);
    time                         = mxGetPr(time_out_m);
  
    cnorm_cccg_E_1_E2_main_out_m = plhs[0] = mxCreateDoubleMatrix(2*M_L+1,1,mxREAL);
    cnorm_cccg_E_1_E2_main       = mxGetPr(cnorm_cccg_E_1_E2_main_out_m);  
       
    
   /*                                                         // essayer demain
	plhs[0] = mxCreateNumericMatrix(2*M_L+1, 1, mxUINT32_CLASS, mxREAL);
	cnorm_cccg_E_1_E2_main = (double *) mxGetPr(plhs[0]);*/
    
      
    
    for(i=0;i< 2*M_L+1 ;i++)
    {  time[i ] = - M_L + i; }

    /* Forward */        
   for(n_lag = 0; n_lag < M_L +1; n_lag ++)
     {
        temp = 0;    
        for( n_x = 0; n_x <  N - n_lag ; n_x ++)  { temp = temp + x[n_x]*y[n_x + n_lag];}
          
           cnorm_cccg_E_1_E2_main[M_L + n_lag] = temp/(1.0*(N - n_lag));
     } 
    
    /* Backwards */
        for(n_lag = 1; n_lag < M_L +1; n_lag ++    )   {                                
           temp = 0;
           for(n_x = n_lag ; n_x < N ; n_x ++    )   {temp = temp + x[n_x ]*y[n_x - n_lag];}

           cnorm_cccg_E_1_E2_main[M_L  - n_lag] =  temp/(1.0*(N - n_lag)); 
        }
 
   /* Normalize */  
       for(n_x = 0; n_x < N ; n_x ++)
                                       {
                                              xsqr_sum = xsqr_sum + x[n_x]*x[n_x];
                                              ysqr_sum = ysqr_sum + y[n_x]*y[n_x];
                                       }
                                              xsqr_sum = xsqr_sum/(1.0*N);
                                              ysqr_sum = ysqr_sum/(1.0*N);

       temp      = sqrt(xsqr_sum*ysqr_sum);
       if ( temp == 0) {
                temp = 1;
                       }
       for(n_lag = 0; n_lag < 2*M_L+1; n_lag ++)  { cnorm_cccg_E_1_E2_main[n_lag] = cnorm_cccg_E_1_E2_main[n_lag]/temp; }
    
  
    
    return;
}
