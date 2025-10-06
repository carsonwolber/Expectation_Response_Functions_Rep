

#include <math.h>
#include <mex.h>
#include <lapack.h>
#include <blas.h>
#include <stdlib.h>
#include <string.h>


/* Input Arguments */
#define	x1_in	prhs[0]
#define	x2_in	prhs[1]
#define	x3_in	prhs[2]
#define	x4_in	prhs[3]
#define	x5_in	prhs[4]
#define	x6_in	prhs[5]
#define x7_in   prhs[6]
#define x8_in   prhs[7]
#define	xx_in	prhs[8]

/* Output Arguments */
#define out0    plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

void insert_sortd(ptrdiff_t m, mwIndex *VV, double *II){
    int i,j;
    mwIndex x;
    double xi;
    for(i=1;i<m;i++){
        x = VV[i];
        xi = II[i];
        
        j = i;
        
        while(j>0&&VV[j-1]>x){
            VV[j] = VV[j-1];
            II[j] = II[j-1];
            j--;
        }
        VV[j] = x;
        II[j] = xi;
    }
    return;
} 

/* Prints an MxN matrix to Screen*/
void insert_sort(ptrdiff_t m, double *VV, mwIndex *II){
    int i,j;
    double x;
    mwIndex xi;
    for(i=1;i<m;i++){
        x = VV[i];
        xi = II[i];
        
        j = i;
        
        while(j>0&&VV[j-1]>x){
            VV[j] = VV[j-1];
            II[j] = II[j-1];
            j--;
        }
        VV[j] = x;
        II[j] = xi;
    }
    return;
}





void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )

{
    
    double *x1,*x2,*x3,*x4,*x5,*x6,*x7,*x8,*xx,*mu;

    double xcord[8];
    
    long ns,ii,idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8;
    long nx1,nx2, nx3, nx4, nx5, nx6, nx7, nx8;
    mwIndex II[8] = {0,1,2,3,4,5,6,7};
    mwIndex S[8]  = {0,0,0,0,0,0,0,0};
    ptrdiff_t one = 1;
    
    /*Check for proper number of arguments*/
    if (nrhs == 9)  { 
    }
    else{
        mexErrMsgTxt("9 inputs arguments required.");
    }
    
    /* Get Dimensions of Input Arguments*/
    nx1 = mxGetN(x1_in);
    nx2 = mxGetN(x2_in);
    nx3 = mxGetN(x3_in);
    nx4 = mxGetN(x4_in);
    nx5 = mxGetN(x5_in);
    nx6 = mxGetN(x6_in);
    nx7 = mxGetN(x7_in);
    nx8 = mxGetN(x8_in);
    ns  = mxGetN(xx_in);
    
    ptrdiff_t nx1ptr = nx1;
    ptrdiff_t nx2ptr = nx2;
    ptrdiff_t nx3ptr = nx3;
    ptrdiff_t nx4ptr = nx4;
    ptrdiff_t nx5ptr = nx5;
    ptrdiff_t nx6ptr = nx6;
    ptrdiff_t nx7ptr = nx7;
    ptrdiff_t nx8ptr = nx8;
    
    ptrdiff_t sixptr = 6;
    ptrdiff_t sevenptr  = 7;
    ptrdiff_t eightptr  = 8;
    ptrdiff_t nineptr  = 9;
    
    /*Create output argument*/
    mwSize mw_ns  = (mwSize)ns;
    mwSize mw_nnx = (mwSize)nx1*nx2*nx3*nx4*nx5*nx6*nx7*nx8;
    mwSize mw_nnz = (mwSize)ns*(8+1);
    
    double  *Qval;
    mwIndex *Qidx;      
    mwIndex *jcs;
        
    out0  = mxCreateSparse(mw_nnx,mw_ns,mw_nnz,mxREAL); 
    Qval  = mxGetPr(out0); 
    Qidx  = mxGetIr(out0); 
    jcs   = mxGetJc(out0);
    
    /* Assign pointers to the input arguments*/
    x1  = mxGetPr(x1_in);
    x2  = mxGetPr(x2_in);
    x3  = mxGetPr(x3_in);
    x4  = mxGetPr(x4_in);
    x5  = mxGetPr(x5_in);
    x6  = mxGetPr(x6_in);
    x7  = mxGetPr(x7_in);
    x8  = mxGetPr(x8_in);
    xx  = mxGetPr(xx_in);

    
  
    double x1ex[nx1];
    double x2ex[nx2];
    double x3ex[nx3];
    double x4ex[nx4];
    double x5ex[nx5];
    double x6ex[nx6];
    double x7ex[nx7];
    double x8ex[nx8];
    double xxtmp[8];
    
    /*Create versions of x and y that have infinity in the end*/
     memcpy(x1ex, x1, nx1*sizeof(double));
     memcpy(x2ex, x2, nx2*sizeof(double));
     memcpy(x3ex, x3, nx3*sizeof(double));
     memcpy(x4ex, x4, nx4*sizeof(double));
     memcpy(x5ex, x5, nx5*sizeof(double));
     memcpy(x6ex, x6, nx6*sizeof(double));
     memcpy(x7ex, x7, nx7*sizeof(double));
     memcpy(x8ex, x8, nx8*sizeof(double));
     
    x1ex[0] = -DBL_MAX;
    x1ex[nx1-1] = DBL_MAX;
    
    x2ex[0] = -DBL_MAX;
    x2ex[nx2-1] = DBL_MAX;
    
    x3ex[0] = -DBL_MAX;
    x3ex[nx3-1] = DBL_MAX;
      
    x4ex[0] = -DBL_MAX;
    x4ex[nx4-1] = DBL_MAX;
    
    x5ex[0] = -DBL_MAX;
    x5ex[nx5-1] = DBL_MAX;
    
    x6ex[0] = -DBL_MAX;
    x6ex[nx6-1] = DBL_MAX;
    
    x7ex[0] = -DBL_MAX;
    x7ex[nx7-1] = DBL_MAX;
    
    x8ex[0] = -DBL_MAX;
    x8ex[nx8-1] = DBL_MAX;
    
    for (ii=0;ii<ns;ii++){
        idx1 = 0;
        idx2 = 0;
        idx3 = 0;
        idx4 = 0;
        idx5 = 0;
        idx6 = 0;
        idx7 = 0;
        idx8 = 0;
        
        /*Copy over the current column of x*/
        memcpy(xxtmp, &xx[8*ii], 8*sizeof(double));
        
        while (x1ex[idx1]<=xxtmp[0]){
            idx1++;
        }
        while (x2ex[idx2]<=xxtmp[1]){
            idx2++;
        }
        while (x3ex[idx3]<=xxtmp[2]){
            idx3++;
        }
        while (x4ex[idx4]<=xxtmp[3]){
            idx4++;
        }
        while (x5ex[idx5]<=xxtmp[4]){
            idx5++;
        }
        while (x6ex[idx6]<=xxtmp[5]){
            idx6++;
        }
        while (x7ex[idx7]<=xxtmp[6]){
            idx7++;
        }
        while (x8ex[idx8]<=xxtmp[7]){
            idx8++;
        }
        
        
        /*xcord(jj) = (xx(jj)-x(idx1(jj)))/(x(idx1(jj)+1) - x(idx1(jj)));*/
        xcord[0] = (xxtmp[0] - x1[idx1-1])/(x1[idx1] - x1[idx1-1]);
        xcord[1] = (xxtmp[1] - x2[idx2-1])/(x2[idx2] - x2[idx2-1]);
        xcord[2] = (xxtmp[2] - x3[idx3-1])/(x3[idx3] - x3[idx3-1]);
        xcord[3] = (xxtmp[3] - x4[idx4-1])/(x4[idx4] - x4[idx4-1]);
        xcord[4] = (xxtmp[4] - x5[idx5-1])/(x5[idx5] - x5[idx5-1]);
        xcord[5] = (xxtmp[5] - x6[idx6-1])/(x6[idx6] - x6[idx6-1]);
        xcord[6] = (xxtmp[6] - x7[idx7-1])/(x7[idx7] - x7[idx7-1]);
        xcord[7] = (xxtmp[7] - x8[idx8-1])/(x8[idx8] - x8[idx8-1]);
        
        /*Sort the x cords*/
        S[0] = 0;
        S[1] = 0;
        S[2] = 0;
        S[3] = 0;
        S[4] = 0;
        S[5] = 0;
        S[6] = 0;
        S[7] = 0;
        
        II[0] = 0;
        II[1] = 1;
        II[2] = 2;
        II[3] = 3;
        II[4] = 4;
        II[5] = 5;
        II[6] = 6;
        II[7] = 7;
        
        insert_sort(eightptr,xcord,II);
        
        
        /*Compute the coordinate vectors*/
        Qidx[9*ii] = (idx1 + S[0]) + nx1*(idx2+S[1]) + nx1*nx2*(idx3+S[2]) + nx1*nx2*nx3*(idx4+S[3]) + nx1*nx2*nx3*nx4*(idx5+S[4]) + nx1*nx2*nx3*nx4*nx5*(idx6+S[5]) + nx1*nx2*nx3*nx4*nx5*nx6*(idx7+S[6]) + nx1*nx2*nx3*nx4*nx5*nx6*nx7*(idx8+S[7]);
        Qval[9*ii] = xcord[0];
        S[II[0]] = -1;
               
        Qidx[9*ii+1] = (idx1 + S[0]) + nx1*(idx2+S[1]) + nx1*nx2*(idx3+S[2]) + nx1*nx2*nx3*(idx4+S[3]) + nx1*nx2*nx3*nx4*(idx5+S[4]) + nx1*nx2*nx3*nx4*nx5*(idx6+S[5]) + nx1*nx2*nx3*nx4*nx5*nx6*(idx7+S[6])+ nx1*nx2*nx3*nx4*nx5*nx6*nx7*(idx8+S[7]);
        Qval[9*ii+1] = xcord[1]-xcord[0];
        S[II[1]] = -1;

        Qidx[9*ii+2] = (idx1 + S[0]) + nx1*(idx2+S[1]) + nx1*nx2*(idx3+S[2]) + nx1*nx2*nx3*(idx4+S[3]) + nx1*nx2*nx3*nx4*(idx5+S[4]) + nx1*nx2*nx3*nx4*nx5*(idx6+S[5]) + nx1*nx2*nx3*nx4*nx5*nx6*(idx7+S[6])+ nx1*nx2*nx3*nx4*nx5*nx6*nx7*(idx8+S[7]);
        Qval[9*ii+2] = xcord[2]-xcord[1];
        S[II[2]] = -1;
        
        Qidx[9*ii+3] = (idx1 + S[0]) + nx1*(idx2+S[1]) + nx1*nx2*(idx3+S[2]) + nx1*nx2*nx3*(idx4+S[3]) + nx1*nx2*nx3*nx4*(idx5+S[4]) + nx1*nx2*nx3*nx4*nx5*(idx6+S[5])+ nx1*nx2*nx3*nx4*nx5*nx6*(idx7+S[6])+ nx1*nx2*nx3*nx4*nx5*nx6*nx7*(idx8+S[7]);
        Qval[9*ii+3] = xcord[3]-xcord[2];
        S[II[3]] = -1;
        
        Qidx[9*ii+4] = (idx1 + S[0]) + nx1*(idx2+S[1]) + nx1*nx2*(idx3+S[2]) + nx1*nx2*nx3*(idx4+S[3]) + nx1*nx2*nx3*nx4*(idx5+S[4]) + nx1*nx2*nx3*nx4*nx5*(idx6+S[5])+ nx1*nx2*nx3*nx4*nx5*nx6*(idx7+S[6])+ nx1*nx2*nx3*nx4*nx5*nx6*nx7*(idx8+S[7]);
        Qval[9*ii+4] = xcord[4]-xcord[3];
        S[II[4]] = -1;
        
        Qidx[9*ii+5] = (idx1 + S[0]) + nx1*(idx2+S[1]) + nx1*nx2*(idx3+S[2]) + nx1*nx2*nx3*(idx4+S[3]) + nx1*nx2*nx3*nx4*(idx5+S[4]) + nx1*nx2*nx3*nx4*nx5*(idx6+S[5])+ nx1*nx2*nx3*nx4*nx5*nx6*(idx7+S[6])+ nx1*nx2*nx3*nx4*nx5*nx6*nx7*(idx8+S[7]);
        Qval[9*ii+5] = xcord[5]-xcord[4];
        S[II[5]] = -1;
        
        Qidx[9*ii+6] = (idx1 + S[0]) + nx1*(idx2+S[1]) + nx1*nx2*(idx3+S[2]) + nx1*nx2*nx3*(idx4+S[3]) + nx1*nx2*nx3*nx4*(idx5+S[4]) + nx1*nx2*nx3*nx4*nx5*(idx6+S[5])+ nx1*nx2*nx3*nx4*nx5*nx6*(idx7+S[6])+ nx1*nx2*nx3*nx4*nx5*nx6*nx7*(idx8+S[7]);
        Qval[9*ii+6] = xcord[6]-xcord[5];  
        S[II[6]] = -1;
        
        Qidx[9*ii+7] = (idx1 + S[0]) + nx1*(idx2+S[1]) + nx1*nx2*(idx3+S[2]) + nx1*nx2*nx3*(idx4+S[3]) + nx1*nx2*nx3*nx4*(idx5+S[4]) + nx1*nx2*nx3*nx4*nx5*(idx6+S[5])+ nx1*nx2*nx3*nx4*nx5*nx6*(idx7+S[6])+ nx1*nx2*nx3*nx4*nx5*nx6*nx7*(idx8+S[7]);
        Qval[9*ii+7] = xcord[7]-xcord[6];  
        S[II[7]] = -1;
        
        Qidx[9*ii+8] = (idx1 + S[0]) + nx1*(idx2+S[1]) + nx1*nx2*(idx3+S[2]) + nx1*nx2*nx3*(idx4+S[3]) + nx1*nx2*nx3*nx4*(idx5+S[4]) + nx1*nx2*nx3*nx4*nx5*(idx6+S[5])+ nx1*nx2*nx3*nx4*nx5*nx6*(idx7+S[6])+ nx1*nx2*nx3*nx4*nx5*nx6*nx7*(idx8+S[7]);
        Qval[9*ii+8] = 1-xcord[7]; 
        
        /*Sort the coordinate vectors*/
        insert_sortd(nineptr,&Qidx[9*ii],&Qval[9*ii]);

        /*row indexes*/
        jcs[ii+1] = 9*(ii+1);
    }
   return;
}

