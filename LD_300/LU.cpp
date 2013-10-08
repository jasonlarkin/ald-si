/*                              LU.cpp                              */
/*                            12/02/2008                            */
/*********************************************************************
*    Program file that contains the LU factorization method and      *
*  matrix algorithms that depend upon LU factorization.              *
*********************************************************************/


/*DECLARE HEADERS*/
#include <cstdlib>
#include <cmath>
#define TINY 1.0e-20


/*DEFINE LU FACTORIZATION ROUTINE*/
int LU(int n, double **A, int *permute=NULL) {
  //DECLARE LOCAL VARIABLES
  int i, j, k, imax, d=1;
  double big, dum, sum, temp;
  double *vv=new double[n];


  //GET SCALE FACTORS AND STORE IN vv
  for(i=0;i<n;i++) {    //Loop over rows to get implicit scaling factor
    big = 0.0;
    for(j=0;j<n;j++) {if((temp=fabs(A[i][j]))>big) {big = temp;}}
    if(big<TINY) {return(0);}  //Singular or nearly singular matrix
    vv[i] = 1.0/big;
  }


  //FIND L AND U MATRICIES
  for(j=0;j<n;j++) {
    for(i=0;i<j;i++) {
      sum = A[i][j];
      for(k=0;k<i;k++) {sum -= A[i][k]*A[k][j];}
      A[i][j] = sum;
    }
    big = -1.0;            //For pivot element
    imax = 0;
    for(i=j;i<n;i++) {
      sum = A[i][j];
      for(k=0;k<j;k++) {sum -= A[i][k]*A[k][j];}
      A[i][j] = sum;
      if((dum=vv[i]*fabs(sum))>big) {big = dum;imax = i;}//Find best pivot
    }
    if(j!=imax) {         //Interchange rows
      for(k=0;k<n;k++) {
        dum = A[imax][k];
        A[imax][k] = A[j][k];
        A[j][k] = dum;
      }
      d = -d;
      vv[imax] = vv[j];//Does this need swaped instead?
    }
    if(permute!=NULL) {permute[j] = imax;}
    if(fabs(A[j][j])<TINY) {if(A[j][j]<0.0) {A[j][j] = -TINY;} else{A[j][j] = TINY;}}
    if(j!=n) {            //Divide by pivot element
      dum = 1.0/A[j][j];
      for(i=j+1;i<n;i++) {A[i][j] *= dum;}
    }
  }


  //DEALLOCATE MEMORY
  delete[] vv;  vv = NULL;

  return(d);                  //Gives the sign of the determinant
}


/*DEFINE SUBROUTINE LU_LinearEqs: FINDS b IN Ax = b*/
void LU_LinearEqs(int n, double **A, double *b, int *permute) {
  //DECLARE LOCAL VARIABLES
  int i, ii=-1, ip, j;
  double sum;

  //BACK SUBSTITUTE FROM LU RESULTS
  //Find Ly = b (y is stored in b)
  for(i=0;i<n;i++) {
    ip = permute[i];              //Find proper row as altered in LU
    sum = b[ip];
    b[ip] = b[i];
    if(ii>=0) {for(j=ii;j<=i-1;j++) {sum -= A[i][j]*b[j];}}
    else if(fabs(sum)>=TINY) {ii = i;}
    b[i] = sum;
  }
  //Find LUx = b -> LUx = Ly -> Ux = y (x is stored in b)
  for(i=n-1;i>=0;i--) {
    sum = b[i];
    for(j=i+1;j<n;j++) {sum -= A[i][j]*b[j];}
    b[i] = sum/A[i][i];
  }

  return;
}

/*DEFINE SUBROUTINE LU_Determinant: USES THE U FACTOR MATRIX TO COMPUTE DET(A)*/
double LU_Determinant(int n, double **A, double sign=1.0) {
  //DET(A) = sign*DET(L)*DET(U)
  for(int i=0;i<n;i++) {sign *= A[i][i];}
  return(sign);
}


/*DEFINE SUBROUTINE LU_Inverse: COMPUTES INVERSE OF A GIVEN LU FACTORS*/
void LU_Inverse(int n, double **A, int *permute) {
  //DECLARE LOCAL VARIABLES
  int i, j;
  double *col=new double[n];
  double **A_inv = new double*[n];
  A_inv[0] = new double[n*n];
  for(i=0;i<n;i++) {if(i<n-1) {A_inv[i+1] = A_inv[i] + n;}}
  //COMPUTE INVERSE
  for(j=0;j<n;j++) {
    for(i=0;i<n;i++) {col[i] = 0.0;}
    col[j] = 1.0;
    LU_LinearEqs(n, A, col, permute);
    for(i=0;i<n;i++) {A_inv[i][j] = col[i];}
  }
  for(i=0;i<n;i++) {for(j=0;j<n;j++) {A[i][j] = A_inv[i][j];}}
  //DEALLOCATE MEMORY AND RETURN
  delete[] col;  col = NULL;
  delete[] A_inv[0];  A_inv[0] = NULL;
  delete[] A_inv;  A_inv = NULL;
  return;
}
