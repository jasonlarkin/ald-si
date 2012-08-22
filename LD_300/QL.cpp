/*							              	QL.cpp		              						*/
/*							              06/25/2007		              					*/
/*********************************************************************
*    Subroutine that solves the eigenvalue problem for a tridiagonal *
*  matrix.  Adapted from numerical recipies: "QL algorithm with		   *
*  implicit shifts, to determine the eigenvalues and eigenvectors of *
*  a real, symmetric, tridiagonal matrix, or of a real, symmetric 	 *
*  matrix previously reduced by tred2 §11.2. On input, d[0..n-1]	 	 *
*  contains the diagonal elements of the tridiagonal matrix. On	  	 *
*  output, it returns the eigenvalues. The vector e[0..n-1] inputs   *
*  the subdiagonal elements of the tridiagonal matrix, with e[0]		 *
*  arbitrary. On output e is destroyed. When finding only the		     *
*  eigenvalues, several lines may be omitted, as noted in the		     *
*  comments. If the eigenvectors of a tridiagonal matrix are desired,*
*  the matrix z[0..n-1][0..n-1] is input as the identity matrix. If	 *
*  the eigenvectors of a matrix that has been reduced by tred2 are   *
*  required, then z is input as the matrix output by tred2. In either*
*  case, the kth column of z returns the normalized eigenvector		   *
*  corresponding to d[k]."											                     *
*********************************************************************/

/*HEADERS*/
#include <cmath>


/*DEFINE SYSTEM DEFINITIONS*/
//#define EIGENVECTOR


/*SUBROUTINES*/
inline double Sign(double a, double b) {
	double c;
	if(b>=0.0) {c = fabs(a);}
	else {c = -fabs(a);}
	return(c);
}

inline double Pythag(double a, double b) {
	double c;
	c = sqrt(a*a + b*b);
	return(c);
}


/*SUBROUTINE QL1*/
void QL1(double *d, double *e, int n, double ***z) {
	//VARIABLE DECLARATION
	int m,l,i;//,k=0;
	double s,r,p,g,f,dd,c,b;


	//QL ALGORITHM
	for (l=0;l<n;l++) {
		do {
			for (m=l;m<n-1;m++) { //Look for a single small subdiagonal element to split the matrix.
				dd=fabs(d[m])+fabs(d[m+1]);
				if (float(fabs(e[m])+dd)==float(dd)) {break;}//if ((float)(fabs(e[m])+dd) == dd) {break;}
			}
      if (m != l) {
				g=(d[l+1]-d[l])/(2.0*e[l]); //Form shift.
				r=Pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+Sign(r,g)); //This is dm - ks.
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) { //A plane rotation as in the original QL, followed by Givens rotations to restore tridiagonal form.
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=Pythag(f,g));
					if (r == 0.0) { //Recover from underflow.
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					/* Next loop can be omitted if eigenvectors not wanted*/
/*#if defined EIGENVECTOR
					for (k=0;k<n;k++) { //Form eigenvectors.
						f=z[k][i+1][1];
						z[k][i+1][1]=s*z[k][i][1]+c*f;
						z[k][i][1]=c*z[k][i][1]-s*f;
						f=z[k][i+1][0];
						z[k][i+1][0]=s*z[k][i][0]+c*f;
						z[k][i][0]=c*z[k][i][0]-s*f;
					}
#endif*/
				}
				if (r == 0.0 && i >= l) {continue;}
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while (m != l);
	}

	return;
}


/*SUBROUTINE QL2*/
void QL2(double *d, double *e, int n, double ***z) {
	//VARIABLE DECLARATION
	int m,l,i,k=0;
	double s,r,p,g,f,dd,c,b;


	//QL ALGORITHM
	for (l=0;l<n;l++) {
		do {
			for (m=l;m<n-1;m++) { //Look for a single small subdiagonal element to split the matrix.
				dd=fabs(d[m])+fabs(d[m+1]);
				if (float(fabs(e[m])+dd)==float(dd)) {break;}
			}
      if (m != l) {
				g=(d[l+1]-d[l])/(2.0*e[l]); //Form shift.
				r=Pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+Sign(r,g)); //This is dm - ks.
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) { //A plane rotation as in the original QL, followed by Givens rotations to restore tridiagonal form.
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=Pythag(f,g));
					if (r == 0.0) { //Recover from underflow.
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					/* Next loop can be omitted if eigenvectors not wanted*/
/*#if defined EIGENVECTOR*/
					for (k=0;k<n;k++) { //Form eigenvectors.
						f=z[k][i+1][1];
						z[k][i+1][1]=s*z[k][i][1]+c*f;
						z[k][i][1]=c*z[k][i][1]-s*f;
						f=z[k][i+1][0];
						z[k][i+1][0]=s*z[k][i][0]+c*f;
						z[k][i][0]=c*z[k][i][0]-s*f;
					}
/*#endif*/
				}
				if (r == 0.0 && i >= l) {continue;}
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while (m != l);
	}

	return;
}
