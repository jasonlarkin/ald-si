/*						              Householder.cpp	              					*/
/*							              06/25/2007				              			*/
/*********************************************************************
*    Subroutine that reduces a hermetian matrix to tridiagonal form	 *
*  via Householder transformations with unitary matricies.  Original *
*  routine from numerical recipies for a real symmetric matrix		   *
*  adapted via prescription described by Lehoucq and used in Lapack	 *
*  for the complex case.										                      	 *
*********************************************************************/

/*HEADERS*/
#include <cmath>


/*DEFINE SYSTEM DEFINITIONS*/
//#define EIGENVECTOR


/*SUBROUTINE Householder*/
void Householder1(double ***a, int n, double *d, double *e) {
	//VARIABLE DECLARATION
	int i, j, k, l = 1;
	double v;
	double zeta_re, zeta_im;
	double sigma_re, sigma_im;
	double scale;
	double c_1, c_2;
	int *flag = new int[n];


	//REDUCE a MATRIX
	for(i=0;i<n-1;i++,l++) {  //Goto n-1 not just n-2 to ensure that matrix is real
		scale = fabs(a[l][i][1]);
		for(j=l+1;j<n;j++) {scale += fabs(a[j][i][0]) + fabs(a[j][i][1]);}
		if(scale<1.0e-20) {
			flag[i] = 0;
			d[i] = a[i][i][0];
			e[i] = a[l][i][0];
		}
		else {
			flag[i] = 1;
			scale += fabs(a[l][i][0]);
			v = 0.0;
			for(j=l;j<n;j++) {
				a[j][i][0] /= scale;
				a[j][i][1] /= scale;
				v += a[j][i][0]*a[j][i][0] + a[j][i][1]*a[j][i][1];
			}
			//zeta = x[i+1] (l=i+1)
			zeta_re = a[l][i][0];
			zeta_im = a[l][i][1];
			//v = +/- Euclidiean norm of x (sqrt(x*x'))
			v = sqrt(v);
			if(zeta_re<0.0) {v = -v;}

			c_1 = zeta_re + v;
			c_2 = c_1*c_1 + zeta_im*zeta_im;
			//sigma' = (zeta' + v)/v
			sigma_re = (zeta_re+v)/v;
			sigma_im = -zeta_im/v;
			for(j=l+1;j<n;j++) {
				//store w[j] in a[i][j] = w = (x + v*e)/(zeta + v)
				a[i][j][0] = (a[j][i][0]*c_1 + a[j][i][1]*zeta_im)/c_2;
				a[i][j][1] = (a[j][i][1]*c_1 - a[j][i][0]*zeta_im)/c_2;
				//store xx[j] in a[j][i] = xx = sigma*w
				a[j][i][0] = sigma_re*a[i][j][0] - sigma_im*a[i][j][1];
				a[j][i][1] = sigma_im*a[i][j][0] + sigma_re*a[i][j][1];
			}
			a[i][l][0] = 1.0;
			a[i][l][1] = 0.0;
			a[l][i][0] = sigma_re;
			a[l][i][1] = sigma_im;

			zeta_re = zeta_im = 0.0;
			for(j=l;j<n;j++) {
				//d(Re)+e(Im) = y = A*w
				d[j] = e[j] = 0.0;
				for(k=l;k<n;k++) {
					d[j] += a[j][k][0]*a[i][k][0] - a[j][k][1]*a[i][k][1];
					e[j] += a[j][k][1]*a[i][k][0] + a[j][k][0]*a[i][k][1];
				}
				//zeta = w'*y = w'*A*w
				zeta_re = zeta_re + a[i][j][0]*d[j] + a[i][j][1]*e[j];
				zeta_im = zeta_im - a[i][j][1]*d[j] + a[i][j][0]*e[j];
			}
			//Update a matrix A = A - Z' - Z + (ww * xx')
			for(j=l;j<n;j++) {
				//sigma = ww = zeta*xx, Z = y*xx'
				sigma_re = zeta_re*a[j][i][0] - zeta_im*a[j][i][1];
				sigma_im = zeta_im*a[j][i][0] + zeta_re*a[j][i][1];
				for(k=l;k<=j;k++) {
					a[j][k][0] = a[j][k][0] + (sigma_re*a[k][i][0] + sigma_im*a[k][i][1])
						- (d[j]*a[k][i][0] + e[j]*a[k][i][1])
						- (d[k]*a[j][i][0] + e[k]*a[j][i][1]);
					a[j][k][1] = a[j][k][1] + (sigma_im*a[k][i][0] - sigma_re*a[k][i][1])
						- (e[j]*a[k][i][0] - d[j]*a[k][i][1])
						+ (e[k]*a[j][i][0] - d[k]*a[j][i][1]);
					a[k][j][0] = a[j][k][0];
					a[k][j][1] = -a[j][k][1];
				}
			}
			e[i] = -v*scale;
			d[i] = a[i][i][0];
		}
	}
	d[n-1] = a[n-1][n-1][0];


	//FORM PRODUCT OF TRANSFORMATIONS
/*#if defined EIGENVECTOR
	l = n;
	flag[n-1] = 0;
	//Q_i = (I-sigma*w*w')*Q_i-1 = (Q_i-1)-w*xx'*(Q_i-1)
	for(i=n-1;i>=0;i--,l--) {
		if(flag[i]) {
			for(j=l;j<n;j++) {
				//sigma = xx'[k]*Q[k][j]
				sigma_re = sigma_im = 0.0;
				for(k=l;k<n;k++) {
					sigma_re += a[k][i][0]*a[k][j][0] + a[k][i][1]*a[k][j][1];
					sigma_im += a[k][i][0]*a[k][j][1] - a[k][i][1]*a[k][j][0];
				}
				//Q_i[k][j] = Q_i-1[k][j] - w[k]*sigma[j]
				for(k=l;k<n;k++) {
					a[k][j][0] -= sigma_re*a[i][k][0] - sigma_im*a[i][k][1];
					a[k][j][1] -= sigma_re*a[i][k][1] + sigma_im*a[i][k][0];
				}
			}
		}
		a[i][i][0] = 1.0;
		a[i][i][1] = 0.0;
		for(j=l;j<n;j++) {a[i][j][0] = a[j][i][0] = a[i][j][1] = a[j][i][1] = 0.0;}
	}
	for(i=1;i<n;i++) {
		for(j=1;j<n;j++) {
			if(fabs(a[i][j][0])<1.0e-8) {a[i][j][0] = 0.0;}
			if(fabs(a[i][j][1])<1.0e-8) {a[i][j][1] = 0.0;}
		}
	}
#endif*/
	delete[] flag;

	return;
}


/*SUBROUTINE Householder2*/
void Householder2(double ***a, int n, double *d, double *e) {
	//VARIABLE DECLARATION
	int i, j, k, l = 1;
	double v;
	double zeta_re, zeta_im;
	double sigma_re, sigma_im;
	double scale;
	double c_1, c_2;
	int *flag = new int[n];


	//REDUCE a MATRIX
	for(i=0;i<n-1;i++,l++) {  //Goto n-1 not just n-2 to ensure that matrix is real
		scale = fabs(a[l][i][1]);
		for(j=l+1;j<n;j++) {scale += fabs(a[j][i][0]) + fabs(a[j][i][1]);}
		if(scale<1.0e-20) {
			flag[i] = 0;
			d[i] = a[i][i][0];
			e[i] = a[l][i][0];
		}
		else {
			flag[i] = 1;
			scale += fabs(a[l][i][0]);
			v = 0.0;
			for(j=l;j<n;j++) {
				a[j][i][0] /= scale;
				a[j][i][1] /= scale;
				v += a[j][i][0]*a[j][i][0] + a[j][i][1]*a[j][i][1];
			}
			//zeta = x[i+1] (l=i+1)
			zeta_re = a[l][i][0];
			zeta_im = a[l][i][1];
			//v = +/- Euclidiean norm of x (sqrt(x*x'))
			v = sqrt(v);
			if(zeta_re<0.0) {v = -v;}

			c_1 = zeta_re + v;
			c_2 = c_1*c_1 + zeta_im*zeta_im;
			//sigma' = (zeta' + v)/v
			sigma_re = (zeta_re+v)/v;
			sigma_im = -zeta_im/v;
			for(j=l+1;j<n;j++) {
				//store w[j] in a[i][j] = w = (x + v*e)/(zeta + v)
				a[i][j][0] = (a[j][i][0]*c_1 + a[j][i][1]*zeta_im)/c_2;
				a[i][j][1] = (a[j][i][1]*c_1 - a[j][i][0]*zeta_im)/c_2;
				//store xx[j] in a[j][i] = xx = sigma*w
				a[j][i][0] = sigma_re*a[i][j][0] - sigma_im*a[i][j][1];
				a[j][i][1] = sigma_im*a[i][j][0] + sigma_re*a[i][j][1];
			}
			a[i][l][0] = 1.0;
			a[i][l][1] = 0.0;
			a[l][i][0] = sigma_re;
			a[l][i][1] = sigma_im;

			zeta_re = zeta_im = 0.0;
			for(j=l;j<n;j++) {
				//d(Re)+e(Im) = y = A*w
				d[j] = e[j] = 0.0;
				for(k=l;k<n;k++) {
					d[j] += a[j][k][0]*a[i][k][0] - a[j][k][1]*a[i][k][1];
					e[j] += a[j][k][1]*a[i][k][0] + a[j][k][0]*a[i][k][1];
				}
				//zeta = w'*y = w'*A*w
				zeta_re = zeta_re + a[i][j][0]*d[j] + a[i][j][1]*e[j];
				zeta_im = zeta_im - a[i][j][1]*d[j] + a[i][j][0]*e[j];
			}
			//Update a matrix A = A - Z' - Z + (ww * xx')
			for(j=l;j<n;j++) {
				//sigma = ww = zeta*xx, Z = y*xx'
				sigma_re = zeta_re*a[j][i][0] - zeta_im*a[j][i][1];
				sigma_im = zeta_im*a[j][i][0] + zeta_re*a[j][i][1];
				for(k=l;k<=j;k++) {
					a[j][k][0] = a[j][k][0] + (sigma_re*a[k][i][0] + sigma_im*a[k][i][1])
						- (d[j]*a[k][i][0] + e[j]*a[k][i][1])
						- (d[k]*a[j][i][0] + e[k]*a[j][i][1]);
					a[j][k][1] = a[j][k][1] + (sigma_im*a[k][i][0] - sigma_re*a[k][i][1])
						- (e[j]*a[k][i][0] - d[j]*a[k][i][1])
						+ (e[k]*a[j][i][0] - d[k]*a[j][i][1]);
					a[k][j][0] = a[j][k][0];
					a[k][j][1] = -a[j][k][1];
				}
			}
			e[i] = -v*scale;
			d[i] = a[i][i][0];
		}
	}
	d[n-1] = a[n-1][n-1][0];


	//FORM PRODUCT OF TRANSFORMATIONS
/*#if defined EIGENVECTOR*/
	l = n;
	flag[n-1] = 0;
	//Q_i = (I-sigma*w*w')*Q_i-1 = (Q_i-1)-w*xx'*(Q_i-1)
	for(i=n-1;i>=0;i--,l--) {
		if(flag[i]) {
			for(j=l;j<n;j++) {
				//sigma = xx'[k]*Q[k][j]
				sigma_re = sigma_im = 0.0;
				for(k=l;k<n;k++) {
					sigma_re += a[k][i][0]*a[k][j][0] + a[k][i][1]*a[k][j][1];
					sigma_im += a[k][i][0]*a[k][j][1] - a[k][i][1]*a[k][j][0];
				}
				//Q_i[k][j] = Q_i-1[k][j] - w[k]*sigma[j]
				for(k=l;k<n;k++) {
					a[k][j][0] -= sigma_re*a[i][k][0] - sigma_im*a[i][k][1];
					a[k][j][1] -= sigma_re*a[i][k][1] + sigma_im*a[i][k][0];
				}
			}
		}
		a[i][i][0] = 1.0;
		a[i][i][1] = 0.0;
		for(j=l;j<n;j++) {a[i][j][0] = a[j][i][0] = a[i][j][1] = a[j][i][1] = 0.0;}
	}
/*	for(i=1;i<n;i++) {
		for(j=1;j<n;j++) {
			if(fabs(a[i][j][0])<1.0e-10) {a[i][j][0] = 0.0;}
			if(fabs(a[i][j][1])<1.0e-10) {a[i][j][1] = 0.0;}
		}
	}*/
/*#endif*/
	delete[] flag;

	return;
}
