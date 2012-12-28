/*							               Eigen.cpp		              					*/
/*							              06/25/2007	              						*/
/*********************************************************************
*    Subroutine that computes the eigenvalues and vectors of a	  	 *
*  matrix.													                            		 *
*********************************************************************/

/*DEFINE HEADERS*/
#include <cstdlib>


/*DEFINE SYSTEM DEFINITIONS*/
//#define EIGENVECTOR


/*SUBROUTINE Eigen1*/
void Eigen1(int size, double *Freq, double ***D) {
	//DECLARE LOCAL VARIABLES
	int i, j, n;
	void Householder1(double ***a, int n, double *d, double *e);
	void QL1(double *d, double *e, int size, double ***a);
	double *e, *d;
	double swap, trial;


	//ALLOCATE MEMORY
	e = (double*) calloc(size, sizeof(double));
	d = (double*) calloc(size, sizeof(double));


	//CALL SUBROUTINES
	Householder1(D, size, d, e);
	QL1(d, e, size, D);


	//ORDER EIGENVALUES FROM LARGEST TO SMALLEST
	for(i=0;i<size;i++) {
		swap = d[i];
		n = i;
		for(j=i+1;j<size;j++) {
			trial = d[j];
			if(swap<trial) {swap=trial;n=j;}
		}
		if(n!=i) {
			d[n] = d[i];
			d[i] = swap;
/*#if defined EIGENVECTOR
			for(j=0;j<size;j++) {
				swap = D[j][i][0];
				D[j][i][0] = D[j][n][0];
				D[j][n][0] = swap;
				swap = D[j][i][1];
				D[j][i][1] = D[j][n][1];
				D[j][n][1] = swap;
			}
#endif*/
		}
		Freq[i] = d[i];
	}


	//FREE DYNAMIC ARRAYS
	free(e);
	free(d);

	return;
}


/*SUBROUTINE Eigen*/
void Eigen2(int size, double *Freq, double ***D) {
	//DECLARE LOCAL VARIABLES
	int i, j, n;
	void Householder2(double ***a, int n, double *d, double *e);
	void QL2(double *d, double *e, int size, double ***a);
	double *e, *d;
	double swap, trial;


  //ALLOCATE MEMORY
	e = (double*) calloc(size, sizeof(double));
	d = (double*) calloc(size, sizeof(double));


  //CALL SUBROUTINES
	Householder2(D, size, d, e);
	QL2(d, e, size, D);


	//ORDER EIGENVALUES FROM LARGEST TO SMALLEST
	for(i=0;i<size;i++) {
		swap = d[i];
		n = i;
		for(j=i+1;j<size;j++) {
			trial = d[j];
			if(swap<trial) {swap=trial;n=j;}
		}
		if(n!=i) {
			d[n] = d[i];
			d[i] = swap;
/*#if defined EIGENVECTOR*/
			for(j=0;j<size;j++) {
				swap = D[j][i][0];
				D[j][i][0] = D[j][n][0];
				D[j][n][0] = swap;
				swap = D[j][i][1];
				D[j][i][1] = D[j][n][1];
				D[j][n][1] = swap;
			}
/*#endif*/
		}
		Freq[i] = d[i];
	}


	//FREE DYNAMIC ARRAYS
	free(e);
	free(d);

	return;
}
