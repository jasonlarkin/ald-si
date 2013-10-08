/*                           Velocity.cpp                           */
/*                            12/02/2008                            */
/*********************************************************************
*    Subroutine that computes the phonon group velocity v_g=(dw/dk). *
*********************************************************************/

/*DEFINE HEADERS*/
#include "LDCode.h"
#include <cmath>


#define SMALL 0.0001


/*SUBROUTINE Velocity*/
void Velocity(PD_HARMONIC **FC, double *k, double *Freq, double **V, double ***D) {
	//DECLARE LOCAL VARIABLES
	int j, dim;		  					    //Counter
	int *permute = new int[DIM];  //
	double **A;                   //
	double dk;		      	        //Step UC_DOF
	double *k2=new double[DIM];		//Use new array so that k is preserved
	double *Freq_0, *Freq_1;		  //Frequency array at k-dk (V used for k+dk)
	double *f0, *f1;				      //Pointers to frequency arrays
	void Dynamical_Matrix(PD_HARMONIC **FC, double *k, double ***D, int dim=-1);//Computes dynamical matrix
	void Eigen1(int UC_DOF, double *Freq, double ***D);//Computes eigenvalues/vectors of D
	int LU(int n, double **A, int *permute=NULL);
	void LU_LinearEqs(int n, double **A, double *b, int *permute);


	//ALLOCATE MEMORY
	Freq_0 = new double[UC_DOF];
	Freq_1 = new double[UC_DOF];
	A = new double*[DIM];
	A[0] = new double[DIM*DIM];
	for(j=0;j<DIM;j++) {
	  if(j<DIM-1) {A[j+1] = A[j] + DIM;}
	  for(dim=0;dim<DIM;dim++) {A[j][dim] = Lattice->b[dim][j];}
	}
	LU(DIM, A, permute);

	for(j=0;j<UC_DOF;j++) {Freq[j] = Freq[j]*1.0e12*Parameter->time;}//Non-dimensionalize Freq
	for(dim=0;dim<DIM;dim++) {
	  f0 = Freq_0;
	  f1 = Freq_1;

	  //COMPUTE FREQUENCY AT k-dk
	  for(j=0;j<DIM;j++) {k2[j] = 0.0;}
	  k2[dim] = SMALL;
	  LU_LinearEqs(DIM, A, k2, permute);
	  for(j=0;j<DIM;j++) {k2[j] = k[j] - k2[j];}//k-dk
	  dk = 2.0*Parameter->pi*SMALL;
/*	  dk = 0.0;
	  b_mag = 0.0;
	  for(j=0;j<DIM;j++) {
	    k2[j] = k[j] - SMALL*Lattice->a[dim][j]*Lattice->b[dim][j];
	    dk += Lattice->a[dim][j]*Lattice->b[dim][j];
	    b_mag += Lattice->b[dim][j];
    }
    dk *= 2.0*Parameter->pi*SMALL*b_mag;*/
	  Dynamical_Matrix(FC, k2, D);	    //Compute Dynamical Matrix
	  Eigen1(UC_DOF, Freq_0, D);			//Compute Eigenvalues only

	  //COMPUTE FREQUENCY AT k+dk
	  if(k[dim]<0.499999) {//If not on a boundary
      if(fabs(k[dim])>SMALL) {      //If not at k=0
        dk *= 2.0;					        //Central difference
        //for(j=0;j<DIM;j++) {k2[j] = k[j] + SMALL*Lattice->a[dim][j]*Lattice->b[dim][j];}
        for(j=0;j<DIM;j++) {k2[j] = k[j]+k[j]-k2[j];}//k+dk
        //k2[dof] += SMALL*SMALL;		  //k+dk
        Dynamical_Matrix(FC, k2, D);    //Compute Dynamical Matrix
        Eigen1(UC_DOF, Freq_1, D);	//Compute Eigenvalues only
      }
      else {f1=Freq_0;f0=Freq;}		  //Swap pointers for forward difference
	  }
	  else {f1=Freq;}					        //Adjust pointers for backward difference
	  if(f0!=Freq) {
      for(j=0;j<UC_DOF;j++) {
        if(f0[j]<0.0) {f0[j] = 0.0;}
        else {f0[j] = sqrt(f0[j]);}
      }
    }
    if(f1!=Freq) {
        for(j=0;j<UC_DOF;j++) {
        if(f1[j]<0.0) {f1[j] = 0.0;}
        else {f1[j] = sqrt(f1[j]);}
      }
	  }

	  //EVALUATE DERIVATIVE
	  for(j=0;j<UC_DOF;j++) {V[j][dim] = (f1[j]-f0[j])/dk;}

	}


	//DEALLOCATE MEMORY
	delete[] permute; permute=NULL;
	delete[] A[0];    A[0]   =NULL;
	delete[] A;       A      =NULL;
	delete[] k2;      k2     =NULL;
	delete[] Freq_0;  Freq_0 =NULL;
	delete[] Freq_1;  Freq_1 =NULL;

	return;
}

#undef SMALL
