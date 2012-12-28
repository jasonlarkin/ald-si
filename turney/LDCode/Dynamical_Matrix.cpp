/*						           Dynamical_Matrix.cpp			            			*/
/*							              12/02/2008		              					*/
/*********************************************************************
*    Subroutine that computes the dynamical matrix.	        				 *
*********************************************************************/

/*DEFINE HEADERS*/
#include "LDCode.h"
#include <cmath>

#define DM_ZERO 1.0e-14


/*DECLARE SUBROUTINES*/


/*SUBROUTINE Dynamical_Matrix*/
void Dynamical_Matrix(PD_HARMONIC **FC, double *k, double ***D, int dim=-1) {
	//DECLARE LOCAL VARIABLES
	int b0, b1;	                //Counters
	int d0, d1;				          //Cartesian counters
	int n_i, n_j;				        //Atom number
	double *ptr_fc;             //Pointer to force constants
	POT_DER_IDENTIFIER *ptr_id; //Pointer to identifier
	POT_DER_LINK_LIST *P2_01;   //Pointer to array of force constants
	void Update(int row, int column, double *k, int *r, double K, double ***D, int dim=-1);


	//ZERO DYNAMICAL MATRIX
	for(d0=0;d0<UC_DOF;d0++) {
		for(d1=0;d1<UC_DOF;d1++) {D[d0][d1][0] = D[d0][d1][1] = 0.0;}
	}


	//EVALUATE DYNAMICAL MATRIX
	for(b0=0;b0<Unit_Cell->natom;b0++) {
		n_i = DIM*b0;
    //F2_01
    ptr_id = FC[b0]->id[2];
    P2_01 = FC[b0]->F2_01;
		for(b1=0;b1<FC[b0]->N[2];b1++) {
			ptr_fc = P2_01->fc;
			n_j = DIM*ptr_id->b[0];
			for(d0=0;d0<DIM;d0++) {
				for(d1=0;d1<DIM;d1++) {
					Update(n_i+d0, n_j+d1, k, ptr_id->l[0], *(ptr_fc++), D, dim);
				}
			}
			ptr_id = ptr_id->next;
			P2_01 = P2_01->next;
		}  //for b1

		//F2_00 (Self term)
		if((dim<0)||(dim>=DIM)) {
		  ptr_fc = FC[b0]->F2_00->fc;
      for(d0=0;d0<DIM;d0++) {
        for(d1=0;d1<DIM;d1++) {
          D[n_i+d0][n_i+d1][0] = D[n_i+d0][n_i+d1][0] + *(ptr_fc++);
        }
      }
		}

	}  //end for b0


	//SET SMALL VALUES TO ZERO
	for(d0=0;d0<UC_DOF;d0++) {
	  if(fabs(D[d0][d0][0])<DM_ZERO) {D[d0][d0][0] = 0.0;}
    if(fabs(D[d0][d0][1])<DM_ZERO) {D[d0][d0][1] = 0.0;}
		for(d1=d0+1;d1<UC_DOF;d1++) {
		  D[d0][d1][0] = 0.5*(D[d0][d1][0]+D[d1][d0][0]);
		  D[d1][d0][0] = D[d0][d1][0];
		  D[d0][d1][1] = 0.5*(D[d0][d1][1]-D[d1][d0][1]);
		  D[d1][d0][1] =-D[d0][d1][1];
			if(fabs(D[d0][d1][0])<DM_ZERO) {D[d0][d1][0] = D[d1][d0][0] = 0.0;}
			if(fabs(D[d0][d1][1])<DM_ZERO) {D[d0][d1][1] = D[d1][d0][1] = 0.0;}
		}
	}

	return;
}


/*SUBROUTINE Update FOR UPDATING THE DYNAMICAL MATRIX OR ITS DERIVATIVE*/
void Update(int row, int column, double *k, int *r, double K, double ***D, int dim=-1) {
	//VARIABLE DECLARATION
	double e[3][2];			//exp(ik*r)
	double a;				    //Dummy variable
	double pi2 = 8.0*atan(1.0);
	//COMPUTE exp[i(k*r)]
	a = pi2*k[0]*r[0];
	e[0][0] = cos(a);
	e[0][1] = sin(a);
	a = pi2*k[1]*r[1];
	e[1][0] = cos(a);
	e[1][1] = sin(a);
	a = pi2*k[2]*r[2];
	e[2][0] = cos(a);
	e[2][1] = sin(a);
	//COMPUTE (d/dk)exp[i(k*r)]
	if((dim>=0)&&(dim<DIM)) {
    a = r[dim]*e[dim][0];
    e[dim][0] = -r[dim]*e[dim][1];
    e[dim][1] = a;
    a = 0.0;
    for(int i=0;i<DIM;i++) {a += Lattice->a[dim][i];}
    a = sqrt(a);
    e[dim][0] *= a;
    e[dim][1] *= a;
	}
  //UPDATE ELEMENTS OF DYNAMICAL MATRIX
	a = e[2][0]*(e[0][0]*e[1][0] - e[0][1]*e[1][1])
	  - e[2][1]*(e[0][0]*e[1][1] + e[1][0]*e[0][1]);
	D[row][column][0] += K*a;
	a = e[2][0]*(e[0][0]*e[1][1] + e[1][0]*e[0][1])
	  + e[2][1]*(e[0][0]*e[1][0] - e[0][1]*e[1][1]);
	D[row][column][1] += K*a;

	return;
}

#undef DM_ZERO
