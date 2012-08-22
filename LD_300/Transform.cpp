/*                           Transform.cpp                          */
/*                            03/06/2009                            */
/*********************************************************************
*    Header file and code for the TRANSFORM class object.  TRANSFORM *
*  computes and stores intermediate steps and the final transform	   *
*  in the three and four phonon interaction expressions.  G4[][9]	   *
*  and G3[][9] are similar to each other.  They are to be computed	 *
*  after all wave vectors are specified and are defined as:		    	 *
*  G3[jj][DIM*a'+a"] = \sum_{i,a}\sum_{l',l"}e(kj|ia)				         *
*   *\Phi_{a,a',a"}(0i,l'i',l"i")exp[i(k'.r(l'0)+k".r(l"0))],		     *
*  G4[jj][DIM*a"+a~] = \sum_{i,a}\sum_{i',a'}\sum_{l',l",l~}e(kj|ia) *
*   *e(-kj|i'a')\Phi_{a,a',a",a~}(0i,l'i',l"i",l~i~)exp[-ik.r(l'0)]	 *
*   *exp[ik'.(r(l"0)-r(l~0))],										                   *
*  where jj = i'*neighbors+i".  H3[natom][3] is computed after the	 *
*  second eigenvector is specified (after j_in) and is defined by:	 *
*  H3[i"][a"] = \sum_{i',a'}G4[jj][DIM*a'+a"]e(k'j'|i'a').			     *
*  T3 and T4 are the final step and return the transformed values.	 *
*********************************************************************/


/*DEFINE GLOBAL HEADERS*/
#include "LDCode.h"
#include "Transform.h"
#include <cmath>


/*CONSTRUCTOR: ALLOCATES MEMORY AND FULLY DEFINES pointG[] AND listG[]*/
TRANSFORM::TRANSFORM(PD_ANHARMONIC **FC_) {
	//DECLARE LOCAL VARIABLES
	int j, m;                   //Loops over all entries in each list, counter
	int b0, b1, b2, b3;	        //Atom counters
	int n = Unit_Cell->natom;	  //Number of atoms in unit cell
	POT_DER_IDENTIFIER *id_ptr; //Pointer to identifier
	POT_DER_LINK_LIST *list_ptr;//Pointer to force constant


	//ALLOCATE MEMORY AND INITIALIZE
	FC = FC_;
	bool **temp = new bool* [n];//true if atom combination occurs
	temp[0] = new bool [n*n];
	for(b0=0;b0<n;b0++) {
		if(b0<n-1) {temp[b0+1] = temp[b0] + n;}
		for(b1=0;b1<n;b1++) {temp[b0][b1] = false;}
	}
	pi2 = 8.0*atan(1.0);


	//IDENTIFY AND STORE ATOM PAIRS THAT ARE NEIGHBORS AND APPEAR IN DERIVATIVES
	for(b0=0;b0<n;b0++) {

    id_ptr = FC[b0]->id[2];
    list_ptr = FC[b0]->F3_001;
	  for(j=0;j<FC[b0]->N[2];j++) {
      b1 = id_ptr->b[0];
      //F3_000  //F4_0000, F4_0100
      temp[b0][b0] = true;
      //F3_001  //F4_0001, F4_0101
      temp[b0][b1] = true;
      //F3_010  //F4_0010, F4_0110
      temp[b1][b0] = true;
      //F3_011  //F4_0011, F4_0111
      temp[b1][b1] = true;
      id_ptr = id_ptr->next;
      list_ptr = list_ptr->next;
	  }  //for(j)

    id_ptr = FC[b0]->id[3];
    list_ptr = FC[b0]->F3_012;
    for(j=0;j<FC[b0]->N[3];j++) {
      b1 = id_ptr->b[0];
      b2 = id_ptr->b[1];
      //F3_012  //F4_0012, F4_0112, F4_0212
      temp[b1][b2] = true;
      //F3_021  //F4_0021, F4_0121, F4_0221
      temp[b2][b1] = true;
      //F4_0211
      temp[b1][b1] = true;
      //F4_0122
      temp[b2][b2] = true;
      id_ptr = id_ptr->next;
      list_ptr = list_ptr->next;
    }  //for(j)

    id_ptr = FC[b0]->id[4];
    list_ptr = FC[b0]->F4_0123;
    for(j=0;j<FC[b0]->N[4];j++) {
      b1 = id_ptr->b[0];
      b2 = id_ptr->b[1];
      b3 = id_ptr->b[2];
      //F4_0123
      temp[b2][b3] = true;
      //F4_0312
      temp[b1][b2] = true;
      //F4_0231
      temp[b3][b1] = true;
      //F4_0321
      temp[b2][b1] = true;
      //F4_0132
      temp[b3][b2] = true;
      //F4_0213
      temp[b1][b3] = true;
      id_ptr = id_ptr->next;
      list_ptr = list_ptr->next;
    }  //for(j)

  }  //for(i)


	//COUNT TOTAL NUMBER OF ENTRIES TO GO IN listG
	pointG = new int [n+1];
	pointG[0] = 0;
	for(b0=0;b0<n;b0++) {
		pointG[b0+1] = pointG[b0];
		for(b1=0;b1<n;b1++) {
			if(temp[b0][b1]) {pointG[b0+1] += 1;}
		}
	}

	//BUILD listG AND ALLOCATE MEMORY FOR POINTERS
	listG = new int [pointG[n]];
	G4_re = new double* [pointG[n]];
	G4_re[0] = new double [pointG[n]*DIM*DIM];
	G4_im = new double* [pointG[n]];
	G4_im[0] = new double [pointG[n]*DIM*DIM];
	G3_re = new double* [pointG[n]];
	G3_re[0] = new double [pointG[n]*DIM*DIM];
	G3_im = new double* [pointG[n]];
	G3_im[0] = new double [pointG[n]*DIM*DIM];
	H3_re = new double* [n];
	H3_re[0] = new double [n*DIM];
	H3_im = new double* [n];
	H3_im[0] = new double [n*DIM];
	for(b0=0;b0<n;b0++) {
		if(b0<n-1) {
			H3_re[b0+1] = H3_re[b0] + DIM;
			H3_im[b0+1] = H3_im[b0] + DIM;
		}
		m = 0;
		for(j=pointG[b0];j<pointG[b0+1];j++) {
			if(j!=0) {
				G4_re[j] = G4_re[j-1] + DIM*DIM;
				G4_im[j] = G4_im[j-1] + DIM*DIM;
				G3_re[j] = G3_re[j-1] + DIM*DIM;
				G3_im[j] = G3_im[j-1] + DIM*DIM;
			}
			for(b1=m;b1<n;b1++) {
				if(temp[b0][b1]) {listG[j] = b1;m=b1+1;break;}
			}
		}
	}


	//DEALLOCATE MEMORY AND RETURN
	delete[] temp[0]; temp[0]=NULL;
	delete[] temp;	  temp   =NULL;
	return;
}


/*TRANSFORM DESTRUCTOR: DEALLOCATES MEMORY*/
TRANSFORM::~TRANSFORM(void) {
	//DEALLOCATE DYNAMIC MEMORY
	delete[] pointG;   pointG  =NULL;
	delete[] listG;    listG   =NULL;
	delete[] G3_re[0]; G3_re[0]=NULL;
	delete[] G3_re;    G3_re   =NULL;
	delete[] G3_im[0]; G3_im[0]=NULL;
	delete[] G3_im;    G3_im   =NULL;
	delete[] H3_re[0]; H3_re[0]=NULL;
	delete[] H3_re;    H3_re   =NULL;
	delete[] H3_im[0]; H3_im[0]=NULL;
	delete[] H3_im;    H3_im   =NULL;
	delete[] G4_re[0]; G4_re[0]=NULL;
	delete[] G4_re;    G4_re   =NULL;
	delete[] G4_im[0]; G4_im[0]=NULL;
	delete[] G4_im;    G4_im   =NULL;
	return;
}


/*FUNCTION Transform3G*/
void TRANSFORM::TransformG3(double **E, double **kvect) {
	//DECLARE LOCAL VARIABLES
	int j, m, mm, mmm, i_DOF;		//Counters
	int b0, b1, b2;             //Atom indicies within the unit cell
	int d0, d1, d2;             //Cartesian coordinate counters
  int G00;                    //Location in G3 to write F3_000,
  int G11, G01, G10;          //F3_011, F3_001, F3_010
  int G12, G21;               //F3_012, F3_021
  double ek1r1_re, ek1r1_im;  //exp(i*k1*r1)
  double ek2r1_re, ek2r1_im;  //exp(i*k2*r1)
  double ek1r2_re, ek1r2_im;  //exp(i*k1*r2)
  double ek2r2_re, ek2r2_im;  //exp(i*k2*r2)
  double D11_re, D11_im;      //exp(i*k1*r1)*exp(i*k2*r1)
  double D12_re, D12_im;      //exp(i*k1*r1)*exp(i*k2*r2)
  double D21_re, D21_im;      //exp(i*k1*r2)*exp(i*k2*r1)
  double E010_re, E010_im;    //E[b0*DIM+d0]*exp(i*k1*r1)
  double E001_re, E001_im;    //E[b0*DIM+d0]*exp(i*k2*r1)
  double E011_re, E011_im;    //E[b0*DIM+d0]*exp(i*k1*r1)*exp(i*k2*r1)
  double E012_re, E012_im;    //E[b0*DIM+d0]*exp(i*k1*r1)*exp(i*k2*r2)
  double E021_re, E021_im;    //E[b0*DIM+d0]*exp(i*k1*r2)*exp(i*k2*r1)
	double fc3;				          //Force constant
	POT_DER_IDENTIFIER *id2;
	POT_DER_LINK_LIST *P3_001;
	POT_DER_LINK_LIST *P3_011;
	POT_DER_IDENTIFIER *id3;
	POT_DER_LINK_LIST *P3_012;


	//INITIALIZE G3_re AND G3_im TO ZERO
	for(b0=0;b0<pointG[Unit_Cell->natom];b0++) {
    for(j=0;j<DIM*DIM;j++) {G3_re[b0][j] = G3_im[b0][j] = 0.0;}
	}
	//Initialize integers to eliminate warnings
	G00 = G01 = G10 = G11 = 0;
	G12 = G21 = 0;


	//COMPUTE FOURIER TRANSFORM OF THIRD ORDER DERIVATIVE
	for(b0=0;b0<Unit_Cell->natom;b0++) {
	  i_DOF = b0*DIM;
	  G00 = pointG[b0];
	  while(b0!=listG[G00]) {G00+=1;}
	  m = 0;
	  for(d0=0;d0<DIM;d0++) {
      mmm = i_DOF+d0;
      mm = 0;
      for(d1=0;d1<DIM;d1++) {
        for(d2=0;d2<DIM;d2++) {
          //F3_000
          fc3 = FC[b0]->F3_000->fc[m];
          G3_re[G00][mm] += fc3*E[mmm][0];
          G3_im[G00][mm] += fc3*E[mmm][1];
          m += 1;
          mm += 1;
        }  //for(d2=0)
      }  //for(d1=0)
	  }  //for(d0=0)

    //TWO ATOMS
    id2 = FC[b0]->id[2];
    P3_001 = FC[b0]->F3_001;
    P3_011 = FC[b0]->F3_011;
    for(j=0;j<FC[b0]->N[2];j++) {
      b1 = id2->b[0];
      for(m=pointG[b0];m<pointG[b0+1];m++) {if(b1==listG[m]) {G01=m;break;}}
      for(m=pointG[b1];m<pointG[b1+1];m++) {if(b0==listG[m]) {G10=m;break;}}
      for(m=pointG[b1];m<pointG[b1+1];m++) {if(b1==listG[m]) {G11=m;break;}}
      //k1*r1, k2*r1
      ek1r1_im = ek2r1_im = 0.0;
      for(d0=0;d0<DIM;d0++) {
        ek1r1_im += kvect[1][d0]*id2->l[0][d0];
        ek2r1_im += kvect[2][d0]*id2->l[0][d0];
      }
      ek1r1_im *= pi2;
      ek2r1_im *= pi2;
      //exp(i*k1*r1), exp(i*k2*r1)
      ek1r1_re = cos(ek1r1_im);
      ek1r1_im = sin(ek1r1_im);
      ek2r1_re = cos(ek2r1_im);
      ek2r1_im = sin(ek2r1_im);
      //exp(i*k1*r1)*exp(i*k2*r1)
      D11_re = ek1r1_re*ek2r1_re - ek1r1_im*ek2r1_im;
      D11_im = ek1r1_im*ek2r1_re + ek1r1_re*ek2r1_im;

      m = 0;
      for(d0=0;d0<DIM;d0++) {
        d1 = i_DOF+d0;
        //E[b0*DIM+d0]*exp(i*k1*r1)*exp(i*k2*r1)
        E011_re = E[d1][0]*D11_re - E[d1][1]*D11_im;
        E011_im = E[d1][1]*D11_re + E[d1][0]*D11_im;
        //E[b0*DIM+d0]*exp(i*k2*r1)
        E001_re = E[d1][0]*ek2r1_re - E[d1][1]*ek2r1_im;
        E001_im = E[d1][1]*ek2r1_re + E[d1][0]*ek2r1_im;
        //E[b0*DIM+d0]*exp(i*k1*r1)
        E010_re = E[d1][0]*ek1r1_re - E[d1][1]*ek1r1_im;
        E010_im = E[d1][1]*ek1r1_re + E[d1][0]*ek1r1_im;
        mm = 0;
        for(d1=0;d1<DIM;d1++) {
          for(d2=0;d2<DIM;d2++) {
            mmm = d2*DIM+d1;
            //F3_011
            fc3 = P3_011->fc[m];
            G3_re[G11][mm] += fc3*E011_re;
            G3_im[G11][mm] += fc3*E011_im;
            //F3_001
            fc3 = P3_001->fc[m];
            G3_re[G01][mm] += fc3*E001_re;
            G3_im[G01][mm] += fc3*E001_im;
            //F3_010
            G3_re[G10][mmm] += fc3*E010_re;
            G3_im[G10][mmm] += fc3*E010_im;
            m += 1;
            mm += 1;
          }  //for(d2=0)
        }  //for(d1=0)
      }  //for(d0=0)
      id2 = id2->next;
      P3_001 = P3_001->next;
      P3_011 = P3_011->next;
	  }  //for(j=0)

    //THREE ATOMS
    id3 = FC[b0]->id[3];
    P3_012 = FC[b0]->F3_012;
	  for(j=0;j<FC[b0]->N[3];j++) {
      b1 = id3->b[0];
      b2 = id3->b[1];
      for(m=pointG[b1];m<pointG[b1+1];m++) {if(b2==listG[m]) {G12=m;break;}}
      for(m=pointG[b2];m<pointG[b2+1];m++) {if(b1==listG[m]) {G21=m;break;}}
      //k1*r1, k2*r1, k1*r2, k2*r2
      ek1r1_im = ek2r1_im = 0.0;
      ek1r2_im = ek2r2_im = 0.0;
      for(d0=0;d0<DIM;d0++) {
        ek1r1_im += kvect[1][d0]*id3->l[0][d0];
        ek2r1_im += kvect[2][d0]*id3->l[0][d0];
        ek1r2_im += kvect[1][d0]*id3->l[1][d0];
        ek2r2_im += kvect[2][d0]*id3->l[1][d0];
      }
      ek1r1_im *= pi2;
      ek2r1_im *= pi2;
      ek1r2_im *= pi2;
      ek2r2_im *= pi2;
      //exp(i*k1*r1), exp(i*k2*r1), exp(i*k1*r2), exp(i*k2*r2)
      ek1r1_re = cos(ek1r1_im);
      ek1r1_im = sin(ek1r1_im);
      ek2r1_re = cos(ek2r1_im);
      ek2r1_im = sin(ek2r1_im);
      ek1r2_re = cos(ek1r2_im);
      ek1r2_im = sin(ek1r2_im);
      ek2r2_re = cos(ek2r2_im);
      ek2r2_im = sin(ek2r2_im);
      //exp(i*k1*r1)*exp(i*k2*r2), exp(i*k1*r2)*exp(i*k2*r1)
      D12_re = ek1r1_re*ek2r2_re - ek1r1_im*ek2r2_im;
      D12_im = ek1r1_im*ek2r2_re + ek1r1_re*ek2r2_im;
      D21_re = ek1r2_re*ek2r1_re - ek1r2_im*ek2r1_im;
      D21_im = ek1r2_im*ek2r1_re + ek1r2_re*ek2r1_im;

      m = 0;
      for(d0=0;d0<DIM;d0++) {
        d1 = i_DOF+d0;
        E012_re = E[d1][0]*D12_re - E[d1][1]*D12_im;
        E012_im = E[d1][1]*D12_re + E[d1][0]*D12_im;
        E021_re = E[d1][0]*D21_re - E[d1][1]*D21_im;
        E021_im = E[d1][1]*D21_re + E[d1][0]*D21_im;
        mm = 0;
        for(d1=0;d1<DIM;d1++) {
          for(d2=0;d2<DIM;d2++) {
            mmm = d2*DIM+d1;
            //F3_012
            fc3 = P3_012->fc[m];
            G3_re[G12][mm] += fc3*E012_re;
            G3_im[G12][mm] += fc3*E012_im;
            //F3_021
            G3_re[G21][mmm] += fc3*E021_re;
            G3_im[G21][mmm] += fc3*E021_im;
            m += 1;
            mm += 1;
          }  //for(d2=0)
        }  //for(d1=0)
      }  //for(d0=0)
      id3 = id3->next;
      P3_012 = P3_012->next;
	  }  //for(j=0)

	}  //for(b0=0)

	return;
}


/*FUNCTION ReduceH3*/
void TRANSFORM::TransformH3(double **E) {
	//DECLARE LOCAL VARIABLES
	int b0, j, jj;				//Atom counters
	int d1, d2;			//Cartesian coordinate counters
	int mm;						//Index
	int i_DOF, i_beta;
	double E_re, E_im;


	//INITIALIZE H3_re AND H3_im TO ZERO
	for(b0=0;b0<Unit_Cell->natom;b0++) {
		for(j=0;j<DIM;j++) {H3_re[b0][j] = H3_im[b0][j] = 0.0;}
	}


	//COMPUTE FOURIER TRANSFORM OF THIRD ORDER DERIVATIVE
	for(b0=0;b0<Unit_Cell->natom;b0++) {
	  i_DOF = b0*DIM;

	  for(j=pointG[b0];j<pointG[b0+1];j++) {
      jj = listG[j];
      mm = 0;
      for(d1=0;d1<DIM;d1++) {
        i_beta = i_DOF+d1;
        E_re = E[i_beta][0];
        E_im = E[i_beta][1];
        for(d2=0;d2<DIM;d2++) {
          //F3_000, F3_001, F3_010, F3_011
          H3_re[jj][d2] += G3_re[j][mm]*E_re - G3_im[j][mm]*E_im;
          H3_im[jj][d2] += G3_im[j][mm]*E_re + G3_re[j][mm]*E_im;
          mm += 1;
        }  //for(d2=0)
      }  //for(d1=0)
	  }  //for(j=0)
	}  //for(b0=0)

	return;
}


/*SUBROUTINE ReduceT3*/
double TRANSFORM::TransformT3(double **E) {
	//DECLARE LOCAL VARIABLES
	int b0;						//Atom counters
	int d2;					//Cartesian coordinate counters
	int i_DOF, i_gamma;
	double T_re = 0.0, T_im = 0.0;	//Fourier transform of derivative
	double E_re, E_im;


	//COMPUTE FOURIER TRANSFORM OF THIRD ORDER DERIVATIVE
	for(b0=0;b0<Unit_Cell->natom;b0++) {
	  i_DOF = b0*DIM;
	  for(d2=0;d2<DIM;d2++) {
      i_gamma = i_DOF+d2;
      E_re = E[i_gamma][0];
      E_im = E[i_gamma][1];
      //F3_000, F3_001, F3_010, F3_011, F3_012, F3_021
      T_re += H3_re[b0][d2]*E_re - H3_im[b0][d2]*E_im;
      T_im += H3_im[b0][d2]*E_re + H3_re[b0][d2]*E_im;
	  }  //for(d2=0)
	}  //for(b0=0)

	return(T_re*T_re+T_im*T_im);
}


/*FUNCTION ReduceG4 COMPUTES G4[][] FOR Phi^(4)(kj,-kj,k'j',-k'j')*/
void TRANSFORM::TransformG4(double **E, double **kvect) {
	//DECLARE LOCAL VARIABLES
	int j, m, mm;				        //Counters
	int b0, b1, b2, b3;         //Atom indicies within the unit cell
	int d0, d1, d2, d3;         //Cartesian coordinate counters
  int G00;                    //Location in G4 to write F4_0n00,
  int G11, G01, G10;            //F4_0n11, F4_0n01, F4_0n10
  int G02, G20, G12, G21, G22;  //F4_0n02, F4_0n20, F4_0n12, F4_0n21, F4_0n22
  int G13, G31, G23, G32;       //etc.
  int d0_3;                   //d0*DIM*DIM*DIM
	int d1_2, d2_2, d3_2;       //dn*DIM*DIM + d0_3
	int d1_1, d2_1, d3_1;       //dn*DIM
  int DIM2, DIM3;             //DIM^2, DIM^3
  double x_re, x_im;          //Generic doubles for real and imaginary values
  double ek0r1_re, ek0r1_im;  //exp(i*k0*r1)
  double ek1r1_re, ek1r1_im;  //exp(i*k1*r1)
  double ek0r2_re, ek0r2_im;  //exp(i*k0*r2)
  double ek1r2_re, ek1r2_im;  //exp(i*k1*r2)
  double ek0r3_re, ek0r3_im;  //exp(i*k0*r3)
  double ek1r3_re, ek1r3_im;  //exp(i*k1*r3)
  double D12_re, D12_im;      //exp(i*k1*r1)*exp(-i*k1*r2)
  double D13_re, D13_im;      //exp(i*k1*r1)*exp(-i*k1*r3)
  double D23_re, D23_im;      //exp(i*k1*r2)*exp(-i*k1*r3)
  double E00_re, E00_im;      //E[b0*DIM+d0]*E'[b0*DIM+d1]*exp(-i*k0*r0)
  double E01_re, E01_im;      //E[b0*DIM+d0]*E'[b1*DIM+d1]*exp(-i*k0*r1)
  double E02_re, E02_im;      //E[b0*DIM+d0]*E'[b2*DIM+d1]*exp(-i*k0*r2)
  double E03_re, E03_im;      //E[b0*DIM+d0]*E'[b3*DIM+d1]*exp(-i*k0*r3)
	double fc4;				          //Force constant
	POT_DER_IDENTIFIER *id2;
	POT_DER_LINK_LIST *P4_0001;
	POT_DER_LINK_LIST *P4_0011;
  POT_DER_LINK_LIST *P4_0111;
  POT_DER_IDENTIFIER *id3;
  POT_DER_LINK_LIST *P4_0012;
  POT_DER_LINK_LIST *P4_0112;
  POT_DER_LINK_LIST *P4_0122;
  POT_DER_IDENTIFIER *id4;
  POT_DER_LINK_LIST *P4_0123;


	//INITIALIZE G4_re AND G4_im TO ZERO
  for(d0=0;d0<pointG[Unit_Cell->natom];d0++) {
    for(d1=0;d1<DIM*DIM;d1++) {G4_re[d0][d1] = G4_im[d0][d1] = 0.0;}
	}
	//Initialize integers to eliminate warnings
	G00 = G01 = G10 = G11 = 0;
	G02 = G20 = G12 = G21 = G22 = 0;
	G13 = G31 = G23 = G32 = 0;
	//
	DIM2 = DIM*DIM;
	DIM3 = DIM2*DIM;


	//COMPUTE FOURIER TRANSFORM OF FOURTH DERIVATIVE
	for(b0=0;b0<Unit_Cell->natom;b0++) {

	  G00 = pointG[b0];
	  while(b0!=listG[G00]) {G00+=1;}
	  m = 0;
	  for(d0=0;d0<DIM;d0++) {
      for(d1=0;d1<DIM;d1++) {
        d2 = b0*DIM+d0;
        //E[b0*DIM+d0]*E'[b1*DIM+d1]*exp(-i*k0*r0)
        d3 = b0*DIM+d1;
        E00_re = E[d2][0]*E[d3][0] + E[d2][1]*E[d3][1];
        E00_im = E[d2][1]*E[d3][0] - E[d2][0]*E[d3][1];
        mm = 0;
        for(d2=0;d2<DIM;d2++) {
          for(d3=0;d3<DIM;d3++) {
            //F4_0000
            fc4 = FC[b0]->F4_0000->fc[m];
            G4_re[G00][mm] += fc4*E00_re;
            G4_im[G00][mm] += fc4*E00_im;
            m += 1;
            mm += 1;
          }  //for(d3=0)
        }  //for(d2=0)
      }  //for(d1=0)
	  }  //for(d0=0)

    //TWO ATOMS
    id2 = FC[b0]->id[2];
    P4_0001 = FC[b0]->F4_0001;
    P4_0011 = FC[b0]->F4_0011;
    P4_0111 = FC[b0]->F4_0111;
	  for(j=0;j<FC[b0]->N[2];j++) {
      b1 = id2->b[0];
      for(m=pointG[b0];m<pointG[b0+1];m++) {if(b1==listG[m]) {G01=m;break;}}
      for(m=pointG[b1];m<pointG[b1+1];m++) {if(b0==listG[m]) {G10=m;break;}}
      for(m=pointG[b1];m<pointG[b1+1];m++) {if(b1==listG[m]) {G11=m;break;}}
      //k0*r1, k1*r1
      ek0r1_im = ek1r1_im = 0.0;
      for(d0=0;d0<DIM;d0++) {
        ek0r1_im += kvect[0][d0]*id2->l[0][d0];
        ek1r1_im += kvect[1][d0]*id2->l[0][d0];
      }
      ek0r1_im *= pi2;
      ek1r1_im *= pi2;
      //exp(i*k0*r1), exp(i*k1*r1)
      ek0r1_re = cos(ek0r1_im);
      ek0r1_im = sin(ek0r1_im);
      ek1r1_re = cos(ek1r1_im);
      ek1r1_im = sin(ek1r1_im);

      //m = 0;
      d0_3 = 0;
      for(d0=0;d0<DIM;d0++) {
        d1_2 = d0_3;
        d1_1 = 0;
        for(d1=0;d1<DIM;d1++) {
          d2 = b0*DIM+d0;
          //E[b0*DIM+d0]*E'[b1*DIM+d1]*exp(-i*k0*r0)
          d3 = b0*DIM+d1;
          E00_re = E[d2][0]*E[d3][0] + E[d2][1]*E[d3][1];
          E00_im = E[d2][1]*E[d3][0] - E[d2][0]*E[d3][1];
          //E[b0*DIM+d0]*E'[b1*DIM+d1]
          d3 = b1*DIM+d1;
          x_re = E[d2][0]*E[d3][0] + E[d2][1]*E[d3][1];
          x_im = E[d2][1]*E[d3][0] - E[d2][0]*E[d3][1];
          //E[b0*DIM+d0]*E'[b1*DIM+d1]*exp(-i*k0*r1)
          E01_re = x_re*ek0r1_re + x_im*ek0r1_im;
          E01_im = x_im*ek0r1_re - x_re*ek0r1_im;
          mm = 0;
          d2_2 = d0_3;
          d2_1 = 0;
          for(d2=0;d2<DIM;d2++) {
            d3_2 = d0_3;
            d3_1 = 0;
            for(d3=0;d3<DIM;d3++) {
              //F4_0111   //d1-->d1, d2-->d2, d3-->d3
              fc4 = P4_0111->fc[d1_2+mm];//[((d0*DIM+d1)*DIM+d2)*DIM+d3];
              G4_re[G11][mm] += fc4*(E01_re);
              G4_im[G11][mm] += fc4*(E01_im);

              //F4_0011   //d1-->d1, d2-->d2, d3-->d3
              fc4 = P4_0011->fc[d1_2+mm];//[((d0*DIM+d1)*DIM+d2)*DIM+d3];
              G4_re[G11][mm] += fc4*(E00_re);
              G4_im[G11][mm] += fc4*(E00_im);
              //F4_0101   //d1-->d2, d2-->d1, d3-->d3
              fc4 = P4_0011->fc[d2_2+d1_1+d3];//[((d0*DIM+d2)*DIM+d1)*DIM+d3];
              G4_re[G01][mm] += fc4*(E01_re*ek1r1_re + E01_im*ek1r1_im);
              G4_im[G01][mm] += fc4*(E01_im*ek1r1_re - E01_re*ek1r1_im);
              //F4_0110   //d1-->d3, d2-->d1, d3-->d2
              fc4 = P4_0011->fc[d3_2+d1_1+d2];//[((d0*DIM+d3)*DIM+d1)*DIM+d2];
              G4_re[G10][mm] += fc4*(E01_re*ek1r1_re - E01_im*ek1r1_im);
              G4_im[G10][mm] += fc4*(E01_im*ek1r1_re + E01_re*ek1r1_im);

              //F4_0001   //d1-->d1, d2-->d2, d3-->d3
              fc4 = P4_0001->fc[d1_2+mm];//[((d0*DIM+d1)*DIM+d2)*DIM+d3];
              G4_re[G01][mm] += fc4*(E00_re*ek1r1_re + E00_im*ek1r1_im);
              G4_im[G01][mm] += fc4*(E00_im*ek1r1_re - E00_re*ek1r1_im);
              //F4_0010   //d1-->d1, d2-->d3, d3-->d2
              fc4 = P4_0001->fc[d1_2+d3_1+d2];//[((d0*DIM+d1)*DIM+d3)*DIM+d2];
              G4_re[G10][mm] += fc4*(E00_re*ek1r1_re - E00_im*ek1r1_im);
              G4_im[G10][mm] += fc4*(E00_im*ek1r1_re + E00_re*ek1r1_im);
              //F4_0100   //d1-->d2, d2-->d3, d3-->d1
              fc4 = P4_0001->fc[d2_2+d3_1+d1];//[((d0*DIM+d2)*DIM+d3)*DIM+d1];
              G4_re[G00][mm] += fc4*(E01_re);
              G4_im[G00][mm] += fc4*(E01_im);
              mm += 1;
              d3_2 += DIM2;
              d3_1 += DIM;
            }  //for(d3=0)
            d2_2 += DIM2;
            d2_1 += DIM;
          }  //for(d2=0)
          //m += DIM2;
          d1_2 += DIM2;
          d1_1 += DIM;
        }  //for(d1=0)
        d0_3 += DIM3;
      }  //for(d0=0)
      id2 = id2->next;
      P4_0001 = P4_0001->next;
      P4_0011 = P4_0011->next;
      P4_0111 = P4_0111->next;
	  }  //for(j=0)

	  //THREE ATOMS
	  id3 = FC[b0]->id[3];
	  P4_0012 = FC[b0]->F4_0012;
	  P4_0112 = FC[b0]->F4_0112;
	  P4_0122 = FC[b0]->F4_0122;
    for(j=0;j<FC[b0]->N[3];j++) {
      b1 = id3->b[0];
      b2 = id3->b[1];
      for(m=pointG[b0];m<pointG[b0+1];m++) {if(b1==listG[m]) {G01=m;break;}}
      for(m=pointG[b1];m<pointG[b1+1];m++) {if(b0==listG[m]) {G10=m;break;}}
      for(m=pointG[b1];m<pointG[b1+1];m++) {if(b1==listG[m]) {G11=m;break;}}
      for(m=pointG[b0];m<pointG[b0+1];m++) {if(b2==listG[m]) {G02=m;break;}}
      for(m=pointG[b2];m<pointG[b2+1];m++) {if(b0==listG[m]) {G20=m;break;}}
      for(m=pointG[b1];m<pointG[b1+1];m++) {if(b2==listG[m]) {G12=m;break;}}
      for(m=pointG[b2];m<pointG[b2+1];m++) {if(b1==listG[m]) {G21=m;break;}}
      for(m=pointG[b2];m<pointG[b2+1];m++) {if(b2==listG[m]) {G22=m;break;}}
      //k0*r1, k1*r1, k0*r2, k1*r2
      ek0r1_im = ek1r1_im = 0.0;
      ek0r2_im = ek1r2_im = 0.0;
      for(d0=0;d0<DIM;d0++) {
        ek0r1_im += kvect[0][d0]*id3->l[0][d0];
        ek1r1_im += kvect[1][d0]*id3->l[0][d0];
        ek0r2_im += kvect[0][d0]*id3->l[1][d0];
        ek1r2_im += kvect[1][d0]*id3->l[1][d0];
      }
      ek0r1_im *= pi2;
      ek1r1_im *= pi2;
      ek0r2_im *= pi2;
      ek1r2_im *= pi2;
      //exp(i*k0*r1), exp(i*k1*r1), exp(i*k0*r2), exp(i*k1*r2)
      ek0r1_re = cos(ek0r1_im);
      ek0r1_im = sin(ek0r1_im);
      ek1r1_re = cos(ek1r1_im);
      ek1r1_im = sin(ek1r1_im);
      ek0r2_re = cos(ek0r2_im);
      ek0r2_im = sin(ek0r2_im);
      ek1r2_re = cos(ek1r2_im);
      ek1r2_im = sin(ek1r2_im);
      //exp(i*k1*r1)*exp(-i*k1*r2)
      D12_re = ek1r1_re*ek1r2_re + ek1r1_im*ek1r2_im;
      D12_im = ek1r1_im*ek1r2_re - ek1r1_re*ek1r2_im;

      d0_3 = 0;
      for(d0=0;d0<DIM;d0++) {
        d1_2 = d0_3;
        d1_1 = 0;
        for(d1=0;d1<DIM;d1++) {
          d2 = b0*DIM+d0;
          //E[b0*DIM+d0]*E'[b1*DIM+d1]*exp(-i*k0*r0)
          d3 = b0*DIM+d1;
          E00_re = E[d2][0]*E[d3][0] + E[d2][1]*E[d3][1];
          E00_im = E[d2][1]*E[d3][0] - E[d2][0]*E[d3][1];
          //E[b0*DIM+d0]*E'[b1*DIM+d1]
          d3 = b1*DIM+d1;
          x_re = E[d2][0]*E[d3][0] + E[d2][1]*E[d3][1];
          x_im = E[d2][1]*E[d3][0] - E[d2][0]*E[d3][1];
          //E[b0*DIM+d0]*E'[b1*DIM+d1]*exp(-i*k0*r1)
          E01_re = x_re*ek0r1_re + x_im*ek0r1_im;
          E01_im = x_im*ek0r1_re - x_re*ek0r1_im;
          //E[b0*DIM+d0]*E'[b2*DIM+d1]
          d3 = b2*DIM+d1;
          x_re = E[d2][0]*E[d3][0] + E[d2][1]*E[d3][1];
          x_im = E[d2][1]*E[d3][0] - E[d2][0]*E[d3][1];
          //E[b0*DIM+d0]*E'[b2*DIM+d1]*exp(-i*k0*r2)
          E02_re = x_re*ek0r2_re + x_im*ek0r2_im;
          E02_im = x_im*ek0r2_re - x_re*ek0r2_im;
          mm = 0;
          d2_2 = d0_3;
          d2_1 = 0;
          for(d2=0;d2<DIM;d2++) {
            d3_2 = d0_3;
            d3_1 = 0;
            for(d3=0;d3<DIM;d3++) {
              //F4_0112   //d1-->d1, d2-->d2, d3-->d3
              fc4 = P4_0112->fc[d1_2+mm];//[((d0*DIM+d1)*DIM+d2)*DIM+d3];
              G4_re[G12][mm] += fc4*(E01_re*D12_re - E01_im*D12_im);
              G4_im[G12][mm] += fc4*(E01_im*D12_re + E01_re*D12_im);
              //F4_0121   //d1-->d1, d2-->d3, d3-->d2
              fc4 = P4_0112->fc[d1_2+d3_1+d2];//[((d0*DIM+d1)*DIM+d3)*DIM+d2];
              G4_re[G21][mm] += fc4*(E01_re*D12_re + E01_im*D12_im);
              G4_im[G21][mm] += fc4*(E01_im*D12_re - E01_re*D12_im);
              //F4_0211   //d1-->d2, d2-->d3, d3-->d1
              fc4 = P4_0112->fc[d2_2+d3_1+d1];//[((d0*DIM+d2)*DIM+d3)*DIM+d1];
              G4_re[G11][mm] += fc4*(E02_re);
              G4_im[G11][mm] += fc4*(E02_im);

              //F4_0122   //d1-->d1, d2-->d2, d3-->d3
              fc4 = P4_0122->fc[d1_2+mm];//[((d0*DIM+d1)*DIM+d2)*DIM+d3];
              G4_re[G22][mm] += fc4*(E01_re);
              G4_im[G22][mm] += fc4*(E01_im);
              //F4_0212   //d1-->d2, d2-->d1, d3-->d3
              fc4 = P4_0122->fc[d2_2+d1_1+d3];//[((d0*DIM+d2)*DIM+d1)*DIM+d3];
              G4_re[G12][mm] += fc4*(E02_re*D12_re - E02_im*D12_im);
              G4_im[G12][mm] += fc4*(E02_im*D12_re + E02_re*D12_im);
              //F4_0221   //d1-->d3, d2-->d1, d3-->d2
              fc4 = P4_0122->fc[d3_2+d1_1+d2];//[((d0*DIM+d3)*DIM+d1)*DIM+d2];
              G4_re[G21][mm] += fc4*(E02_re*D12_re + E02_im*D12_im);
              G4_im[G21][mm] += fc4*(E02_im*D12_re - E02_re*D12_im);

              //F4_0012   //d1-->d1, d2-->d2, d3-->d3
              fc4 = P4_0012->fc[d1_2+mm];//[((d0*DIM+d1)*DIM+d2)*DIM+d3];
              G4_re[G12][mm] += fc4*(E00_re*D12_re - E00_im*D12_im);
              G4_im[G12][mm] += fc4*(E00_im*D12_re + E00_re*D12_im);
              //F4_0102   //d1-->d2, d2-->d1, d3-->d3
              fc4 = P4_0012->fc[d2_2+d1_1+d3];//[((d0*DIM+d2)*DIM+d1)*DIM+d3];
              G4_re[G02][mm] += fc4*(E01_re*ek1r2_re + E01_im*ek1r2_im);
              G4_im[G02][mm] += fc4*(E01_im*ek1r2_re - E01_re*ek1r2_im);
              //F4_0120   //d1-->d3, d2-->d1, d3-->d2
              fc4 = P4_0012->fc[d3_2+d1_1+d2];//[((d0*DIM+d3)*DIM+d1)*DIM+d2];
              G4_re[G20][mm] += fc4*(E01_re*ek1r2_re - E01_im*ek1r2_im);
              G4_im[G20][mm] += fc4*(E01_im*ek1r2_re + E01_re*ek1r2_im);
              //F4_0021   //d1-->d1, d2-->d3, d3-->d2
              fc4 = P4_0012->fc[d1_2+d3_1+d2];//[((d0*DIM+d1)*DIM+d3)*DIM+d2];
              G4_re[G21][mm] += fc4*(E00_re*D12_re + E00_im*D12_im);
              G4_im[G21][mm] += fc4*(E00_im*D12_re - E00_re*D12_im);
              //F4_0201   //d1-->d2, d2-->d3, d3-->d1
              fc4 = P4_0012->fc[d2_2+d3_1+d1];//[((d0*DIM+d2)*DIM+d3)*DIM+d1];
              G4_re[G01][mm] += fc4*(E02_re*ek1r1_re + E02_im*ek1r1_im);
              G4_im[G01][mm] += fc4*(E02_im*ek1r1_re - E02_re*ek1r1_im);
              //F4_0210   //d1-->d3, d2-->d2, d3-->d1
              fc4 = P4_0012->fc[d3_2+d2_1+d1];//[((d0*DIM+d3)*DIM+d2)*DIM+d1];
              G4_re[G10][mm] += fc4*(E02_re*ek1r1_re - E02_im*ek1r1_im);
              G4_im[G10][mm] += fc4*(E02_im*ek1r1_re + E02_re*ek1r1_im);
              mm += 1;
              d3_2 += DIM2;
              d3_1 += DIM;
            }  //for(d3=0)
            d2_2 += DIM2;
            d2_1 += DIM;
          }  //for(d2=0)
          d1_2 += DIM2;
          d1_1 += DIM;
        }  //for(d1=0)
        d0_3 += DIM3;
      }  //for(d0=0)
      id3 = id3->next;
      P4_0012 = P4_0012->next;
      P4_0112 = P4_0112->next;
      P4_0122 = P4_0122->next;
    }  //for(j)

    //FOUR ATOMS
    id4 = FC[b0]->id[4];
    P4_0123 = FC[b0]->F4_0123;
    for(j=0;j<FC[b0]->N[4];j++) {
      b1 = id4->b[0];
      b2 = id4->b[1];
      b3 = id4->b[2];
      for(m=pointG[b1];m<pointG[b1+1];m++) {if(b2==listG[m]) {G12=m;break;}}
      for(m=pointG[b2];m<pointG[b2+1];m++) {if(b1==listG[m]) {G21=m;break;}}
      for(m=pointG[b1];m<pointG[b1+1];m++) {if(b3==listG[m]) {G13=m;break;}}
      for(m=pointG[b3];m<pointG[b3+1];m++) {if(b1==listG[m]) {G31=m;break;}}
      for(m=pointG[b3];m<pointG[b3+1];m++) {if(b2==listG[m]) {G32=m;break;}}
      for(m=pointG[b2];m<pointG[b2+1];m++) {if(b3==listG[m]) {G23=m;break;}}
      ek0r1_im = ek1r1_im = 0.0;
      ek0r2_im = ek1r2_im = 0.0;
      ek0r3_im = ek1r3_im = 0.0;
      for(d0=0;d0<DIM;d0++) {
        ek0r1_im += kvect[0][d0]*id4->l[0][d0];
        ek1r1_im += kvect[1][d0]*id4->l[0][d0];
        ek0r2_im += kvect[0][d0]*id4->l[1][d0];
        ek1r2_im += kvect[1][d0]*id4->l[1][d0];
        ek0r3_im += kvect[0][d0]*id4->l[2][d0];
        ek1r3_im += kvect[1][d0]*id4->l[2][d0];
      }
      ek0r1_im *= pi2;
      ek1r1_im *= pi2;
      ek0r2_im *= pi2;
      ek1r2_im *= pi2;
      ek0r3_im *= pi2;
      ek1r3_im *= pi2;
      ek0r1_re = cos(ek0r1_im);
      ek0r1_im = sin(ek0r1_im);
      ek1r1_re = cos(ek1r1_im);
      ek1r1_im = sin(ek1r1_im);
      ek0r2_re = cos(ek0r2_im);
      ek0r2_im = sin(ek0r2_im);
      ek1r2_re = cos(ek1r2_im);
      ek1r2_im = sin(ek1r2_im);
      ek0r3_re = cos(ek0r3_im);
      ek0r3_im = sin(ek0r3_im);
      ek1r3_re = cos(ek1r3_im);
      ek1r3_im = sin(ek1r3_im);
      //exp(i*k1*r1)*exp(-i*k1*r2)
      D12_re = ek1r1_re*ek1r2_re + ek1r1_im*ek1r2_im;
      D12_im = ek1r1_im*ek1r2_re - ek1r1_re*ek1r2_im;
      //exp(i*k1*r1)*exp(-i*k1*r3)
      D13_re = ek1r1_re*ek1r3_re + ek1r1_im*ek1r3_im;
      D13_im = ek1r1_im*ek1r3_re - ek1r1_re*ek1r3_im;
      //exp(i*k1*r2)*exp(-i*k1*r3)
      D23_re = ek1r2_re*ek1r3_re + ek1r2_im*ek1r3_im;
      D23_im = ek1r2_im*ek1r3_re - ek1r2_re*ek1r3_im;

      //m = 0;
      d0_3 = 0;
      for(d0=0;d0<DIM;d0++) {
        d1_2 = d0_3;
        d1_1 = 0;
        for(d1=0;d1<DIM;d1++) {
          d2 = b0*DIM+d0;
          //E[b0*DIM+d0]*E'[b1*DIM+d1]
          d3 = b1*DIM+d1;
          x_re = E[d2][0]*E[d3][0] + E[d2][1]*E[d3][1];
          x_im = E[d2][1]*E[d3][0] - E[d2][0]*E[d3][1];
          //E[b0*DIM+d0]*E'[b1*DIM+d1]*exp(-i*k0*r1)
          E01_re = x_re*ek0r1_re + x_im*ek0r1_im;
          E01_im = x_im*ek0r1_re - x_re*ek0r1_im;
          //E[b0*DIM+d0]*E'[b2*DIM+d1]
          d3 = b2*DIM+d1;
          x_re = E[d2][0]*E[d3][0] + E[d2][1]*E[d3][1];
          x_im = E[d2][1]*E[d3][0] - E[d2][0]*E[d3][1];
          //E[b0*DIM+d0]*E'[b2*DIM+d1]*exp(-i*k0*r2)
          E02_re = x_re*ek0r2_re + x_im*ek0r2_im;
          E02_im = x_im*ek0r2_re - x_re*ek0r2_im;
          //E[b0*DIM+d0]*E'[b3*DIM+d1]
          d3 = b3*DIM+d1;
          x_re = E[d2][0]*E[d3][0] + E[d2][1]*E[d3][1];
          x_im = E[d2][1]*E[d3][0] - E[d2][0]*E[d3][1];
          //E[b0*DIM+d0]*E'[b3*DIM+d1]*exp(-i*k0*r3)
          E03_re = x_re*ek0r3_re + x_im*ek0r3_im;
          E03_im = x_im*ek0r3_re - x_re*ek0r3_im;
          mm = 0;
          d2_2 = d0_3;
          d2_1 = 0;
          for(d2=0;d2<DIM;d2++) {
            d3_2 = d0_3;
            d3_1 = 0;
            for(d3=0;d3<DIM;d3++) {
              //F4_0123   //d1-->d1, d2-->d2, d3-->d3
              fc4 = P4_0123->fc[d1_2+mm];//[((d0*DIM+d1)*DIM+d2)*DIM+d3];
              G4_re[G23][mm] += fc4*(E01_re*D23_re - E01_im*D23_im);
              G4_im[G23][mm] += fc4*(E01_im*D23_re + E01_re*D23_im);
              //F4_0132   //d1-->d1, d2-->d3, d3-->d2
              fc4 = P4_0123->fc[d1_2+d3_1+d2];//[((d0*DIM+d1)*DIM+d3)*DIM+d2];
              G4_re[G32][mm] += fc4*(E01_re*D23_re + E01_im*D23_im);
              G4_im[G32][mm] += fc4*(E01_im*D23_re - E01_re*D23_im);
              //F4_0213   //d1-->d2, d2-->d1, d3-->d3
              fc4 = P4_0123->fc[d2_2+d1_1+d3];//[((d0*DIM+d2)*DIM+d1)*DIM+d3];
              G4_re[G13][mm] += fc4*(E02_re*D13_re - E02_im*D13_im);
              G4_im[G13][mm] += fc4*(E02_im*D13_re + E02_re*D13_im);
              //F4_0231   //d1-->d3, d2-->d1, d3-->d2
              fc4 = P4_0123->fc[d3_2+d1_1+d2];//[((d0*DIM+d3)*DIM+d1)*DIM+d2];
              G4_re[G31][mm] += fc4*(E02_re*D13_re + E02_im*D13_im);
              G4_im[G31][mm] += fc4*(E02_im*D13_re - E02_re*D13_im);
              //F4_0312   //d1-->d2, d2-->d3, d3-->d1
              fc4 = P4_0123->fc[d2_2+d3_1+d1];//[((d0*DIM+d2)*DIM+d3)*DIM+d1];
              G4_re[G12][mm] += fc4*(E03_re*D12_re - E03_im*D12_im);
              G4_im[G12][mm] += fc4*(E03_im*D12_re + E03_re*D12_im);
              //F4_0321   //d1-->d3, d2-->d2, d3-->d1
              fc4 = P4_0123->fc[d3_2+d2_1+d1];//[((d0*DIM+d3)*DIM+d2)*DIM+d1];
              G4_re[G21][mm] += fc4*(E03_re*D12_re + E03_im*D12_im);
              G4_im[G21][mm] += fc4*(E03_im*D12_re - E03_re*D12_im);
              //m += 1;
              mm += 1;
              d3_2 += DIM2;
              d3_1 += DIM;
            }  //for(d3=0)
            d2_2 += DIM2;
            d2_1 += DIM;
          }  //for(d2=0)
          d1_2 += DIM2;
          d1_1 += DIM;
        }  //for(d1=0)
        d0_3 += DIM3;
      }  //for(d0=0)
      id4 = id4->next;
      P4_0123 = P4_0123->next;
    }  //for(j)

	}  //for(b0=0)

	return;
}


/*SUBROUTINE TransformT4*/
double TRANSFORM::TransformT4(double **E) {
	//DECLARE LOCAL VARIABLES
	int i, j, jj;				        //Atom counters
	int d2, d3;			            //Cartesian coordinate counters
	int mm;						          //Index
	double T_re = 0.0;
	//double T_im = 0.0;


	//COMPUTE FOURIER TRANSFORM OF FOURTH ORDER DERIVATIVE
	for(i=0;i<Unit_Cell->natom;i++) {
	  for(j=pointG[i];j<pointG[i+1];j++) {
      jj = listG[j];
      mm = 0;
      for(d2=0;d2<DIM;d2++) {
        for(d3=0;d3<DIM;d3++) {
          //All fourth order derivatives
          T_re += G4_re[j][mm]*(E[i*DIM+d2][0]*E[jj*DIM+d3][0]+E[i*DIM+d2][1]*E[jj*DIM+d3][1])
                - G4_im[j][mm]*(E[i*DIM+d2][1]*E[jj*DIM+d3][0]-E[i*DIM+d2][0]*E[jj*DIM+d3][1]);
          /*T_im += G4_im[j][mm]*(E[i*DIM+d2][0]*E[jj*DIM+d3][0]+E[i*DIM+d2][1]*E[jj*DIM+d3][1])
                + G4_re[j][mm]*(E[i*DIM+d2][1]*E[jj*DIM+d3][0]-E[i*DIM+d2][0]*E[jj*DIM+d3][1]);*/
          mm += 1;
        }  //for(d3=0)
      }  //for(d2=0)
    }  //for(j=0)
	}  //for(i=0)

	return(T_re);
}


/*SUBROUTINE Transform_3_0 COMPUTES Phi^(3)(kj,-kj,k''j'')*Phi^(3)(k'j',-k'j',-k''j'') where k''j'' must be a zone center optical mode*/
double TRANSFORM::Transform_3_0(double ***E, double **kvect) {
	//DECLARE LOCAL VARIABLES
	int b0, j, jj;				//Atom counters
	int d0, d1, d2;		//Cartesian coordinate counters
	int m, mm, mmm;				//Indices
	double fc3;				//Constant*ALD_Neighbor->d3_b2_l2[j][m]
	double dot_j0_re, dot_j0_im;//Dot product of atom j and kvect
	double dot_j1_re, dot_j1_im;//Dot product of atom j and kvect'
	double T_re = 0.0, T_im = 0.0;	//Fourier transform of derivative
	double U_re = 0.0, U_im = 0.0;	//Fourier transform of derivative
	double K_iii_re, K_iii_im;	//Intermediate values
	double K_iij_re, K_iij_im;	//Intermediate values
	double K_iji_re, K_iji_im;	//Intermediate values
	double K_ijj_re, K_ijj_im;	//Intermediate values
	double L_iii_re, L_iii_im;	//Intermediate values
	double L_iij_re, L_iij_im;	//Intermediate values
	double L_iji_re, L_iji_im;	//Intermediate values
	double L_ijj_re, L_ijj_im;	//Intermediate values
	POT_DER_IDENTIFIER *id2;
	POT_DER_LINK_LIST *P3_001;
	POT_DER_LINK_LIST *P3_011;
	POT_DER_IDENTIFIER *id3;
	POT_DER_LINK_LIST *P3_012;


	//COMPUTE FOURIER TRANSFORM OF THIRD ORDER DERIVATIVE
	for(b0=0;b0<Unit_Cell->natom;b0++) {
	  m = 0;
	  mm = 0;
	  for(d0=0;d0<DIM;d0++) {
      mm = b0*DIM+d0;
      for(d1=0;d1<DIM;d1++) {
        mmm = b0*DIM+d1;
        K_iii_re = E[0][mm][0]*E[0][mmm][0]+E[0][mm][1]*E[0][mmm][1];
        K_iii_im = E[0][mm][1]*E[0][mmm][0]-E[0][mm][0]*E[0][mmm][1];
        L_iii_re = E[1][mm][0]*E[1][mmm][0]+E[1][mm][1]*E[1][mmm][1];
        L_iii_im = E[1][mm][1]*E[1][mmm][0]-E[1][mm][0]*E[1][mmm][1];
        for(d2=0;d2<DIM;d2++) {
          mmm = b0*DIM+d2;
          //F3_000
          fc3 = FC[b0]->F3_000->fc[m];
          /*Phi^(3)(kj,-kj,k''j'')*/
          T_re += fc3*( K_iii_re*E[2][mmm][0] - K_iii_im*E[2][mmm][1] );
          T_im += fc3*( K_iii_im*E[2][mmm][0] + K_iii_re*E[2][mmm][1] );

          /*Phi^(3)(k'j',-k'j',-k''j'')*/
          U_re += fc3*( L_iii_re*E[2][mmm][0] + L_iii_im*E[2][mmm][1] );
          U_im += fc3*( L_iii_im*E[2][mmm][0] - L_iii_re*E[2][mmm][1] );
          m += 1;
        }  //for(d2=0)
      }  //for(d1=0)
	  }  //for(d0=0)

    id2 = FC[b0]->id[2];
    P3_001 = FC[b0]->F3_001;
    P3_011 = FC[b0]->F3_011;
	  for(j=0;j<FC[b0]->N[2];j++) {
      jj = id2->b[0];
      dot_j0_im = dot_j1_im = 0.0;
      for(d0=0;d0<DIM;d0++) {
        dot_j0_im += kvect[0][d0]*id2->l[0][d0];
        dot_j1_im += kvect[1][d0]*id2->l[0][d0];
      }
      dot_j0_im *= pi2;
      dot_j1_im *= pi2;
      dot_j0_re = cos(dot_j0_im);
      dot_j0_im = sin(dot_j0_im);
      dot_j1_re = cos(dot_j1_im);
      dot_j1_im = sin(dot_j1_im);

      m = 0;
      for(d0=0;d0<DIM;d0++) {
        mm = b0*DIM+d0;
        K_iii_re = E[0][mm][0]*dot_j0_re + E[0][mm][1]*dot_j0_im;
        K_iii_im = E[0][mm][1]*dot_j0_re - E[0][mm][0]*dot_j0_im;
        L_iii_re = E[1][mm][0]*dot_j1_re + E[1][mm][1]*dot_j1_im;
        L_iii_im = E[1][mm][1]*dot_j1_re - E[1][mm][0]*dot_j1_im;
        for(d1=0;d1<DIM;d1++) {
          mmm = b0*DIM+d1;
          K_iij_re = E[0][mm][0]*E[0][mmm][0]+E[0][mm][1]*E[0][mmm][1];
          K_iij_im = E[0][mm][1]*E[0][mmm][0]-E[0][mm][0]*E[0][mmm][1];
          L_iij_re = E[1][mm][0]*E[1][mmm][0]+E[1][mm][1]*E[1][mmm][1];
          L_iij_im = E[1][mm][1]*E[1][mmm][0]-E[1][mm][0]*E[1][mmm][1];

          K_iji_re = K_iii_re*E[2][mmm][0] - K_iii_im*E[2][mmm][1];
          K_iji_im = K_iii_im*E[2][mmm][0] + K_iii_re*E[2][mmm][1];
          L_iji_re = L_iii_re*E[2][mmm][0] + L_iii_im*E[2][mmm][1];
          L_iji_im = L_iii_im*E[2][mmm][0] - L_iii_re*E[2][mmm][1];

          mmm = jj*DIM+d1;
          K_ijj_re = K_iii_re*E[0][mmm][0] + K_iii_im*E[0][mmm][1];
          K_ijj_im = K_iii_im*E[0][mmm][0] - K_iii_re*E[0][mmm][1];
          L_ijj_re = L_iii_re*E[1][mmm][0] + L_iii_im*E[1][mmm][1];
          L_ijj_im = L_iii_im*E[1][mmm][0] - L_iii_re*E[1][mmm][1];
          for(d2=0;d2<DIM;d2++) {
            mmm = jj*DIM+d2;
            /*Phi^(3)(kj,-kj,k''j'')*/
            //F3_011
            fc3 = P3_011->fc[m];
            T_re += fc3*(K_ijj_re*E[2][mmm][0] - K_ijj_im*E[2][mmm][1]);
            T_im += fc3*(K_ijj_im*E[2][mmm][0] + K_ijj_re*E[2][mmm][1]);
            //F3_001
            fc3 = P3_001->fc[m];
            T_re += fc3*(K_iij_re*E[2][mmm][0] - K_iij_im*E[2][mmm][1]);
            T_im += fc3*(K_iij_im*E[2][mmm][0] + K_iij_re*E[2][mmm][1]);
            //F3_010
            T_re += fc3*(K_iji_re*E[0][mmm][0] + K_iji_im*E[0][mmm][1]);
            T_im += fc3*(K_iji_im*E[0][mmm][0] - K_iji_re*E[0][mmm][1]);

            /*Phi^(3)(k'j',-k'j',-k''j'')*/
            //F3_011
            fc3 = P3_011->fc[m];
            U_re += fc3*(L_ijj_re*E[2][mmm][0] + L_ijj_im*E[2][mmm][1]);
            U_im += fc3*(L_ijj_im*E[2][mmm][0] - L_ijj_re*E[2][mmm][1]);
            //F3_001
            fc3 = P3_001->fc[m];
            U_re += fc3*(L_iij_re*E[2][mmm][0] + L_iij_im*E[2][mmm][1]);
            U_im += fc3*(L_iij_im*E[2][mmm][0] - L_iij_re*E[2][mmm][1]);
            //F3_010
            U_re += fc3*(L_iji_re*E[1][mmm][0] + L_iji_im*E[1][mmm][1]);
            U_im += fc3*(L_iji_im*E[1][mmm][0] - L_iji_re*E[1][mmm][1]);
            m += 1;
          }  //for(d2=0)
        }  //for(d1=0)
      }  //for(d0=0)
      id2 = id2->next;
      P3_001 = P3_001->next;
      P3_011 = P3_011->next;
    }  //for(j=0)


    int kk;
    double dot_k0_re, dot_k0_im;
    double dot_k1_re, dot_k1_im;
    double K_ijk_re, K_ijk_im;		//Intermediate values
    double K_ikj_re, K_ikj_im;		//Intermediate values
    double L_ijk_re, L_ijk_im;		//Intermediate values
    double L_ikj_re, L_ikj_im;		//Intermediate values
    id3 = FC[b0]->id[3];
    P3_012 = FC[b0]->F3_012;
    for(j=0;j<FC[b0]->N[3];j++) {
      jj = id3->b[0];
      kk = id3->b[1];
      dot_j0_im = dot_j1_im = 0.0;
      dot_k0_im = dot_k1_im = 0.0;
      for(d0=0;d0<DIM;d0++) {
        dot_j0_im += kvect[0][d0]*id3->l[0][d0];
        dot_j1_im += kvect[1][d0]*id3->l[0][d0];
        dot_k0_im += kvect[0][d0]*id3->l[1][d0];
        dot_k1_im += kvect[1][d0]*id3->l[1][d0];
      }
      dot_j0_im *= pi2;
      dot_j1_im *= pi2;
      dot_k0_im *= pi2;
      dot_k1_im *= pi2;
      dot_j0_re = cos(dot_j0_im);
      dot_j0_im = sin(dot_j0_im);
      dot_j1_re = cos(dot_j1_im);
      dot_j1_im = sin(dot_j1_im);
      dot_k0_re = cos(dot_k0_im);
      dot_k0_im = sin(dot_k0_im);
      dot_k1_re = cos(dot_k1_im);
      dot_k1_im = sin(dot_k1_im);

      m = 0;
      for(d0=0;d0<DIM;d0++) {
        mm = b0*DIM+d0;
        K_iii_re = E[0][mm][0]*dot_j0_re + E[0][mm][1]*dot_j0_im;
        K_iii_im = E[0][mm][1]*dot_j0_re - E[0][mm][0]*dot_j0_im;
        L_iii_re = E[1][mm][0]*dot_j1_re + E[1][mm][1]*dot_j1_im;
        L_iii_im = E[1][mm][1]*dot_j1_re - E[1][mm][0]*dot_j1_im;

        K_iji_re = E[0][mm][0]*dot_k0_re + E[0][mm][1]*dot_k0_im;
        K_iji_im = E[0][mm][1]*dot_k0_re - E[0][mm][0]*dot_k0_im;
        L_iji_re = E[1][mm][0]*dot_k1_re + E[1][mm][1]*dot_k1_im;
        L_iji_im = E[1][mm][1]*dot_k1_re - E[1][mm][0]*dot_k1_im;
        for(d1=0;d1<DIM;d1++) {
          mmm = jj*DIM+d1;
          K_ijk_re = K_iii_re*E[0][mmm][0] + K_iii_im*E[0][mmm][1];
          K_ijk_im = K_iii_im*E[0][mmm][0] - K_iii_re*E[0][mmm][1];
          L_ijk_re = L_iii_re*E[1][mmm][0] + L_iii_im*E[1][mmm][1];
          L_ijk_im = L_iii_im*E[1][mmm][0] - L_iii_re*E[1][mmm][1];

          K_ikj_re = K_iji_re*E[2][mmm][0] - K_iji_im*E[2][mmm][1];
          K_ikj_im = K_iji_im*E[2][mmm][0] + K_iji_re*E[2][mmm][1];
          L_ikj_re = L_iji_re*E[2][mmm][0] + L_iji_im*E[2][mmm][1];
          L_ikj_im = L_iji_im*E[2][mmm][0] - L_iji_re*E[2][mmm][1];

          for(d2=0;d2<DIM;d2++) {
            mmm = kk*DIM+d2;
            //F3_012
            fc3 = P3_012->fc[m];
            T_re += fc3*(K_ijk_re*E[2][mmm][0] - K_ijk_im*E[2][mmm][1]);
            T_im += fc3*(K_ijk_im*E[2][mmm][0] + K_ijk_re*E[2][mmm][1]);
            //F3_021
            T_re += fc3*(K_ikj_re*E[0][mmm][0] + K_ikj_im*E[0][mmm][1]);
            T_im += fc3*(K_ikj_im*E[0][mmm][0] - K_ikj_re*E[0][mmm][1]);


            //F3_012
            U_re += fc3*(L_ijk_re*E[2][mmm][0] + L_ijk_im*E[2][mmm][1]);
            U_im += fc3*(L_ijk_im*E[2][mmm][0] - L_ijk_re*E[2][mmm][1]);
            //F3_021
            U_re += fc3*(L_ikj_re*E[1][mmm][0] + L_ikj_im*E[1][mmm][1]);
            U_im += fc3*(L_ikj_im*E[1][mmm][0] - L_ikj_re*E[1][mmm][1]);
            m += 1;
          }  //for(d2=0)
        }  //for(d1=0)
      }  //for(d0=0)
      id3 = id3->next;
      P3_012 = P3_012->next;
    }  //for(j)

	}  //for(b0=0)

	return(T_re*U_re-T_im*U_im);
}


/*NEED ABOUT 10 INTERACTING ATOMS IN THE UNIT CELL FOR THIS TO BE WORTHWHILE*/
int TRANSFORM::Location(int b0, int b1) {
  //DECLARE LOCAL VARIABLES
  int bot, n, top;
  int test;

  //LOCATE POSITION BY BISECTION METHOD
  bot = pointG[b0];
  top = pointG[b0+1];
  while(bot!=top) {
    n = (bot + top) >> 1;
    test = listG[n];
    if(b1==test) {return(n);}
    else if(b1<test) {top=n;}
    else {bot=n;}
  }
  return(bot);
}
