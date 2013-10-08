/*                           Potential.cpp                          */
/*                            12/02/2008                            */
/*********************************************************************
*    Subroutines for the POTENTIAL class.                            *
*********************************************************************/

/*DECLARE HEADERS*/
#include "LDCode.h"


/*DEFINE SUBROUTINES*/
string Read_Next(ifstream &Input);


/*POTENTIAL CLASS Constructor*/
POTENTIAL::POTENTIAL(int n_mat_, int *mat_) {
  next = NULL;
  n_mat = n_mat_;
  n_expect = ((n_mat+1)*n_mat)/2;
  int i, j;
  mat = new int[n_mat];
  for(i=0;i<n_mat;i++) {
    mat[i] = mat_[i];
    for(j=0;j<i;j++) {
      if(mat[i]==mat[j]) {
        Log <<"When defining potentials the material list is for unique materials."<<endl;
        exit(0);
      }
    }
  }
  return;
}


/*POTENTIAL FUNCTION Fill: ASSIGNS VALUES READ TO PASSED PARAMETER*/
void POTENTIAL::Fill(double *Para, ifstream &Input) {
  for(int i=0;i<n_expect;i++) {Para[i] = atof(Read_Next(Input).c_str());}
  return;
}


/*POTENTIAL::Num_Der COMPUTES n^th ORDER NUMERICAL DERIVATIVE OF POTENTIAL ENERGY*/
#define SMALL 0.008
double POTENTIAL::Num_Der(int n, double **X, double **DX) {
	//EVALUATE THE n^th ORDER DERIVATIVE BY RECURSIVELY CALCULATING FIRST
	  //ORDER DERIVARIVES BY CENTRAL DIFFERENCE OR FIVE-POINT STENCIL
	if(n>=1) {
	  n -= 1;
	  double df_dx;                 //df_dx
	  /*//Central Finite Difference (Order 2^n computations)//
	  DX[n][0] += SMALL;            //x+dx
	  df_dx  = Num_Der(n, X, DX);   //+f(x+dx)
    DX[n][0] -= SMALL+SMALL;      //x-dx
    df_dx -= Num_Der(n, X, DX);   //-f(x-dx)
    DX[n][0] += SMALL;            //x
	  return(df_dx/(SMALL+SMALL));  //df/dx = [f(x+dx)-f(x-dx)]/[2dx]*/
	  //Five-Point Stencil (Order 4^n computations)//
	  DX[n][0] -= SMALL + SMALL;    //x-2dx
	  df_dx  = Num_Der(n, X, DX);   //+f(x-2dx)
	  DX[n][0] += SMALL;            //x-dx
	  df_dx -= 8.*Num_Der(n, X, DX);//-8f(x-dx)
	  DX[n][0] += SMALL + SMALL;    //x+dx
	  df_dx += 8.*Num_Der(n, X, DX);//+8f(x+dx)
	  DX[n][0] += SMALL;            //x+2dx
	  df_dx -= Num_Der(n, X, DX);   //-f(x+2dx)
	  DX[n][0] -= SMALL + SMALL;    //x
	  df_dx = df_dx/(12.0*SMALL);
	  //if( (df_dx<1.0e-10)&&(df_dx>-1.0e-10) ) {return(0.0);}
	  return(df_dx);   //df/dx = [f(x-2dx)-8f(x-dx)+8f(x+dx)-f(x+2dx)]/[12dx]
  }
	else if(n==0) {return(Energy(X));}
	else {Log <<"Error computing numerical derivative."<<endl;exit(0);}
}
#undef SMALL


/*POTENTIAL CLASS Destructor*/
POTENTIAL::~POTENTIAL() {
  delete   next;      next      =NULL;
  delete[] mat;       mat       =NULL;
  delete[] mat_pos;   mat_pos   =NULL;
  delete[] Pot_Symb;  Pot_Symb  =NULL;
  return;
}
