/*						              POT_TERSOFF.cpp              						*/
/*							              04/12/2009	              						*/
/*********************************************************************
*    File that contains the POT_TERSOFF class function definitions.	 *
*  POT_TERSOFF is a derivative of the POTENTIAL class.			       	 *
*********************************************************************/

/*DEFINE HEADERS*/
#include "Tersoff.h"
#include <cmath>


/*DEFINE SUBROUTINES*/
string Read_Next(ifstream &Input);


/*POT_TERSOFF CONSTRUCTOR: STORES POTENTIAL AND POTENTIAL PARAMETER NAMES*/
POT_TERSOFF::POT_TERSOFF(int n_mat_, int *mat_) : POTENTIAL(n_mat_, mat_) {
	Pot_Symb = new char[8];
  strcpy(Pot_Symb, "TERSOFF");
  if(n_mat!=1) {Log <<"Currently only one unique material allowed for "<<Pot_Symb<<endl;exit(0);}
  Log <<"Warning.  Potential "<<Pot_Symb<<" is untested."<<endl;
  mat_pos = new int[n_mat];
  A      = new double[n_expect];
  B      = new double[n_expect];
  lambda = new double[n_expect];
  mu     = new double[n_expect];
  beta   = new double[n_expect];
  n      = new double[n_expect];
  c      = new double[n_expect];
  d      = new double[n_expect];
  h      = new double[n_expect];
  R      = new double[n_expect];
  cutoff = new double[n_expect];
  for(int i=0;i<n_expect;i++) {
    A[i]=B[i]=lambda[i]=mu[i]=beta[i]=n[i]=c[i]=d[i]=h[i]=R[i]=0.0;
  }
  return;
}


/*POT_TERSOFF ReadInput: READS AND STORES POTENTIAL PARAMETER VALUES*/
string POT_TERSOFF::ReadInput(ifstream &Input) {
  string str;
  while(!Input.eof()) {
    str = Read_Next(Input);
    if(str.compare(0, 3, "cutoff", 3)==0) {Fill(cutoff, Input);} else
    if(str.compare(0, 1, "A"     , 1)==0) {Fill(A     , Input);} else
    if(str.compare(0, 1, "B"     , 1)==0) {Fill(B     , Input);} else
    if(str.compare(0, 3, "lambda", 3)==0) {Fill(lambda, Input);} else
    if(str.compare(0, 2, "mu"    , 2)==0) {Fill(mu    , Input);} else
    if(str.compare(0, 3, "beta"  , 3)==0) {Fill(beta  , Input);} else
    if(str.compare(0, 1, "n"     , 1)==0) {Fill(n     , Input);} else
    if(str.compare(0, 1, "c"     , 1)==0) {Fill(c     , Input);} else
    if(str.compare(0, 1, "d"     , 1)==0) {Fill(d     , Input);} else
    if(str.compare(0, 1, "h"     , 1)==0) {Fill(h     , Input);} else
    if(str.compare(0, 1, "R"     , 1)==0) {Fill(R     , Input);} else
    {break;}
  }
  return(str);
}


/*POT_TERSOFF FUNCTION GetScale: GETS THE LENGTH AND ENERGY SCALES*/
void POT_TERSOFF::GetScale(double &e, double &l) {
  e = 0.0;
  l = 0.0;
	for(int i=0;i<n_expect;i++) {
	  e += A[i] + B[i];
	  l += R[i];
	}
	e /= 2.0*n_expect;
	l /= n_expect;
	return;
}


/*POT_TERSOFF FUNCTION Initialize: INITIALIZES THE POTENTIAL PARAMETERS*/
void POT_TERSOFF::Initialize() {
	for(int i=0;i<n_expect;i++) {
	  A[i]      /= Parameter->energy;
	  B[i]      /= Parameter->energy;
	  cutoff[i] /= Parameter->length;
	  R[i]      /= Parameter->length;
	  lambda[i] *= Parameter->length;
	  mu[i]     *= Parameter->length;
	}
	return;
}


/*POT_TERSOFF FUNCTION Print: PRINTS POTENTIAL PARAMETERS*/
void POT_TERSOFF::Print() {
  int i, j;
  Log <<"\nPOTENTIAL "<<Pot_Symb<<" "<<n_mat;
  for(i=0;i<n_mat;i++) {Log <<" "<<Material->symb[mat[i]];}
  Log <<"  \n  cutoff =";
  for(i=0;i<n_expect;i++) {Log <<" "<<cutoff[i];}
  Log <<" m\n  A      =";
  for(i=0;i<n_expect;i++) {Log <<" "<<A[i];}
  Log <<" J\n  B      =";
  for(i=0;i<n_expect;i++) {Log <<" "<<B[i];}
  Log <<" J\n  lambda =";
  for(i=0;i<n_expect;i++) {Log <<" "<<lambda[i];}
  Log <<" 1/m\n  mu     =";
  for(i=0;i<n_expect;i++) {Log <<" "<<mu[i];}
  Log <<" 1/m\n  beta   =";
  for(i=0;i<n_expect;i++) {Log <<" "<<beta[i];}
  Log <<"  \n  n      =";
  for(i=0;i<n_expect;i++) {Log <<" "<<n[i];}
  Log <<"  \n  c      =";
  for(i=0;i<n_expect;i++) {Log <<" "<<c[i];}
  Log <<"  \n  d      =";
  for(i=0;i<n_expect;i++) {Log <<" "<<d[i];}
  Log <<"  \n  h      =";
  for(i=0;i<n_expect;i++) {Log <<" "<<h[i];}
  Log <<"  \n  R      =";
  for(i=0;i<n_expect;i++) {Log <<" "<<R[i];}
  Log <<" m"<<endl;
  return;
}


/*POT_TERSOFF FUNCTION Copy: RETURNS A POINTER TO A NEW POTENTIAL WITH A DUPLICATE OF THE DATA*/
POTENTIAL *POT_TERSOFF::Copy() {
  int i, j;
  POT_TERSOFF *P_new = new POT_TERSOFF(n_mat, mat);
  for(i=0;i<n_expect;i++) {
    P_new->cutoff[i] = cutoff[i];
    P_new->A[i]      = A[i];
    P_new->B[i]      = B[i];
    P_new->lambda[i] = lambda[i];
    P_new->mu[i]     = mu[i];
    P_new->beta[i]   = beta[i];
    P_new->n[i]      = n[i];
    P_new->c[i]      = c[i];
    P_new->d[i]      = d[i];
    P_new->h[i]      = h[i];
    P_new->R[i]      = R[i];
  }
  return(P_new);
}


/*
  POT_TERSOFF FUNCTION Energy: COMPUTES THE POTENTIAL ENERGY
Pass NULL as the last augument to indicate that the list of atoms
has ended.
*/
double POT_TERSOFF::Energy(double **X) {
  //DECLARE LOCAL VARIABLES
  int i, j, k;
  int dim, m;
  double b, g, zeta, f_C;
	double f_R, f_A, phi;
	double *r_ij = new double[DIM];
	double *r_ik = new double[DIM];
	double r_ij_mag, r_ik_mag;
	double theta_ijk;


	//DETERMINE WHICH CONSTANTS ARE NEEDED
	i = mat_pos[0];
	j = k = mat_pos[1];
	dim = (j>k) ? j+((k+1)*k)/2 : k+((j+1)*j)/2;
	j = (i>j) ? i+((j+1)*j)/2 : j+((i+1)*i)/2;
	k = (i>k) ? i+((k+1)*k)/2 : k+((i+1)*i)/2;


	//COMPUTE ENERGY
	r_ij_mag = 0.0;
	for(dim=0;dim<DIM;dim++) {
	  r_ij[dim] = X[0][dim] - X[1][dim];
	  r_ij_mag += r_ij[dim]*r_ij[dim];
	}
	r_ij_mag = sqrt(r_ij_mag);
	if(r_ij_mag>cutoff[j]) {return(0.0);}
	m=2;
	zeta = 0.0;
	while(X[m]!=NULL) {
	  r_ik_mag = theta_ijk = 0.0;
	  for(dim=0;dim<DIM;dim++) {
      r_ik[dim] = X[0][dim] - X[m][dim];
      r_ik_mag += r_ik[dim]*r_ik[dim];
      theta_ijk += r_ij[dim]*r_ik[dim];
    }
    m += 1;
    r_ik_mag = sqrt(r_ik_mag);
    if(r_ik_mag>cutoff[k]) {continue;}
    theta_ijk /= r_ij_mag*r_ik_mag;
    if(r_ik_mag<R[k]) {f_C = 1.0;}
    else {f_C = 0.5 + 0.5*cos(3.141592654*(r_ij_mag-R[j])/(cutoff[j]-R[j]));}
    b = h[i] - theta_ijk;
    g = 1.0 + c[j]*c[j]/(d[j]*d[j]) - c[j]*c[j]/(d[j]*d[j] + b*b);
    zeta += f_C*g;
	}
	g = pow(beta[j]*zeta, n[j]);
	b = pow(1.0+g, -0.5/n[j]);
	if(r_ij_mag<R[j]) {f_C = 1.0;}
	else {f_C = 0.5 + 0.5*cos(3.141592654*(r_ij_mag-R[j])/(cutoff[j]-R[j]));}
	f_R = A[j]*exp(-lambda[j]*r_ij_mag);
	f_A = -B[j]*exp(-mu[j]*r_ij_mag);
	phi = f_C*(f_R + b*f_A);


  //DEALLOCATE MEMORY AND RETURN ENERGY
  delete[] r_ij; r_ij=NULL;
  delete[] r_ik; r_ik=NULL;
	return(phi);
}


/*POT_TERSOFF::ForceConsants: COMPUTES THE HARMONIC FORCE CONSTANTS*/
void POT_TERSOFF::ForceConstants(int b0, POT_DER_X *F) {
  //DECLARE LOCAL VARIABLES
  int n0, n1, n2;             //Material's position in array mat[]
  int b1, b2, i;              //Counters
  int n_neig1, n_neig2;       //Number of neighbors in lists
  int *B, **L;                //Stores atom and unit cell numbers
  int *l0;                    //Array of zeros (center unit cell)
  double cutoff_;             //Stores proper cutoff
  double **X;                 //Array of atomic positions used in potential
  N_LIST *N_list1;            //Linked list of neighbors to b0
  N_LIST *N_list1_ptr;        //Pointer to neighbors of b0
  N_LIST *N_list2;            //Linked list of neighbors to b1
  N_LIST *N_list2_ptr;        //Pointer to neighbors of b1
  N_LIST *BuildNeighborList_R(int &, int, double);
  void ManageNeighborList(int&, N_LIST*, int, int*, int, int*);


  //INITIALIZE VARIABLES
  X = new double*[6];
  X[0] = new double[6*DIM];
  X[1] = X[0] + DIM;
  X[2] = X[1] + DIM;
  X[3] = X[2] + DIM;
  X[4] = X[3] + DIM;
  X[5] = X[4] + DIM;
  B = new int[6];
  L = new int*[6];
  l0 = new int[DIM];
  for(i=0;i<DIM;i++) {l0[i] = 0;}
  for(n0=0;n0<n_mat;n0++) {if(mat[n0]==Unit_Cell->mat[b0]) {break;}}


  //LOOP OVER ALL MATERIAL INTERACTIONS
  for(n1=0;n1<n_mat;n1++) {


/** b0 at the center **/

    //BUILD NEIGHBOR LISTS
    cutoff_ = cutoff[(n0<=n1) ? n0*((n1+1)*n1)/2 : n1*((n0+1)*n0)/2];
    n_neig1 = b0;
    N_list1 = BuildNeighborList_R(n_neig1, mat[n1], cutoff_);
    n_neig1 += 2;
    mat_pos = new int[n_neig1];
    X = new double*[n_neig1];
    X[0] = new double[n_neig1*DIM];
    B = new int[n_neig1];
    L = new int*[n_neig1];
    n_neig1 -= 2;
    for(b1=0;b1<n_neig1;b1++) {
      mat_pos[b1+1] = n1;
      X[b1+1] = X[b1] + DIM;
    }
    mat_pos[0] = n0;
    mat_pos[n_neig1+1] = -1;
    X[n_neig1+1] = NULL;
    L[n_neig1+1] = NULL;


    //COMPUTE FORCE CONSTANTS AND STORE IN LINKED LIST
    for(i=0;i<DIM;i++) {X[0][i] = Unit_Cell->X[b0][i];}
    B[0] = b0;
    L[0] = l0;
    N_list1_ptr = N_list1;
    for(b1=1;b1<=n_neig1;b1++) {
      for(i=0;i<DIM;i++) {X[b1][i] = N_list1_ptr->x[i];}
      B[b1] = N_list1_ptr->b;
      L[b1] = N_list1_ptr->l;
      N_list1_ptr = N_list1_ptr->next;
    }
    N_list1_ptr = N_list1;
    for(b1=1;b1<=n_neig1;b1++) {
      for(i=0;i<DIM;i++) {
        X[b1][i] = X[1][i];
        X[1][i] = N_list1_ptr->x[i];
      }
      B[b1] = B[1];
      B[1]  = N_list1_ptr->b;
      L[b1] = L[1];
      L[1]  = N_list1_ptr->l;

      //Compute the derivatives
      F->Compute(0, n_neig1+1, B, L, X, this);

      N_list1_ptr = N_list1_ptr->next;
    }


/** b0 on the end **/

    N_list1_ptr = N_list1;
    for(b1=0;b1<n_neig1;b1++) {
      //DEALLOCATE MEMROY
      delete[] B;       B      = NULL;
      delete[] L;       L      = NULL;
      delete[] mat_pos; mat_pos= NULL;
      delete[] X[0];    X[0]   = NULL;
      delete[] X;       X      = NULL;


      //BUILD NEIGHBOR LISTS
      n_neig2 = N_list1_ptr->b;
      N_list2 = BuildNeighborList_R(n_neig2, mat[n0], cutoff_);//All atoms of type defined by Mat[1] and within cutoff radius
      //Correct Higher Order Neighbor Lists And Remove Duplicates
      ManageNeighborList(n_neig2, N_list2, N_list1_ptr->b, N_list1_ptr->l, b0, l0);


      //ALLOCATE MEMORY
      n_neig2 += 3;
      mat_pos = new int[n_neig2];
      X = new double*[n_neig2];
      X[0] = new double[n_neig2*DIM];
      B = new int[n_neig2];
      L = new int*[n_neig2];
      n_neig2 -= 3;
      for(b2=0;b2<=n_neig2;b2++) {
        mat_pos[b2+1] = n0;
        X[b2+1] = X[b2] + DIM;
      }
      mat_pos[0] = n1;
      mat_pos[n_neig2+2] = -1;
      X[n_neig2+2] = NULL;
      L[n_neig2+2] = NULL;


      //COMPUTE FORCE CONSTANTS AND STORE IN LINKED LIST
      for(i=0;i<DIM;i++) {
        X[0][i] = N_list1_ptr->x[i];
        X[1][i] = Unit_Cell->X[b0][i];
      }
      B[0] = N_list1_ptr->b;
      B[1] = b0;
      L[0] = N_list1_ptr->l;
      L[1] = l0;
      N_list2_ptr = N_list2;
      for(b2=2;b2<=n_neig2+1;b2++) {
        for(i=0;i<DIM;i++) {X[b2][i] = N_list2_ptr->x[i];}
        B[b2] = N_list2_ptr->b;
        L[b2] = N_list2_ptr->l;
        N_list2_ptr = N_list2_ptr->next;
      }
      //Compute the derivatives
      F->Compute(1, n_neig2+2, B, L, X, this);

      N_list2_ptr = N_list2;
      for(b2=2;b2<=n_neig2+1;b2++) {
        for(i=0;i<DIM;i++) {
          X[b2][i] = X[1][i];
          X[1][i] = N_list2_ptr->x[i];
        }
        B[b2] = B[1];
        B[1]  = N_list2_ptr->b;
        L[b2] = L[1];
        L[1]  = N_list2_ptr->l;

        //Compute the derivatives
        F->Compute(2, n_neig2+2, B, L, X, this);

        N_list2_ptr = N_list2_ptr->next;
      }

      N_list1_ptr = N_list1_ptr->next;
      delete N_list2; N_list2=NULL;
    }  //for(b1)

    //DEALLOCATE MEMROY
    delete   N_list1; N_list1= NULL;
    delete[] B;       B      = NULL;
    delete[] L;       L      = NULL;
    delete[] mat_pos; mat_pos= NULL;
    delete[] X[0];    X[0]   = NULL;
    delete[] X;       X      = NULL;
    delete[] l0;      l0     = NULL;

  }  //for(n1)

  return;
}


/*POT_TERSOFF DESTRUCTOR: DEALLOCATES MEMORY*/
POT_TERSOFF::~POT_TERSOFF() {
  delete[] cutoff; cutoff=NULL;
  delete[] A     ; A     =NULL;
  delete[] B     ; B     =NULL;
  delete[] lambda; lambda=NULL;
  delete[] mu    ; mu    =NULL;
  delete[] beta  ; beta  =NULL;
  delete[] n     ; n     =NULL;
  delete[] c     ; c     =NULL;
  delete[] d     ; d     =NULL;
  delete[] h     ; h     =NULL;
  delete[] R     ; R     =NULL;
  return;
}
