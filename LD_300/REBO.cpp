/*						                 REBO.cpp                 						*/
/*							              01/27/2009	              						*/
/*********************************************************************
*    File that contains the POT_REBO class function definitions.  	 *
*  POT_REBO is a derivative of the POTENTIAL class.			          	 *
*********************************************************************/

/*DEFINE HEADERS*/
#include "REBO.h"
#include <cmath>


/*DEFINE SUBROUTINES*/
string Read_Next(ifstream &Input);


/*POT_REBO2 CONSTRUCTOR: STORES POTENTIAL AND POTENTIAL PARAMETER NAMES*/
POT_REBO2::POT_REBO2(int n_mat_, int *mat_) : POTENTIAL(n_mat_, mat_) {
	Pot_Symb = new char[6];
  strcpy(Pot_Symb, "REBO2");
  mat_pos = new int[2];
  cutoff = new double[n_expect];
  Q      = new double[n_expect];
  A      = new double[n_expect];
  alpha  = new double[n_expect];
  for(int i=0;i<n_expect;i++) {cutoff[i]=Q[i]=A[i]=alpha[i]=0.0;}
  return;
}


/*POT_REBO2 ReadInput: READS AND STORES POTENTIAL PARAMETER VALUES*/
string POT_REBO2::ReadInput(ifstream &Input) {
  string str;
  while(!Input.eof()) {
    str = Read_Next(Input);
    if(str.compare(0, 3, "cutoff", 3)==0) {Fill(cutoff, Input);} else
    if(str.compare(0, 1, "Q"     , 1)==0) {Fill(Q     , Input);} else
    if(str.compare(0, 1, "A"     , 1)==0) {Fill(A     , Input);} else
    if(str.compare(0, 3, "alpha" , 3)==0) {Fill(alpha , Input);} else
    {break;}
  }
  return(str);
}


/*POT_REBO2 FUNCTION GetScale: GETS THE LENGTH AND ENERGY SCALES*/
void POT_REBO2::GetScale(double &e, double &l) {
  e = 0.0;
  l = 0.0;
	for(int i=0;i<n_expect;i++) {
	  e += A[i];
	  l += Q[i] + 1.0/alpha[i];
	}
	e /= n_expect;
	l /= 2.0*n_expect;
	return;
}


/*POT_REBO2 FUNCTION Initialize: INITIALIZES THE POTENTIAL PARAMETERS*/
void POT_REBO2::Initialize() {
	for(int i=0;i<n_expect;i++) {
	  cutoff[i] /= Parameter->length;
	  A[i]      /= Parameter->energy;
	  Q[i]      /= Parameter->length;
	  alpha[i]  *= Parameter->length;
	}
	return;
}


/*POT_REBO2 FUNCTION Print: PRINTS POTENTIAL PARAMETERS*/
void POT_REBO2::Print() {
  int i;
  Log <<"\nPOTENTIAL "<<Pot_Symb<<" "<<n_mat;
  for(i=0;i<n_mat;i++) {Log <<" "<<Material->symb[mat[i]];}
  Log <<endl<<"  cutoff =";
  for(i=0;i<n_expect;i++) {Log <<" "<<cutoff[i];}
  Log <<" m\n  Q      =";
  for(i=0;i<n_expect;i++) {Log <<" "<<Q[i];}
  Log <<" m\n  A      =";
  for(i=0;i<n_expect;i++) {Log <<" "<<A[i];}
  Log <<" J\n  alpha  =";
  for(i=0;i<n_expect;i++) {Log <<" "<<alpha[i];}
  Log <<" 1/m"<<endl;
  return;
}


/*POT_REBO2 FUNCTION Copy: RETURNS A POINTER TO A NEW POTENTIAL WITH A DUPLICATE OF THE DATA*/
POTENTIAL *POT_REBO2::Copy() {
  POT_REBO2 *P_new = new POT_REBO2(n_mat, mat);
  for(int i=0;i<n_expect;i++) {
    P_new->cutoff[i] = cutoff[i];
    P_new->Q[i]      = Q[i];
    P_new->A[i]      = A[i];
    P_new->alpha[i]  = alpha[i];
  }
  return(P_new);
}


/*POT_REBO2 FUNCTION Energy: COMPUTES THE POTENTIAL ENERGY*/
double POT_REBO2::Energy(double **X) {
  int dim, i, j;
	double phi;
	double r_ij = 0.0;
	//
  i = mat_pos[0];
  j = mat_pos[1];
	j = (i>j) ? i+((j+1)*j)/2 : j+((i+1)*i)/2;
	//
	for(dim=0;dim<DIM;dim++) {
		phi = X[0][dim] - X[1][dim];
		r_ij += phi*phi;
	}
	r_ij = sqrt(r_ij);
  phi = A[j]*(1.0 + Q[j]/r_ij)*exp(-alpha[j]*r_ij);
	return(phi);
}


/*POT_REBO2::ForceConsants: COMPUTES THE HARMONIC FORCE CONSTANTS*/
void POT_REBO2::ForceConstants(int b0, POT_DER_X *F) {
  //DECLARE LOCAL VARIABLES
  int n0, n1;                 //Material's position in array mat[]
  int i;                      //Generic counter
  int b1;                     //Atom counter for 1st neighbor
  int n_neig;                 //Number of neighbors
  int *B, **L;                //List of atomic and unit cell indicies
  double cutoff_;             //Cutoff radius
  double **X;                 //List of positions
  N_LIST *N_list;             //Linked list of neighbors
  N_LIST *N_list_ptr;         //Pointer to entry in linked list
  N_LIST *BuildNeighborList_R(int &, int, double);//Computes the neighbor list


  //INITIALIZE VARIABLES
  B = new int[2];
  L = new int*[2];
  X = new double*[2];
  X[0] = new double[2*DIM];
  X[1] = X[0] + DIM;
  for(n0=0;n0<n_mat;n0++) {if(mat[n0]==Unit_Cell->mat[b0]) {break;}}
  mat_pos[0] = n0;
  for(i=0;i<DIM;i++) {X[0][i] = Unit_Cell->X[b0][i];}
  B[0] = b0;
  L[0] = NULL;


  //LOOP OVER ALL MATERIAL INTERACTIONS
  for(n1=0;n1<n_mat;n1++) {
    //GET CUTOFF
    mat_pos[1] = n1;
    cutoff_ = cutoff[(n0<=n1) ? n0*((n1+1)*n1)/2 : n1*((n0+1)*n0)/2];
    //BUILD NEIGHBOR LISTS
    n_neig = b0;
    N_list = BuildNeighborList_R(n_neig, mat[n1], cutoff_);//All atoms of type defined by Mat[1] and within cutoff radius
    //COMPUTE FORCE CONSTANTS AND STORE IN LINKED LIST
    N_list_ptr = N_list;
    for(b1=0;b1<n_neig;b1++) {
      //Assign identifiers in correct order
      for(i=0;i<DIM;i++) {X[1][i] = N_list_ptr->x[i];}
      B[1] = N_list_ptr->b;
      L[1] = N_list_ptr->l;
      //Compute derivatives (pos of b0 in lists, total number of atoms in interaction, ...)
      F->Compute(0, 2, B, L, X, this);
      //Increment neighbor list
      N_list_ptr = N_list_ptr->next;
    }
    delete   N_list;  N_list = NULL;
  }


  //DEALLOCATE MEMORY AND RETURN
  delete[] B;       B      = NULL;
  delete[] L;       L      = NULL;
  delete[] X[0];    X[0]   = NULL;
  delete[] X;       X      = NULL;
  return;
}


/*POT_REBO2 DESTRUCTOR: DEALLOCATES MEMORY*/
POT_REBO2::~POT_REBO2() {
  delete[] cutoff; cutoff=NULL;
  delete[] Q     ; Q     =NULL;
  delete[] A     ; A     =NULL;
  delete[] alpha ; alpha =NULL;
  return;
}


/*POT_REBO3 CONSTRUCTOR: STORES POTENTIAL AND POTENTIAL PARAMETER NAMES*/
POT_REBO3::POT_REBO3(int n_mat_, int *mat_) : POTENTIAL(n_mat_, mat_) {
  int i, j;
	Pot_Symb = new char[6];
  strcpy(Pot_Symb, "REBO3");
  if(n_mat!=1) {Log <<"Currently only one unique material allowed for "<<Pot_Symb<<endl;exit(0);}
  mat_pos = new int[6];
  cutoff = new double[n_expect];
  B1     = new double[n_expect];
  B2     = new double[n_expect];
  B3     = new double[n_expect];
  beta1  = new double[n_expect];
  beta2  = new double[n_expect];
  beta3  = new double[n_expect];
  Q1     = new double[n_expect];
  Q2     = new double[n_expect];
  Q3     = new double[n_expect];
  for(i=0;i<n_expect;i++) {
    cutoff[i]=B1[i]=B2[i]=B3[i]=beta1[i]=beta2[i]=beta3[i]=Q1[i]=Q2[i]=Q3[i]=0.0;
  }
  return;
}


/*POT_REBO3 ReadInput: READS AND STORES POTENTIAL PARAMETER VALUES*/
string POT_REBO3::ReadInput(ifstream &Input) {
  string str;
  while(!Input.eof()) {
    str = Read_Next(Input);
    if(str.compare(0, 3, "cutoff", 3)==0) {Fill(cutoff, Input);} else
    if(str.compare(0, 2, "B1"    , 2)==0) {Fill(B1    , Input);} else
    if(str.compare(0, 2, "B2"    , 2)==0) {Fill(B2    , Input);} else
    if(str.compare(0, 2, "B3"    , 2)==0) {Fill(B3    , Input);} else
    if(str.compare(0, 5, "beta1" , 5)==0) {Fill(beta1 , Input);} else
    if(str.compare(0, 5, "beta2" , 5)==0) {Fill(beta2 , Input);} else
    if(str.compare(0, 5, "beta3" , 5)==0) {Fill(beta3 , Input);} else
    if(str.compare(0, 2, "Q1"    , 2)==0) {Fill(Q1    , Input);} else
    if(str.compare(0, 2, "Q2"    , 2)==0) {Fill(Q2    , Input);} else
    if(str.compare(0, 2, "Q3"    , 2)==0) {Fill(Q3    , Input);} else
    {break;}
  }
  return(str);
}


/*POT_REBO3 FUNCTION GetScale: GETS THE LENGTH AND ENERGY SCALES*/
void POT_REBO3::GetScale(double &e, double &l) {
  e = 0.0;
  l = 0.0;
	for(int i=0;i<n_expect;i++) {
	  e += B1[i] + B2[i] + B3[i];
	  l += 1.0/beta1[i] + 1.0/beta2[i] + 1.0/beta3[i];
	}
	e /= 3.0*n_expect;
	l /= 3.0*n_expect;
	return;
}


/*POT_REBO3 FUNCTION Initialize: INITIALIZES THE POTENTIAL PARAMETERS*/
void POT_REBO3::Initialize() {
	for(int i=0;i<n_expect;i++) {
	  cutoff[i] /= Parameter->length;
	  B1[i]     /= Parameter->energy;
	  B2[i]     /= Parameter->energy;
	  B3[i]     /= Parameter->energy;
	  beta1[i]  *= Parameter->length;
	  beta2[i]  *= Parameter->length;
	  beta3[i]  *= Parameter->length;
	}
	return;
}


/*POT_REBO3 FUNCTION Print: PRINTS POTENTIAL PARAMETERS*/
void POT_REBO3::Print() {
  int i, j;
  Log <<"\nPOTENTIAL "<<Pot_Symb<<" "<<n_mat;
  for(i=0;i<n_mat;i++) {Log <<" "<<Material->symb[mat[i]];}
  Log <<"  \n  cutoff =";
  for(i=0;i<n_expect;i++) {Log <<" "<<cutoff[i];}
  Log <<" m\n  B1     =";
  for(i=0;i<n_expect;i++) {Log <<" "<<B1[i];}
  Log <<" J\n  B2     =";
  for(i=0;i<n_expect;i++) {Log <<" "<<B2[i];}
  Log <<" J\n  B3     =";
  for(i=0;i<n_expect;i++) {Log <<" "<<B3[i];}
  Log <<" J\n  Q1     =";
  for(i=0;i<n_expect;i++) {Log <<" "<<Q1[i];}
  Log <<"  \n  Q2     =";
  for(i=0;i<n_expect;i++) {Log <<" "<<Q2[i];}
  Log <<"  \n  Q3     =";
  for(i=0;i<n_expect;i++) {Log <<" "<<Q3[i];}
  Log <<"  \n  beta1  =";
  for(i=0;i<n_expect;i++) {Log <<" "<<beta1[i];}
  Log <<" 1/m\n  beta2  =";
  for(i=0;i<n_expect;i++) {Log <<" "<<beta2[i];}
  Log <<" 1/m\n  beta3  =";
  for(i=0;i<n_expect;i++) {Log <<" "<<beta3[i];}
  Log <<" 1/m"<<endl;
  return;
}


/*POT_REBO3 FUNCTION Copy: RETURNS A POINTER TO A NEW POTENTIAL WITH A DUPLICATE OF THE DATA*/
POTENTIAL *POT_REBO3::Copy() {
  int i, j;
  POT_REBO3 *P_new = new POT_REBO3(n_mat, mat);
  for(i=0;i<n_expect;i++) {
    P_new->cutoff[i] = cutoff[i];
    P_new->B1[i]     = B1[i];
    P_new->B1[i]     = B1[i];
    P_new->B1[i]     = B1[i];
    P_new->Q1[i]     = Q1[i];
    P_new->Q2[i]     = Q2[i];
    P_new->Q3[i]     = Q3[i];
    P_new->beta1[i]  = beta1[i];
    P_new->beta2[i]  = beta2[i];
    P_new->beta3[i]  = beta3[i];
  }
  return(P_new);
}


/*POT_REBO3 FUNCTION Energy: COMPUTES THE POTENTIAL ENERGY*/
double POT_REBO3::Energy(double **X) {
	//VARIABLE DECLARATION
	int d;
	int i, j, k;
	double r_ij = 0.0, r_ik1 = 0.0, r_ik2 = 0.0, r_jl1 = 0.0, r_jl2 = 0.0;
	double cos_ijk1 = 0.0,  cos_ijk2 = 0.0, cos_jil1 = 0.0, cos_jil2 = 0.0;
	double b_ij, b_ji;
	double G1, G2;
	double phi;


	//DETERMINE WHICH CONSTANTS ARE NEEDED
	i = mat_pos[0];
	j = mat_pos[1];
	k = mat_pos[2];
	d = (j>k) ? j+((k+1)*k)/2 : k+((j+1)*j)/2;
	j = (i>j) ? i+((j+1)*j)/2 : j+((i+1)*i)/2;
	k = (i>k) ? i+((k+1)*k)/2 : k+((i+1)*i)/2;


  //COMPUTE DISTANCES AND COSINES
	for(int d=0;d<DIM;d++) {
		phi = X[0][d]-X[1][d];
		r_ij += phi*phi;

		b_ij = X[0][d]-X[2][d];
		r_ik1 += b_ij*b_ij;
		cos_ijk1 += phi*b_ij;

		G1 = X[0][d]-X[3][d];
		r_ik2 += G1*G1;
		cos_ijk2 += phi*G1;

		b_ji = X[1][d]-X[4][d];
		r_jl1 += b_ji*b_ji;
		cos_jil1 -= phi*b_ji;

		G2 = X[1][d]-X[5][d];
		r_jl2 += G2*G2;
		cos_jil2 -= phi*G2;
	}
	r_ij  = sqrt(r_ij);
	r_ik1 = sqrt(r_ik1);
	r_ik2 = sqrt(r_ik2);
	r_jl1 = sqrt(r_jl1);
	r_jl2 = sqrt(r_jl2);
	cos_ijk1 /= r_ij*r_ik1;
	cos_ijk2 /= r_ij*r_ik2;
	cos_jil1 /= r_ij*r_jl1;
	cos_jil2 /= r_ij*r_jl2;


	//COMPUTE REBO3 ENERGY
	G1 = (Q1[j]*cos_ijk1 + Q2[j])*cos_ijk1 + Q3[j];
	G2 = (Q1[j]*cos_ijk2 + Q2[j])*cos_ijk2 + Q3[j];
	b_ij = 1.0/sqrt(1 + G1 + G2);
	G1 = (Q1[j]*cos_jil1 + Q2[j])*cos_jil1 + Q3[j];
	G2 = (Q1[j]*cos_jil2 + Q2[j])*cos_jil2 + Q3[j];
	b_ji = 1.0/sqrt(1 + G1 + G2);

	G1 = -0.5*(b_ij + b_ji);
	G2 = B1[j]*exp(-beta1[j]*r_ij) + B2[j]*exp(-beta2[j]*r_ij) + B3[j]*exp(-beta3[j]*r_ij);

	phi = G1*G2;

	return(phi);
}


/*POT_REBO3::ForceConsants: COMPUTES THE HARMONIC FORCE CONSTANTS*/
void POT_REBO3::ForceConstants(int b0, POT_DER_X *F) {
  //DECLARE LOCAL VARIABLES
  int n0, n1, n2;             //Material's position in array mat[]
  int b1, b2, i;              //Counters
  int n_neig1, *n_neig2, **n_neig3;//Number of neighbors in lists
  int *B, **L;                //Stores atom and unit cell numbers
  int *l0;                    //Array of zeros (center unit cell)
  double cutoff_;             //Stores proper cutoff
  double **X;                 //Array of atomic positions used in potential
  N_LIST *N_list1;            //Linked list of neighbors to b0
  N_LIST *N_list1_ptr;        //Pointer to neighbors of b0
  N_LIST **N_list2;           //Linked list of neighbors to b1
  N_LIST *N_list2_ptr;        //Pointer to neighbors of b1
  N_LIST ***N_list3;          //Linked list of neighbors to b2
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


/**mat[n0] AT VERTEX**/

    //BUILD FIRST NEIGHBOR LIST
    cutoff_ = cutoff[(n0<=n1) ? n0*((n1+1)*n1)/2 : n1*((n0+1)*n0)/2];
    n_neig1 = b0;
    N_list1 = BuildNeighborList_R(n_neig1, mat[n1], cutoff_);
    if(n_neig1!=3) {Log <<"REBO potential found "<<n_neig1<<" neighbors for atom "<<b0<<"."<<endl;exit(0);}
    N_list1->next->next->next = N_list1;//Make list cyclic (RELEASE BEFORE DEALLOCATING!!)

    for(n2=n1;n2<n_mat;n2++) {
      //BUILD SECOND NEIGHBOR LIST
      cutoff_ = cutoff[(n0<=n2) ? n0*((n2+1)*n2)/2 : n2*((n0+1)*n0)/2];
      n_neig2 = new int[n_neig1];
      N_list2 = new N_LIST*[n_neig1];
      N_list1_ptr = N_list1;
      for(b1=0;b1<n_neig1;b1++) {
        n_neig2[b1] = N_list1_ptr->b;
        N_list2[b1] = BuildNeighborList_R(n_neig2[b1], mat[n2], cutoff_);
        if(n_neig2[b1]!=3) {Log <<"REBO potential found "<<n_neig2[b1]<<" neighbors for atom "<<N_list1_ptr->b<<"."<<endl;exit(0);}
        N_list1_ptr = N_list1_ptr->next;
      }
      //Correct Higher Order Neighbor Lists And Remove Duplicates
      N_list1_ptr = N_list1;
      for(b1=0;b1<n_neig1;b1++) {
        ManageNeighborList(n_neig2[b1], N_list2[b1], N_list1_ptr->b,
          N_list1_ptr->l, b0, l0);
        //Increment neighbor list
        N_list1_ptr = N_list1_ptr->next;
      }
      //Make list cyclic (RELEASE BEFORE DEALLOCATING!!)
      for(b1=0;b1<n_neig1;b1++) {N_list2[b1]->next->next = N_list2[b1];}


      //COMPUTE FORCE CONSTANTS AND STORE IN LINKED LIST
      mat_pos[0] = n0;
      mat_pos[1] = n1;
      mat_pos[2] = n2;
      B[0] = b0;
      L[0] = l0;
      for(i=0;i<DIM;i++) {X[0][i] = Unit_Cell->X[b0][i];}
      N_list1_ptr = N_list1;
      for(b1=0;b1<n_neig1;b1++) {
        //Assign atomic indentifiers
        B[1] = N_list1_ptr->b;
        B[2] = N_list1_ptr->next->b;
        B[3] = N_list1_ptr->next->next->b;
        B[4] = N_list2[b1]->b;
        B[5] = N_list2[b1]->next->b;
        L[1] = N_list1_ptr->l;
        L[2] = N_list1_ptr->next->l;
        L[3] = N_list1_ptr->next->next->l;
        L[4] = N_list2[b1]->l;
        L[5] = N_list2[b1]->next->l;
        for(i=0;i<DIM;i++) {
          X[1][i] = N_list1_ptr->x[i];
          X[2][i] = N_list1_ptr->next->x[i];
          X[3][i] = N_list1_ptr->next->next->x[i];
          X[4][i] = N_list2[b1]->x[i];
          X[5][i] = N_list2[b1]->next->x[i];
        }
        F->Compute(0, 6, B, L, X, this);
        N_list1_ptr = N_list1_ptr->next;
      }


/**mat[n0] ON END**/

      //BUILD THIRD NEIGHBOR LIST
      cutoff_ = cutoff[(n0<=n2) ? n0*((n2+1)*n2)/2 : n2*((n0+1)*n0)/2];
      n_neig3 = new int*[n_neig1];
      N_list3 = new N_LIST**[n_neig1];
      b2 = 0;
      for(b1=0;b1<n_neig1;b1++) {b2 += n_neig2[b1];}
      n_neig3[0] = new int[b2];
      N_list3[0] = new N_LIST*[b2];
      for(b1=0;b1<n_neig1;b1++) {
        if(b1+1<n_neig1) {
          n_neig3[b1+1] = n_neig3[b1] + n_neig2[b1];
          N_list3[b1+1] = N_list3[b1] + n_neig2[b1];
        }
        N_list2_ptr = N_list2[b1];
        for(b2=0;b2<n_neig2[b1];b2++) {
          n_neig3[b1][b2] = N_list2_ptr->b;
          N_list3[b1][b2] = BuildNeighborList_R(n_neig3[b1][b2], mat[n2], cutoff_);//All atoms of type defined by Mat[2] and within cutoff radius
          if(n_neig3[b1][b2]!=3) {Log <<"REBO potential found "<<n_neig3[b1][b2]<<" neighbors for atom "<<N_list1_ptr->b<<"."<<endl;exit(0);}
          N_list2_ptr = N_list2_ptr->next;
        }
      }

      //Correct Higher Order Neighbor Lists And Remove Duplicates
      N_list1_ptr = N_list1;
      for(b1=0;b1<n_neig1;b1++) {
        N_list2_ptr = N_list2[b1];
        for(b2=0;b2<n_neig2[b1];b2++) {
          ManageNeighborList(n_neig3[b1][b2], N_list3[b1][b2],
            N_list2_ptr->b, N_list2_ptr->l, N_list1_ptr->b, N_list1_ptr->l);
          //Increment neighbor list
          N_list2_ptr = N_list2_ptr->next;
        }
        //Increment neighbor list
        N_list1_ptr = N_list1_ptr->next;
      }
      //Make list cyclic (RELEASE BEFORE DEALLOCATING!!)
      for(b1=0;b1<n_neig1;b1++) {
        for(b2=0;b2<n_neig2[b1];b2++) {N_list3[b1][b2]->next->next = N_list3[b1][b2];}
      }


      //COMPUTE FORCE CONSTANTS AND STORE IN LINKED LIST
      mat_pos[0] = n0;
      mat_pos[1] = n1;
      mat_pos[2] = n2;
      B[5] = b0;
      L[5] = l0;
      for(i=0;i<DIM;i++) {X[5][i] = Unit_Cell->X[b0][i];}
      N_list1_ptr = N_list1;
      for(b1=0;b1<n_neig1;b1++) {
        N_list2_ptr = N_list2[b1];
        for(b2=0;b2<n_neig2[b1];b2++) {
          //Assign atomic indentifiers
          B[0] = N_list2_ptr->b;
          B[1] = N_list1_ptr->b;
          B[2] = N_list3[b1][b2]->b;
          B[3] = N_list3[b1][b2]->next->b;
          B[4] = N_list2_ptr->next->b;
          L[0] = N_list2_ptr->l;
          L[1] = N_list1_ptr->l;
          L[2] = N_list3[b1][b2]->l;
          L[3] = N_list3[b1][b2]->next->l;
          L[4] = N_list2_ptr->next->l;
          for(i=0;i<DIM;i++) {
            X[0][i] = N_list2_ptr->x[i];
            X[1][i] = N_list1_ptr->x[i];
            X[2][i] = N_list3[b1][b2]->x[i];
            X[3][i] = N_list3[b1][b2]->next->x[i];
            X[4][i] = N_list2_ptr->next->x[i];
          }
          F->Compute(5, 6, B, L, X, this);
          //Increment neighbor list
          N_list2_ptr = N_list2_ptr->next;
        }
        //Increment neighbor list
        N_list1_ptr = N_list1_ptr->next;
      }
      //Deallocate second and third neighbor lists
      for(b1=0;b1<n_neig1;b1++) {
        for(b2=0;b2<n_neig2[b1];b2++) {N_list3[b1][b2]->next->next=NULL;}
      }
      for(b1=0;b1<n_neig1;b1++) {
        N_list2[b1]->next->next = NULL;   //Break cyclic list
      }
      delete[] n_neig3[0]; n_neig3[0]=NULL;
      delete[] N_list3[0]; N_list3[0]=NULL;
      delete[] n_neig3;    n_neig3   =NULL;
      delete[] N_list3;    N_list3   =NULL;
      delete[] n_neig2;    n_neig2   =NULL;
      delete[] N_list2;    N_list2   =NULL;

    } //for(n2)
    //Deallocate first neighbor list
    N_list1->next->next->next = NULL;   //Break cyclic list
    delete   N_list1;    N_list1   =NULL;

  } //for(n1)


  //DEALLOCATE MEMORY AND RETURN
  delete[] B;          B         =NULL;
  delete[] L;          L         =NULL;
  delete[] X[0];       X[0]      =NULL;
  delete[] X;          X         =NULL;
  delete[] l0;         l0        =NULL;
  return;
}


/*POT_REBO3 DESTRUCTOR: DEALLOCATES MEMORY*/
POT_REBO3::~POT_REBO3() {
  delete[] cutoff; cutoff=NULL;
  delete[] B1    ; B1    =NULL;
  delete[] B2    ; B2    =NULL;
  delete[] B3    ; B3    =NULL;
  delete[] beta1 ; beta1 =NULL;
  delete[] beta2 ; beta2 =NULL;
  delete[] beta3 ; beta3 =NULL;
  delete[] Q1    ; Q1    =NULL;
  delete[] Q2    ; Q2    =NULL;
  delete[] Q3    ; Q3    =NULL;
  return;
}
