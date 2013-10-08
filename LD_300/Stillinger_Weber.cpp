/*						           Stillinger_Weber.cpp           						*/
/*							              12/02/2008	              						*/
/*********************************************************************
*    File that contains the POT_SW class function definitions.	  	 *
*  POT_SW is a derivative of the POTENTIAL class.				          	 *
*********************************************************************/

/*DEFINE HEADERS*/
#include "Stillinger_Weber.h"
#include <cmath>


/*DEFINE SUBROUTINES*/
string Read_Next(ifstream &Input);


/*POT_SW2 CONSTRUCTOR: STORES POTENTIAL AND POTENTIAL PARAMETER NAMES*/
POT_SW2::POT_SW2(int n_mat_, int *mat_) : POTENTIAL(n_mat_, mat_) {
	Pot_Symb = new char[4];
  strcpy(Pot_Symb, "SW2");
  mat_pos = new int[2];
  epsilon = new double[n_expect];
  sigma   = new double[n_expect];
  a       = new double[n_expect];
  A       = new double[n_expect];
  B       = new double[n_expect];
  p       = new double[n_expect];
  q       = new double[n_expect];
  for(int i=0;i<n_expect;i++) {
    epsilon[i]=sigma[i]=a[i]=A[i]=B[i]=p[i]=q[i]=0.0;
  }
  return;
}


/*POT_SW2 ReadInput: READS AND STORES POTENTIAL PARAMETER VALUES*/
string POT_SW2::ReadInput(ifstream &Input) {
  string str;
  while(!Input.eof()) {
    str = Read_Next(Input);
    if(str.compare(0, 7, "epsilon", 7)==0) {Fill(epsilon, Input);} else
    if(str.compare(0, 5, "sigma"  , 5)==0) {Fill(sigma  , Input);} else
    if(str.compare(0, 1, "a"      , 1)==0) {Fill(a      , Input);} else
    if(str.compare(0, 1, "A"      , 1)==0) {Fill(A      , Input);} else
    if(str.compare(0, 1, "B"      , 1)==0) {Fill(B      , Input);} else
    if(str.compare(0, 1, "p"      , 1)==0) {Fill(p      , Input);} else
    if(str.compare(0, 1, "q"      , 1)==0) {Fill(q      , Input);} else
    {break;}
  }
  return(str);
}


/*POT_SW2 FUNCTION GetScale: GETS THE LENGTH AND ENERGY SCALES*/
void POT_SW2::GetScale(double &e, double &l) {
  e = 0.0;
  l = 0.0;
	for(int i=0;i<n_expect;i++) {
	  e += epsilon[i];
	  l += sigma[i];
	}
	e /= n_expect;
	l /= n_expect;
	return;
}


/*POT_SW2 FUNCTION Initialize: INITIALIZES THE POTENTIAL PARAMETERS*/
void POT_SW2::Initialize() {
	for(int i=0;i<n_expect;i++) {
	  epsilon[i] /= Parameter->energy;
	  sigma[i]   /= Parameter->length;
	  a[i]       *= sigma[i];
	}
	return;
}


/*POT_SW2 FUNCTION Print: PRINTS POTENTIAL PARAMETERS*/
void POT_SW2::Print() {
  int i;
  Log <<"\nPOTENTIAL "<<Pot_Symb<<" "<<n_mat;
  for(i=0;i<n_mat;i++) {Log <<" "<<Material->symb[mat[i]];}
  Log <<endl<<"  epsilon =";
  for(i=0;i<n_expect;i++) {Log <<" "<<epsilon[i];}
  Log <<" J\n  sigma   =";
  for(i=0;i<n_expect;i++) {Log <<" "<<sigma[i];}
  Log <<" m\n  a       =";
  for(i=0;i<n_expect;i++) {Log <<" "<<a[i];}
  Log <<"  \n  A       =";
  for(i=0;i<n_expect;i++) {Log <<" "<<A[i];}
  Log <<"  \n  B       =";
  for(i=0;i<n_expect;i++) {Log <<" "<<B[i];}
  Log <<"  \n  p       =";
  for(i=0;i<n_expect;i++) {Log <<" "<<p[i];}
  Log <<"  \n  q       =";
  for(i=0;i<n_expect;i++) {Log <<" "<<q[i];}
  Log <<endl;
  return;
}


/*POT_SW2 FUNCTION Copy: RETURNS A POINTER TO A NEW POTENTIAL WITH A DUPLICATE OF THE DATA*/
POTENTIAL *POT_SW2::Copy() {
  POT_SW2 *P_new = new POT_SW2(n_mat, mat);
  for(int i=0;i<n_expect;i++) {
    P_new->epsilon[i] = epsilon[i];
    P_new->sigma[i]   = sigma[i];
    P_new->a[i]       = a[i];
    P_new->A[i]       = A[i];
    P_new->B[i]       = B[i];
    P_new->p[i]       = p[i];
    P_new->q[i]       = q[i];
  }
  return(P_new);
}


/*POT_SW2 FUNCTION Energy: COMPUTES THE POTENTIAL ENERGY*/
double POT_SW2::Energy(double **X) {
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
	if(r_ij>=a[j]) {return(0.0);}
	r_ij /= sigma[j];
  phi = epsilon[j]*A[j]*(B[j]/pow(r_ij,p[j]) - 1.0/pow(r_ij,q[j]))*exp(sigma[j]/(r_ij*sigma[j]-a[j]));
	return(phi);
}


/*POT_SW2::ForceConsants: COMPUTES THE HARMONIC FORCE CONSTANTS*/
void POT_SW2::ForceConstants(int b0, POT_DER_X *F) {
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
    cutoff_ = a[(n0<=n1) ? n0*((n1+1)*n1)/2 : n1*((n0+1)*n0)/2];
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


/*POT_SW2 DESTRUCTOR: DEALLOCATES MEMORY*/
POT_SW2::~POT_SW2() {
  delete[] epsilon; epsilon=NULL;
  delete[] sigma  ; sigma  =NULL;
  delete[] a      ; a      =NULL;
  delete[] A      ; A      =NULL;
  delete[] B      ; B      =NULL;
  delete[] p      ; p      =NULL;
  delete[] q      ; q      =NULL;
  return;
}


/*POT_SW3 CONSTRUCTOR: STORES POTENTIAL AND POTENTIAL PARAMETER NAMES*/
POT_SW3::POT_SW3(int n_mat_, int *mat_) : POTENTIAL(n_mat_, mat_) {
  int i, j;
	Pot_Symb = new char[4];
  strcpy(Pot_Symb, "SW3");
  mat_pos = new int[3];
  epsilon = new double*[n_expect];
  sigma   = new double[n_expect];
  a       = new double[n_expect];
  gamma   = new double[n_expect];
  lambda  = new double*[n_expect];
  epsilon[0] = new double[n_expect*n_mat];
  lambda[0]  = new double[n_expect*n_mat];
  for(i=0;i<n_mat;i++) {
    if(i<n_mat-1) {
      epsilon[i+1] = epsilon[i] + n_expect;
      lambda[i+1]  = lambda[i]  + n_expect;
    }
  }
  for(i=0;i<n_expect;i++) {
    for(j=0;j<n_mat;j++) {epsilon[j][i]=lambda[j][i]=0.0;}
    sigma[i]=a[i]=gamma[i]=0.0;
  }
  return;
}


/*POT_SW3 ReadInput: READS AND STORES POTENTIAL PARAMETER VALUES*/
string POT_SW3::ReadInput(ifstream &Input) {
  int n_e=0, n_l=0;
  string str;
  while(!Input.eof()) {
    str = Read_Next(Input);
    if(str.compare(0, 7, "epsilon", 7)==0) {Fill(epsilon[n_e], Input);n_e++;} else
    if(str.compare(0, 5, "sigma"  , 5)==0) {Fill(sigma  , Input);} else
    if(str.compare(0, 1, "a"      , 1)==0) {Fill(a      , Input);} else
    if(str.compare(0, 5, "gamma"  , 5)==0) {Fill(gamma  , Input);} else
    if(str.compare(0, 6, "lambda" , 6)==0) {Fill(lambda[n_l], Input);n_l++;} else
    {break;}
  }
  return(str);
}


/*POT_SW3 FUNCTION GetScale: GETS THE LENGTH AND ENERGY SCALES*/
void POT_SW3::GetScale(double &e, double &l) {
  int i, j;
  e = 0.0;
  l = 0.0;
	for(i=0;i<n_expect;i++) {
	  for(j=0;j<n_mat;j++) {e += epsilon[j][i];}
	  l += sigma[i];
	}
	e /= n_expect*n_mat;
	l /= n_expect;
	return;
}


/*POT_SW3 FUNCTION Initialize: INITIALIZES THE POTENTIAL PARAMETERS*/
void POT_SW3::Initialize() {
  int i, j;
	for(i=0;i<n_expect;i++) {
	  for(j=0;j<n_mat;j++) {epsilon[j][i] /= Parameter->energy;}
	  sigma[i]   /= Parameter->length;
	  a[i]       *= sigma[i];
	}
	return;
}


/*POT_SW3 FUNCTION Print: PRINTS POTENTIAL PARAMETERS*/
void POT_SW3::Print() {
  int i, j;
  Log <<"\nPOTENTIAL "<<Pot_Symb<<" "<<n_mat;
  for(i=0;i<n_mat;i++) {Log <<" "<<Material->symb[mat[i]];}
  for(j=0;j<n_mat;j++) {
    Log <<endl<<"  epsilon =";
    for(i=0;i<n_expect;i++) {Log <<" "<<epsilon[j][i];}
    Log <<" J  %for "<<Material->symb[mat[j]]<<" at center";
  }
  Log <<endl<<"  sigma   =";
  for(i=0;i<n_expect;i++) {Log <<" "<<sigma[i];}
  Log <<" m\n  a       =";
  for(i=0;i<n_expect;i++) {Log <<" "<<a[i];}
  Log <<"  \n  gamma   =";
  for(i=0;i<n_expect;i++) {Log <<" "<<gamma[i];}
  for(j=0;j<n_mat;j++) {
    Log <<" \n  lambda  =";
    for(i=0;i<n_expect;i++) {Log <<" "<<lambda[j][i];}
    Log <<"  %for "<<Material->symb[mat[j]]<<" at center";
  }
  Log <<endl;
  return;
}


/*POT_SW3 FUNCTION Copy: RETURNS A POINTER TO A NEW POTENTIAL WITH A DUPLICATE OF THE DATA*/
POTENTIAL *POT_SW3::Copy() {
  int i, j;
  POT_SW3 *P_new = new POT_SW3(n_mat, mat);
  for(i=0;i<n_expect;i++) {
    for(j=0;j<n_mat;j++) {
      P_new->epsilon[j][i] = epsilon[j][i];
      P_new->lambda[j][i]  = lambda[j][i];
    }
    P_new->sigma[i]   = sigma[i];
    P_new->a[i]       = a[i];
    P_new->gamma[i]   = gamma[i];
  }
  return(P_new);
}


/*POT_SW3 FUNCTION Energy: COMPUTES THE POTENTIAL ENERGY*/
double POT_SW3::Energy(double **X) {
	//VARIABLE DECLARATION
	int dim;
	int i, j, k;
	double eps, lam;
	double r_ij = 0.0, r_ik = 0.0;
	double cos_theta = 0.0;
	double e;
	double phi;


	//DETERMINE WHICH CONSTANTS ARE NEEDED
	i = mat_pos[0];
	j = mat_pos[1];
	k = mat_pos[2];
	dim = (j>k) ? j+((k+1)*k)/2 : k+((j+1)*j)/2;
	eps = epsilon[i][dim];
	lam = lambda[i][dim];
	j = (i>j) ? i+((j+1)*j)/2 : j+((i+1)*i)/2;
	k = (i>k) ? i+((k+1)*k)/2 : k+((i+1)*i)/2;


	//COMPUTE STILLINGER-WEBER THREE BODY TERM
	for(dim=0;dim<DIM;dim++) {
		phi = X[0][dim]-X[1][dim];
		r_ij += phi*phi;
		e = X[0][dim]-X[2][dim];
		r_ik += e*e;
		cos_theta += phi*e;
	}
	r_ij = sqrt(r_ij);
	r_ik = sqrt(r_ik);
	if((r_ij>=a[j])||(r_ik>=a[k])) {return(0.0);}
	cos_theta /= r_ij*r_ik;
  r_ij = gamma[j]*sigma[j]/(r_ij-a[j]);
  r_ik = gamma[k]*sigma[k]/(r_ik-a[k]);
  e = lam * exp(r_ij + r_ik);
  phi = cos_theta + 0.33333333333333333;
  phi = eps * e * phi * phi;

	return(phi);
}


/*POT_SW3::ForceConsants: COMPUTES THE HARMONIC FORCE CONSTANTS*/
void POT_SW3::ForceConstants(int b0, POT_DER_X *F) {
  //DECLARE LOCAL VARIABLES
  int n0, n1, n2;             //Material's position in array mat[]
  int b1, b2, i, end;         //Counters
  int n_neig1, *n_neig2;      //Number of neighbors in lists
  int *B, **L;                //Stores atom and unit cell numbers
  int *l0;                    //Array of zeros (center unit cell)
  double cutoff_;             //Stores proper cutoff
  double **X;                 //Array of atomic positions used in potential
  N_LIST *N_list1;            //Linked list of neighbors to b0
  N_LIST *N_list1_ptr;        //Pointer to neighbors of b0
  N_LIST **N_list2;           //Linked list of neighbors to b1
  N_LIST *N_list2_ptr;        //Pointer to neighbors of b1
  N_LIST *BuildNeighborList_R(int &, int, double);
  void ManageNeighborList(int&, N_LIST*, int, int*, int, int*);


  //INITIALIZE VARIABLES
  B = new int[3];
  L = new int*[3];
  X = new double*[3];
  X[0] = new double[3*DIM];
  X[1] = X[0] + DIM;
  X[2] = X[1] + DIM;
  l0 = new int[DIM];
  for(i=0;i<DIM;i++) {l0[i] = 0;}
  for(n0=0;n0<n_mat;n0++) {if(mat[n0]==Unit_Cell->mat[b0]) {break;}}


  //LOOP OVER ALL MATERIAL INTERACTIONS
  for(n1=0;n1<n_mat;n1++) {


/**mat[n0] AT VERTEX**/

    //BUILD FIRST NEIGHBOR LIST
    cutoff_ = a[(n0<=n1) ? n0*((n1+1)*n1)/2 : n1*((n0+1)*n0)/2];
    n_neig1 = b0;
    N_list1 = BuildNeighborList_R(n_neig1, mat[n1], cutoff_);//All atoms of type defined by Mat[1] and within cutoff radius

    for(n2=n1;n2<n_mat;n2++) {
      //BUILD SECOND NEIGHBOR LIST
      mat_pos[2] = n2;
      cutoff_ = a[(n0<=n2) ? n0*((n2+1)*n2)/2 : n2*((n0+1)*n0)/2];
      n_neig2 = new int;
      N_list2 = new N_LIST*;
      if(n1==n2) {
        n_neig2[0] = n_neig1;
        N_list2[0] = N_list1;
      }
      else {
        n_neig2[0] = b0;
        N_list2[0] = BuildNeighborList_R(n_neig2[0], mat[n2], cutoff_);//All atoms of type defined by Mat[2] and within cutoff radius
      }


      //COMPUTE FORCE CONSTANTS AND STORE IN LINKED LIST
      mat_pos[0] = n0;
      mat_pos[1] = n1;
      for(i=0;i<DIM;i++) {X[0][i] = Unit_Cell->X[b0][i];}
      B[0] = b0;
      L[0] = l0;
      N_list1_ptr = N_list1;
      for(b1=0;b1<n_neig1;b1++) {
        //Assign identifiers in correct order
        for(i=0;i<DIM;i++) {X[1][i] = N_list1_ptr->x[i];}
        B[1] = N_list1_ptr->b;
        L[1] = N_list1_ptr->l;
        N_list2_ptr = N_list2[0];
        if(N_list1==N_list2[0]) {end = b1;}
        else {end = n_neig2[0];}
        for(b2=0;b2<end;b2++) {
          //Assign identifiers in correct order
          for(i=0;i<DIM;i++) {X[2][i] = N_list2_ptr->x[i];}
          B[2] = N_list2_ptr->b;
          L[2] = N_list2_ptr->l;
          //Compute derivatives (pos of b0 in lists, total number of unique atoms, ...)
          F->Compute(0, 3, B, L, X, this);
          //Increment neighbor list
          N_list2_ptr = N_list2_ptr->next;
        }
        //Increment neighbor list
        N_list1_ptr = N_list1_ptr->next;
      }
      if(n1!=n2) {
        delete N_list2;
        delete n_neig2;
      }

    }


/**mat[n0] ON END**/
    for(n2=0;n2<n_mat;n2++) {
      //BUILD NEIGHBOR LISTS
      cutoff_ = a[(n0<=n2) ? n0*((n2+1)*n2)/2 : n2*((n0+1)*n0)/2];
      N_list2 = new N_LIST*[n_neig1];
      n_neig2 = new int[n_neig1];
      N_list1_ptr = N_list1;
      for(b1=0;b1<n_neig1;b1++) {
        n_neig2[b1] = N_list1_ptr->b;
        N_list2[b1] = BuildNeighborList_R(n_neig2[b1], mat[n2], cutoff_);//All atoms of type defined by Mat[2] and within cutoff radius
        N_list1_ptr = N_list1_ptr->next;
      }
      //Correct Higher Order Neighbor Lists And Remove Duplicates
      N_list1_ptr = N_list1;
      for(b1=0;b1<n_neig1;b1++) {
        ManageNeighborList(n_neig2[b1], N_list2[b1], N_list1_ptr->b, N_list1_ptr->l, b0, l0);
        //Increment neighbor list
        N_list1_ptr = N_list1_ptr->next;
      }


      //COMPUTE FORCE CONSTANTS AND STORE IN LINKED LIST
      //NOTE: we have to swap X[0] and X[1] to get center atom right,
      mat_pos[0] = n1;
      mat_pos[1] = n0;
      for(i=0;i<DIM;i++) {X[1][i] = Unit_Cell->X[b0][i];}
      B[1] = b0;
      L[1] = l0;
      N_list1_ptr = N_list1;
      for(b1=0;b1<n_neig1;b1++) {
        //Assign identifiers in correct order
        for(i=0;i<DIM;i++) {X[0][i] = N_list1_ptr->x[i];}
        B[0] = N_list1_ptr->b;
        L[0] = N_list1_ptr->l;
        N_list2_ptr = N_list2[b1];
        for(b2=0;b2<n_neig2[b1];b2++) {
          //Assign identifiers in correct order
          for(i=0;i<DIM;i++) {X[2][i] = N_list2_ptr->x[i];}
          B[2] = N_list2_ptr->b;
          L[2] = N_list2_ptr->l;
          //Compute derivatives (pos of b0 in lists, total number of unique atoms, ...)
          F->Compute(1, 3, B, L, X, this);
          //Increment neighbor list
          N_list2_ptr = N_list2_ptr->next;
        }
        //Increment neighbor list
        N_list1_ptr = N_list1_ptr->next;
      }

      delete[] n_neig2;
      for(b1=0;b1<n_neig1;b1++) {delete N_list2[b1]; N_list2[b1]=NULL;}
      delete[] N_list2; N_list2=NULL;
    }  //end for(n2)
    delete N_list1;
  }  //end for(n1)


  //DEALLOCATE MEMORY AND RETURN
  delete[] B;       B       = NULL;
  delete[] L;       L       = NULL;
  delete[] X[0];    X[0]    = NULL;
  delete[] X;       X       = NULL;
  delete[] l0;      l0      = NULL;
  return;
}


/*POT_SW3 DESTRUCTOR: DEALLOCATES MEMORY*/
POT_SW3::~POT_SW3() {
  delete[] epsilon[0]; epsilon[0]=NULL;
  delete[] lambda[0];  lambda[0] =NULL;
  delete[] epsilon; epsilon=NULL;
  delete[] sigma  ; sigma  =NULL;
  delete[] a      ; a      =NULL;
  delete[] gamma  ; gamma  =NULL;
  delete[] lambda ; lambda =NULL;
  return;
}
