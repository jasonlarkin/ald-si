/*						             Lennard_Jones.cpp		            				*/
/*							              12/02/2008		              					*/
/*********************************************************************
*    File that contains the POT_LJ class function definitions.	  	 *
*  POT_LJ is a derivative of the POTENTIAL class.				          	 *
*********************************************************************/

/*DEFINE HEADERS*/
#include "Lennard_Jones.h"


/*DEFINE SUBROUTINES*/
string Read_Next(ifstream &Input);


/*POT_LJ CONSTRUCTOR: STORES POTENTIAL AND POTENTIAL PARAMETER NAMES*/
POT_LJ::POT_LJ(int n_mat_, int *mat_) : POTENTIAL(n_mat_, mat_) {
  Pot_Symb = new char[3];
  strcpy(Pot_Symb, "LJ");
  mat_pos = new int[2];
  epsilon = new double[n_expect];
  sigma   = new double[n_expect];
  cutoff  = new double[n_expect];
  Acutoff = new double[n_expect];
  for(int i=0;i<n_expect;i++) {
    epsilon[i]=sigma[i]=cutoff[i]=Acutoff[i]=0.0;
  }
  return;
}


/*POT_LJ ReadInput: READS AND STORES POTENTIAL PARAMETER VALUES*/
string POT_LJ::ReadInput(ifstream &Input) {
  string str;
  while(!Input.eof()) {
    str = Read_Next(Input);
    if(str.compare(0, 7, "epsilon", 7)==0) {Fill(epsilon, Input);} else
    if(str.compare(0, 5, "sigma"  , 5)==0) {Fill(sigma  , Input);} else
    if(str.compare(0, 3, "cutoff" , 3)==0) {Fill(cutoff , Input);} else
    if(str.compare(0, 3, "Acutoff", 3)==0) {Fill(Acutoff, Input);} else
    {break;}
  }
  return(str);
}


/*POT_LJ FUNCTION GetScale: INITIALIZES INITIALIZES THE LENGTH AND ENERGY SCALES*/
void POT_LJ::GetScale(double &e, double &l) {
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


/*POT_LJ FUNCTION Initialize: INITIALIZES THE POTENTIAL PARAMETERS*/
void POT_LJ::Initialize() {
	for(int i=0;i<n_expect;i++) {
	  if(Acutoff[i]<=0) {Acutoff[i] = cutoff[i];}
	  epsilon[i] /= Parameter->energy;
	  sigma[i]   /= Parameter->length;
	  cutoff[i]  /= Parameter->length;
	  Acutoff[i] /= Parameter->length;
	}
	return;
}


/*POT_LJ FUNCTION Print: PRINTS POTENTIAL PARAMETERS*/
void POT_LJ::Print() {
  int i;
  Log <<"\nPOTENTIAL "<<Pot_Symb<<" "<<n_mat;
  for(i=0;i<n_mat;i++) {Log <<" "<<Material->symb[mat[i]];}
  Log <<endl<<"  epsilon =";
  for(i=0;i<n_expect;i++) {Log <<" "<<epsilon[i];}
  Log <<" J\n  sigma   =";
  for(i=0;i<n_expect;i++) {Log <<" "<<sigma[i];}
  Log <<" m\n  cutoff  =";
  for(i=0;i<n_expect;i++) {Log <<" "<<cutoff[i];}
  Log <<" m\n  Acutoff =";
  for(i=0;i<n_expect;i++) {Log <<" "<<Acutoff[i];}
  Log <<" m\n"<<endl;
  return;
}


/*POT_LJ FUNCTION Copy: RETURNS A POINTER TO A NEW POTENTIAL WITH A DUPLICATE OF THE DATA*/
POTENTIAL *POT_LJ::Copy() {
  POT_LJ *P_new = new POT_LJ(n_mat, mat);
  for(int i=0;i<n_expect;i++) {
    P_new->epsilon[i] = epsilon[i];
    P_new->sigma[i]   = sigma[i];
    P_new->cutoff[i]  = cutoff[i];
    P_new->Acutoff[i] = Acutoff[i];
  }
  return(P_new);
}


/*POT_LJ FUNCTION Energy: COMPUTES THE POTENTIAL ENERGY*/
double POT_LJ::Energy(double **X) {
  int i, j;
	double r_ij2 = 0.0;
	double r_ij6;
	//
	i = mat_pos[0];
	j = mat_pos[1];
	j = (i>j) ? i+((j+1)*j)/2 : j+((i+1)*i)/2;
	//
	for(i=0;i<DIM;i++) {
		r_ij6 = X[0][i] - X[1][i];
		r_ij2 += r_ij6*r_ij6;
	}
	r_ij2 = r_ij2/(sigma[j]*sigma[j]);
	//if(r_ij2>=cutoff[j]*cutoff[j]) {return(0.0);}//By doing this here we might get weird derivatives
	r_ij6 = 1.0/(r_ij2*r_ij2*r_ij2);
	double phi = 4.0 * epsilon[j] * ( r_ij6*r_ij6 - r_ij6 );

	return(phi);
}


/*POT_LJ::ForceConsants: COMPUTES THE HARMONIC FORCE CONSTANTS*/
void POT_LJ::ForceConstants(int b0, POT_DER_X *F) {
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
    if(typeid(F[0])==typeid(PD_ANHARMONIC)) {
      cutoff_ =Acutoff[(n0<=n1) ? n0*((n1+1)*n1)/2 : n1*((n0+1)*n0)/2];
    }
    else {
      cutoff_ = cutoff[(n0<=n1) ? n0*((n1+1)*n1)/2 : n1*((n0+1)*n0)/2];
    }
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


/*POT_LJ DESTRUCTOR: DEALLOCATES MEMORY*/
POT_LJ::~POT_LJ() {
  delete[] epsilon; epsilon=NULL;
  delete[] sigma;   sigma  =NULL;
  delete[] cutoff;  cutoff =NULL;
  delete[] Acutoff; Acutoff=NULL;
  return;
}
