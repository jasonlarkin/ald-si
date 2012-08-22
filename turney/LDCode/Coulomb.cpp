/*						                Coulomb.cpp		            			  	  */
/*							              12/23/2009		              					*/
/*********************************************************************
*    File that contains the POT_COULOMB class function definitions.	 *
*  POT_COULOMB is a derivative of the POTENTIAL class.		        	 *
*********************************************************************/

/*DEFINE HEADERS*/
#include "Coulomb.h"
#include <cmath>


/*DEFINE SUBROUTINES*/
string Read_Next(ifstream &Input);


/*POT_COULOMB CONSTRUCTOR: STORES POTENTIAL AND POTENTIAL PARAMETER NAMES*/
POT_COULOMB::POT_COULOMB(int n_mat_, int *mat_) : POTENTIAL(n_mat_, mat_) {
  int i;
  n_expect = n_mat;
  Pot_Symb = new char[8];
  strcpy(Pot_Symb, "COULOMB");
  mat_pos = new int[2];
  q       = new double[n_expect];
  for(i=0;i<n_expect;i++) {q[i]=0.0;}
  //Parameters for Ewald summations
  cutoff  = new double[1];
  Acutoff = new double[1];
  *cutoff  = 0.0;
  *Acutoff = 0.0;
  alpha    = 0.0;
  erfccut = a_pi = E_const = F_const = 0.0;
  return;
}


/*POT_COULOMB ReadInput: READS AND STORES POTENTIAL PARAMETER VALUES*/
string POT_COULOMB::ReadInput(ifstream &Input) {
  string str;
  while(!Input.eof()) {
    str = Read_Next(Input);
    if(str.compare(      "q"         )==0) {Fill(q, Input);} else
    if(str.compare(0, 3, "alpha"  , 3)==0) {alpha   =atof(Read_Next(Input).c_str());} else
    if(str.compare(0, 3, "cutoff" , 3)==0) {*cutoff =atof(Read_Next(Input).c_str());} else
    if(str.compare(0, 3, "Acutoff", 3)==0) {*Acutoff=atof(Read_Next(Input).c_str());} else
    {break;}
  }
  return(str);
}


/*POT_COULOMB FUNCTION GetScale: INITIALIZES INITIALIZES THE LENGTH AND ENERGY SCALES*/
void POT_COULOMB::GetScale(double &e, double &l) {
  e = 1.60217646e-19; // J/eV (C/electron)
  l = 5.0e-10;        // 5 angstroms
	return;
}


/*POT_COULOMB FUNCTION Initialize: INITIALIZES THE POTENTIAL PARAMETERS*/
void POT_COULOMB::Initialize() {
  int i;
	for(i=0;i<n_expect;i++) {
	  q[i] *= 1.0/sqrt(4.0*4.0*atan(1.0)*8.85418782e-12*Parameter->energy*Parameter->length);//q/sqrt(4*pi*epsilon_0)
	}
	if(alpha<=0.0) {alpha = 2.4/(*cutoff);}//Set alpha based on figs 7 and 8 in Wolf paper
  if((*Acutoff)<=0.0) {*Acutoff = (*cutoff);}
  *cutoff  /= Parameter->length;
  *Acutoff /= Parameter->length;
  alpha    *= Parameter->length;
  double a_cut = alpha*(*cutoff);
  erfccut = erfc(a_cut)/(*cutoff);
  a_pi = 2.0*alpha/sqrt(Parameter->pi);
  E_const = 0.5*(erfccut+a_pi);
  F_const = (erfccut+a_pi*exp(-a_cut*a_cut))/((*cutoff)*(*cutoff));
	return;
}


/*POT_COULOMB FUNCTION Print: PRINTS POTENTIAL PARAMETERS*/
void POT_COULOMB::Print() {
  int i;
  Log <<"\nPOTENTIAL "<<Pot_Symb<<" "<<n_mat;
  for(i=0;i<n_mat;i++) {Log <<" "<<Material->symb[mat[i]];}
  Log <<endl<<"  q       =";
  for(i=0;i<n_expect;i++) {Log <<" "<<q[i];}
  Log <<" C^2"<<endl;
  Log <<"  alpha   = "<<alpha<<" 1/m"<<endl;
  Log <<"  cutoff  = "<<(*cutoff)<<" m"<<endl;
  Log <<"  Acutoff = "<<(*Acutoff)<<" m"<<endl;
  Log <<endl;
  return;
}


/*POT_COULOMB FUNCTION Copy: RETURNS A POINTER TO A NEW POTENTIAL WITH A DUPLICATE OF THE DATA*/
POTENTIAL *POT_COULOMB::Copy() {
  int i;
  POT_COULOMB *P_new = new POT_COULOMB(n_mat, mat);
  for(i=0;i<n_expect;i++) {P_new->q[i] = q[i];}
  P_new->cutoff[0]  = (*cutoff);
  P_new->Acutoff[0] = (*Acutoff);
  P_new->alpha      = alpha;
  P_new->E_const    = E_const;
  P_new->F_const    = F_const;
  P_new->erfccut    = erfccut;
  P_new->a_pi       = a_pi;
  return(P_new);
}


/*POT_COULOMB FUNCTION Energy: COMPUTES THE POTENTIAL ENERGY*/
double POT_COULOMB::Energy(double **X) {
  int i;
  double phi;
	double r_mag = 0.0;
  for(i=0;i<DIM;i++) {
    phi = X[0][i] - X[1][i];
    r_mag += phi*phi;
  }
  if(r_mag<1.0e-10) {phi = -2.0*E_const;}//Correct for self-interactions
  else {                      //Sum over direct space
    r_mag = sqrt(r_mag);
    phi = erfc(alpha*r_mag)/r_mag - erfccut;
  }
  phi *= q[mat_pos[0]]*q[mat_pos[1]];

	return(phi);
}


/*POT_COULOMB::ForceConsants: COMPUTES THE HARMONIC FORCE CONSTANTS*/
void POT_COULOMB::ForceConstants(int b0, POT_DER_X *F) {
  //DECLARE LOCAL VARIABLES
  int n0, n1;                 //Material's position in array mat[]
  int i;                      //Generic counters
  int b1;                     //Atom counter for 1st neighbor
  int n_neig;                 //Number of neighbors
  int *B, **L;                //List of atomic and unit cell indicies
  double cutoff_;             //Cutoff radius
  double **X;                 //List of positions
  N_LIST *N_list;             //Linked list of neighbors
  N_LIST *N_list_ptr;         //Pointer to entry in linked list
  N_LIST *BuildNeighborList_R(int &, int, double);//Computes the neighbor list

cout <<b0<<endl;
  //INITIALIZE VARIABLES
  B = new int[2];
  L = new int*[2];
  int *L1 = new int[DIM];
  X = new double*[2];
  X[0] = new double[2*DIM];
  X[1] = X[0] + DIM;
  for(n0=0;n0<n_mat;n0++) {if(mat[n0]==Unit_Cell->mat[b0]) {break;}}
  mat_pos[0] = n0;
  for(i=0;i<DIM;i++) {X[0][i] = Unit_Cell->X[b0][i];}
  B[0] = b0;
  L[0] = NULL;


  //LOOP OVER ALL MATERIAL INTERACTIONS
  //Correct for self-interaction
  if(typeid(F[0])==typeid(PD_ENERGY_FORCE)) {
    for(i=0;i<DIM;i++) {
      X[1][i] = X[0][i];
      L1[i] = 0;
    }
    B[1] = b0;
    L[1] = L1;
    mat_pos[1] = n0;
    F->Compute(0, 2, B, L, X, this);
  }
  //Do all other interactions
  for(n1=0;n1<n_mat;n1++) {
    //GET CUTOFF
    mat_pos[1] = n1;
    if(typeid(F[0])==typeid(PD_ANHARMONIC)) {cutoff_=(*Acutoff);}
    else {cutoff_=(*cutoff);}
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
  delete[] B;    B    = NULL;
  delete[] L1;   L1   = NULL;
  delete[] L;    L    = NULL;
  delete[] X[0]; X[0] = NULL;
  delete[] X;    X    = NULL;

  return;
}


/*POT_BUCK DESTRUCTOR: DEALLOCATES MEMORY*/
POT_COULOMB::~POT_COULOMB() {
  delete[] q;       q      =NULL;
  delete[] cutoff;  cutoff =NULL;
  delete[] Acutoff; Acutoff=NULL;
  return;
}
