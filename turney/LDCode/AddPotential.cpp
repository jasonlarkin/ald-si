/*                         AddPotential.cpp                         */
/*                            12/02/2008                            */
/*********************************************************************
*    This file contains all of the subroutines that must be altered  *
*  when adding additional potentials.  Add new potentials as classes *
*  (see pre-defined potentials for examples).                        *
*********************************************************************/

/*DECLARE HEADERS*/
#include "LDCode.h"
#include <cstring>


/*DEFINE PREPROCESSOR COMMANDS*/


/**THIS FUNCTION CONTAINS ALL OPERATIONS THAT NEED ALTERED WHEN ADDING POTENTIALS**/
/*SUBROUTINE PotentialKeywords: ALLOCATES MEMORY FOR POTENTIALS*/
POTENTIAL *PotentialKeywords(char *symb, int n_mat, int *mat) {
  //CONVERT symb TO UPPER CASE
  int i=0;
  while(symb[i]) {symb[i]=toupper(symb[i]);i++;}

  //COMPARE POTENTIAL TO THOSE DEFINED AND ALLOCATE MEMORY OF CORRECT TYPE
  /*****Add new potential keywords and class types here (upper case)*****/
  if(strcmp(symb, "LJ"       )==0) {return(new POT_LJ       (n_mat, mat));}
  if(strcmp(symb, "LJ_SPLINE")==0) {return(new POT_LJ_SPLINE(n_mat, mat));}
  if(strcmp(symb, "INV12"    )==0) {return(new POT_INV12    (n_mat, mat));}
  if(strcmp(symb, "BUCK"     )==0) {return(new POT_BUCK     (n_mat, mat));}
  if(strcmp(symb, "SW2"      )==0) {return(new POT_SW2      (n_mat, mat));}
  if(strcmp(symb, "SW3"      )==0) {return(new POT_SW3      (n_mat, mat));}
  if(strcmp(symb, "REBO2"    )==0) {return(new POT_REBO2    (n_mat, mat));}
  if(strcmp(symb, "REBO3"    )==0) {return(new POT_REBO3    (n_mat, mat));}
  if(strcmp(symb, "TERSOFF"  )==0) {return(new POT_TERSOFF  (n_mat, mat));}
  if(strcmp(symb, "COULOMB"  )==0) {return(new POT_COULOMB  (n_mat, mat));}
  Log <<"Error. Unknown potential: "<<symb<<"\nValid potentials (identifiers) include:"<<endl;
  Log <<"12-6 Lennard-Jones (LJ)"<<endl;
  Log <<"12-6 Lennard-Jones with quintic spline (LJ_SPLINE)"<<endl;
  Log <<"Inverse 12th order (INV12)"<<endl;
  Log <<"Buckingham (BUCK)"<<endl;
  Log <<"Stillinger-Weber, 2-body part (SW2)"<<endl;
  Log <<"Stillinger-Weber, 3-body part (SW3)"<<endl;
  Log <<"Brenner REBO, 2-body part (REBO2)"<<endl;
  Log <<"Brenner REBO, 3-body part (REBO3)"<<endl;
  Log <<"Tersoff (TERSOFF)"<<endl;
  Log <<"Coulomb (COULOMB)"<<endl;
  exit(0);
  return(NULL);
}


/** THESE FUNCTIONS SHOULD NOT NEED ALTERED WHEN ADDING POTENTIALS **/
/*SUBROUTINE GetPotential: IDENTIFIES THE POTENTIAL TO USE*/
string GetPotential(ifstream &Input, POTENTIAL **Pot) {
  //DECLARE LOCAL VARIABLES
  int i, n_mat, *mat;
  char *cstr;
  POTENTIAL *P_new, *P_ptr;
  string str;
  string Read_Next(ifstream &Input);
  //READ INITAL DATA FOR DEFINING POTENTIAL TYPE
  str = Read_Next(Input);
  n_mat = int(atof(Read_Next(Input).c_str())+0.001);
  mat = new int[n_mat];
  for(i=0;i<n_mat;i++) {mat[i] = Material->GetMaterial(Read_Next(Input));}
  //COMPARE POTENTIAL TO THOSE DEFINED
  cstr = new char[str.size()+1];
  strcpy(cstr, str.c_str());
  P_new = PotentialKeywords(cstr, n_mat, mat);
  //ASSIGN LOCATION OF MEMORY AND READ DATA
  P_ptr = Pot[mat[0]];
  if(P_ptr==NULL) {Pot[mat[0]] = P_new;}
  else {
    while(P_ptr->next!=NULL) {P_ptr = P_ptr->next;}
    P_ptr->next = P_new;
  }
  str = P_new->ReadInput(Input);
  delete[] mat;  mat=NULL;
  delete[] cstr; cstr=NULL;

  return(str);
}


/*SUBROUTINE InitializePotential: COPIES POTENTIALS SO THAT EACH ATOM HAS
  ALL NEEDED POTENTIALS ASSOCIATED WITH IT AND NON-DIMENSIONALIZES*/
void InitializePotential(POTENTIAL **Pot) {
  //DECLARE LOCAL VARIABLES
  int i, j;
  POTENTIAL *P_ptr, *P_new;
  //INITIALIZE (NON-DIMENSIONALIZE) POTENTIAL PARAMETERS
  for(i=0;i<Material->n_mat;i++) {
    P_ptr = Pot[i];
    while(P_ptr!=NULL) {
      P_ptr->Print();         //Print potential parameters
      P_ptr->Initialize();    //Non-dimensionalizes
      P_ptr = P_ptr->next;
    }
  }
  //COPY POTENTIALS SO Pot[i] HAS ALL POTENTIALS NEEDED FOR ATOM i
  for(i=0;i<Material->n_mat;i++) {
    P_ptr = Pot[i];
    while(P_ptr!=NULL) {
      if(i!=P_ptr->mat[0]) {P_ptr=P_ptr->next;continue;}//Skip if already considered
      for(j=1;j<P_ptr->n_mat;j++) {
        P_new = Pot[P_ptr->mat[j]];
        if(P_new==NULL) {Pot[P_ptr->mat[j]] = P_ptr->Copy();}
        else {
          while(P_new->next!=NULL) {P_new=P_new->next;}
          P_new->next = P_ptr->Copy();
        }
      }
      P_ptr = P_ptr->next;
    }
  }
  return;
}
