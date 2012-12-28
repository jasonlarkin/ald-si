/*                           Material.cpp                           */
/*                            12/02/2008                            */
/*********************************************************************
*    Source code file for functions in the MATERIAL class.           *
*********************************************************************/

/*DECLARE HEADERS*/
#include "LDCode.h"
#include "Material.h"


/*DEFINE PREPROCESSOR VARIABLES*/


MATERIAL::MATERIAL() {return;}


MATERIAL::~MATERIAL() {
  delete[] symb;
  symb = NULL;
  delete[] mass;
  mass = NULL;
  return;
}


string MATERIAL::Define(ifstream &Input, POTENTIAL ***Pot) {
  //DECLARE LOCAL VARIABLES
  int i, j;
  string Read_Next(ifstream &Input);


  //DEALLOCATE Pot IF DEFINED
  if((*Pot)!=NULL) {
    Log  <<"Warning redefinition of MATERIAL all POTENTIAL data lost."<<endl;
    cout <<"Warning redefinition of MATERIAL all POTENTIAL data lost."<<endl;
    for(i=0;i<n_mat;i++) {delete Pot[0][i]; Pot[0][i]=NULL;}
    delete[] Pot[0]; Pot[0]=NULL;
  }


  //PARSE THROUGH FILE
  n_mat = int(atof(Read_Next(Input).c_str())+0.0001);
  symb = new string[n_mat];
  mass = new double[n_mat];
  //Read material symbols
  for(i=0;i<n_mat;i++) {
    if(Input.eof()) {Log <<"Error. Too few materials listed."<<endl;exit(0);}
    symb[i] = Read_Next(Input);
    mass[i] = atof(Read_Next(Input).c_str());
  }


  //ALLOCATE Pot
  (*Pot) = new POTENTIAL*[n_mat];
  for(i=0;i<n_mat;i++) {Pot[0][i] = NULL;}

  return(Read_Next(Input));
}


/*MATERIAL FUNCTION GetMaterial: RETURNS MATERIAL ID GIVEN THE SYMBOL*/
int MATERIAL::GetMaterial(string str) {
  for(int i=0;i<n_mat;i++) {if(str.compare(symb[i])==0) {return(i);}}
  Log <<"Unknown material "<<str<<endl;exit(0);
  return(-1);
}


/*MATERIAL FUNCTION Output: WRITES DATA TO LOG FILE*/
void MATERIAL::Output() {
  int i;
  Log <<"\nMATERIAL = "<<n_mat<<endl;
  for(i=0;i<n_mat;i++) {Log <<"  "<<symb[i]<<" = "<<mass[i]<<endl;}
  return;
}


/*MATERIAL FUNCTION Initialize: INITALIZES PARAMETERS AND NON-DIMENSIONALIZES*/
void MATERIAL::Initialize() {
  int i;
  for(i=0;i<n_mat;i++) {mass[i] /= Parameter->mass;}
  return;
}
