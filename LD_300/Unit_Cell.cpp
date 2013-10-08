/*                           Unit_Cell.cpp                          */
/*                            12/02/2008                            */
/*********************************************************************
*    Source code file for functions in the UNIT_CELL class.          *
*********************************************************************/

/*DECLARE HEADERS*/
#include "LDCode.h"
#include "Unit_Cell.h"
#include <cmath>


/*DEFINE PREPROCESSOR VARIABLES*/


UNIT_CELL::UNIT_CELL() {
  coordinate = CARTESIAN;
  return;
}


UNIT_CELL::~UNIT_CELL() {
  delete[] X[0];
  X[0] = NULL;
  delete[] X;
  X = NULL;
  return;
}


string UNIT_CELL::Define(ifstream &Input) {
  //DECLARE LOCAL VARIABLES
  int i, atom_n;
  string str, temp;
  string Read_Next(ifstream &Input);
  void Coordinate(int coord, double *x, ifstream &Input);
  //void TransformVector(int coord, double *x);


  //PARSE THROUGH FILE
  natom = int(atof(Read_Next(Input).c_str())+0.0001);
  UC_DOF = DIM*natom;
  mat = new int[natom];
  X = new double*[natom];
  X[0] = new double[natom*DIM];
  for(i=0;i<natom;i++) {if(i<natom-1) {X[i+1] = X[i] + DIM;}}
  for(atom_n=0;atom_n<natom;atom_n++) {
    //identify keyword in LATTICE category
    if(Input.eof()) {Log <<"Error. Too few atoms in unit cell."<<endl;exit(0);}
    str = Read_Next(Input);
    temp = str;
    i = 0;while(temp[i]) {temp[i]=tolower(temp[i]);i++;}
    if(temp.compare(0, 4, "coordinate", 4)==0) {
      coordinate = KeywordToInt(Input);
      atom_n -= 1;
    }
    else {
      mat[atom_n] = Material->GetMaterial(str);
      Coordinate(coordinate, X[atom_n], Input);
    }
  }

  return(Read_Next(Input));
}


/*UNIT_CELL FUNCTION KeywordToInt: IDENTIFIES KEYWORDS USED FOR THE UNIT CELL*/
int UNIT_CELL::KeywordToInt(ifstream &Input) {
  int i=0;
  string str;
  string Read_Next(ifstream &Input);
  str = Read_Next(Input);
  while(str[i]) {str[i]=tolower(str[i]);i++;}
  if(str.compare(0, 3, "cartesian",  3)==0) {return(CARTESIAN);}//x, y, z
//  if(str.compare(0, 3, "cylindrical",3)==0) {return(CYLINDRICAL);}//r, theta, z
//  if(str.compare(0, 3, "spherical",  3)==0) {return(SPHERICAL);}//r, theta (from x axis), phi (from z axis)
  if(str.compare(0, 3, "direct",     3)==0) {return(DIRECT);}//a1, a2, a3
  if(str.compare(0, 3, "reciprocal", 3)==0) {return(RECIPROCAL);}//b1, b2, b3
  Log <<"Error. Unknown keyword for UNIT_CELL: "<<str<<endl;
  exit(0);
}


/*UNIT_CELL FUNCTION Output: WRITES DATA TO LOG FILE*/
void UNIT_CELL::Output() {
  int i, j;
  Log <<"\nUNIT_CELL = "<<natom<<endl;
  for(i=0;i<natom;i++) {
    Log <<"  "<<Material->symb[mat[i]]<<" ";
    for(j=0;j<DIM;j++) {Log <<" "<<X[i][j];}
    Log <<endl;
  }
  return;
}


/*UNIT_CELL FUNCTION Initialize: INITALIZES PARAMETERS AND NON-DIMENSIONALIZES*/
void UNIT_CELL::Initialize() {
  int i, j;
  for(i=0;i<natom;i++) {
    for(j=0;j<DIM;j++) {X[i][j] /= Parameter->length;}
  }
  return;
}
