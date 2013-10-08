/*                            Lattice.cpp                           */
/*                            12/02/2008                            */
/*********************************************************************
*    Subroutines for the LATTICE class.                              *
*********************************************************************/

/*DECLARE HEADERS*/
#include "LDCode.h"
#include "Lattice.h"
#include <cmath>


/*DEFINE PREPROCESSOR VARIABLES*/


/*DECLARE GLOBAL VARIABLES*/


/*LATTICE CLASS CONSTRUCTOR*/
LATTICE::LATTICE() {
  int i, j;
  a = new double*[3];
  a[0] = new double[3*3];
  b = new double*[3];
  b[0] = new double[3*3];
  for(i=0;i<3;i++) {
    if(i<3-1) {
      a[i+1] = a[i] + 3;
      b[i+1] = b[i] + 3;
    }
    for(j=0;j<3;j++) {a[i][j] = b[i][j] = 0.0;}
    a[i][i] = 1.0;
    b[i][i] = 1.0;
  }
  BC = new int[3];
  N = new int[3];
  lower = new int[3];
  upper = new int[3];
  for(i=0;i<3;i++) {
    BC[i] = FREE;
    N[i] = 1;
    lower[i] = upper[i] = 0;
  }
  return;
}


/*LATTICE CLASS DESTRUCTOR*/
LATTICE::~LATTICE() {
  delete[] a[0];  a[0] = NULL;
  delete[] a;  a = NULL;
  delete[] b[0];  b[0] = NULL;
  delete[] b;  b = NULL;
  delete[] BC;  BC = NULL;
  delete[] N;  N = NULL;
  delete[] lower;  lower = NULL;
  delete[] upper;  upper = NULL;
  return;
}


/*LATTICE CLASS SUBROUTINE Define: USED TO STORE DATA*/
string LATTICE::Define(ifstream &Input) {
  //DECLARE LOCAL VARIABLES
  int i;
  string str;
  string Read_Next(ifstream &Input);


  //PARSE THROUGH FILE
  str = Read_Next(Input);
  i = int(atof(str.c_str())+0.0001);
  if(i!=DIM) {Log <<"Dimension set to constant value of '"<<DIM<<"' in program."<<endl;}
  while(!Input.eof()) {
    //identify keyword in LATTICE category
    str = Read_Next(Input);
    i = 0;
    while(str[i]) {str[i]=tolower(str[i]);i++;}
    if     (str.compare(0, 2, "a1",  2)==0) {for(i=0;i<DIM;i++) {a[0][i] = atof(Read_Next(Input).c_str());}}
    else if(str.compare(0, 2, "a2",  2)==0) {for(i=0;i<DIM;i++) {a[1][i] = atof(Read_Next(Input).c_str());}}
    else if(str.compare(0, 2, "a3",  2)==0) {for(i=0;i<DIM;i++) {a[2][i] = atof(Read_Next(Input).c_str());}}
    else if(str.compare(0, 3, "bc1", 3)==0) {BC[0] = KeywordToInt(Input);}
    else if(str.compare(0, 3, "bc2", 3)==0) {BC[1] = KeywordToInt(Input);}
    else if(str.compare(0, 3, "bc3", 3)==0) {BC[2] = KeywordToInt(Input);}
    else if(str.compare(0, 2, "n1",  2)==0) {N[0] = int(atof(Read_Next(Input).c_str())+0.0001);}
    else if(str.compare(0, 2, "n2",  2)==0) {N[1] = int(atof(Read_Next(Input).c_str())+0.0001);}
    else if(str.compare(0, 2, "n3",  2)==0) {N[2] = int(atof(Read_Next(Input).c_str())+0.0001);}
    else{break;}
  }


  //CHECK CONSISTENCY
  for(i=0;i<DIM;i++) {
    if((BC[i]==2)&&(N[i]>1)) {Log <<"Warning Inconsistent BC and N for direction "<<i+1<<endl;}
  }


  //COMPUTE VOLUME OF UNIT CELL, RECIPROCAL LATTICE VECTORS, AND BOUNDS OF WAVE VECTOR
  int j;
  int LU(int, double**, int*);
  void LU_Inverse(int, double**, int*);
  //Volume
  double LU_Determinant(int n, double **A, double sign);//Computes determinant of A
  double **temp = new double*[DIM];
  temp[0] = new double[DIM*DIM];
  for(i=0;i<DIM;i++) {
    if(i<DIM-1) {temp[i+1] = temp[i] + DIM;}
    for(j=0;j<DIM;j++) {temp[i][j] = Lattice->a[i][j];}
  }
  i = LU(DIM, temp, NULL);
  V = fabs(LU_Determinant(DIM, temp, double(i)));//Determinant gives the hyper-volume
  delete[] temp[0];  temp[0] = NULL;
  delete[] temp;  temp = NULL;
  //Reciprocal lattice vectors
  int *permute = new int[DIM];
  for(i=0;i<DIM;i++) {
    lower[i] = (1-N[i])/2;
		upper[i] = N[i]/2;
    for(j=0;j<DIM;j++) {b[i][j] =  a[j][i];}
  }
  LU(DIM, b, permute);
  LU_Inverse(DIM, b, permute);
  delete[] permute;  permute = NULL;


  return(str);
}


/*LATTICE FUNCTION KeywordToInt: IDENTIFIES KEYWORDS USED FOR THE LATTICE*/
int LATTICE::KeywordToInt(ifstream &Input) {
  int i=0;
  string str;
  string Read_Next(ifstream &Input);
  str = Read_Next(Input);
  while(str[i]) {str[i]=tolower(str[i]);i++;}
  if(str.compare(0, 3, "periodic",   3)==0) {return(PERIODIC);}
  if(str.compare(0, 3, "scattering", 3)==0) {return(SCATTERING);}
  if(str.compare(0, 3, "free",       3)==0) {return(FREE);}
  Log <<"Error. Unknown keyword for LATTICE: "<<str<<endl;
  exit(0);
}


/*LATTICE FUNCTION Output: WRITES DATA TO LOG FILE*/
void LATTICE::Output() {
  int i, j;
  Log <<"\nLATTICE = "<<DIM<<endl;
  for(i=0;i<DIM;i++) {
    Log <<"  a"<<i+1<<" =";
    for(j=0;j<DIM;j++) {Log <<" "<<a[i][j];}
    Log <<endl;
  }
  for(i=0;i<DIM;i++) {Log <<"  BC"<<i+1<<" = "<<BC[i]<<endl;}
  for(i=0;i<DIM;i++) {Log <<"  N"<<i+1<<"  = "<<N[i] <<endl;}
  Log <<"  Volume = "<<V<<endl;
  for(i=0;i<DIM;i++) {
    Log <<"  b"<<i+1<<" =";
    for(j=0;j<DIM;j++) {
      Log <<" "<<b[i][j];
    }
    Log <<endl;
  }
  return;
}


/*LATTICE FUNCTION Initialize: INITALIZES PARAMETERS AND NON-DIMENSIONALIZES*/
void LATTICE::Initialize() {
  int i, j;
  for(i=0;i<DIM;i++) {
    for(j=0;j<DIM;j++) {
      a[i][j] /= Parameter->length;
      b[i][j] *= Parameter->length;
    }
    V /= Parameter->length;
  }
  return;
}
