/*                          Read_Input.cpp                          */
/*                            12/02/2008                            */
/*********************************************************************
*    Subroutine that reads the data from the input file	             *
*  "LD_Input.txt".                                                   *
*********************************************************************/


/*DECLARE HEADERS*/
#include "LDCode.h"
#include "Variables.h"
#include "PreProcess.h"
#include "Dispersion.h"
#include "Harmonic.h"
#include "Anharmonic.h"
#include "PostProcess.h"
#include <cstdio>
#include <cstring>
#include <cmath>


/*DEFINE PREPROCESSOR COMMANDS*/


/*DECLARE GLOBAL VARIABLES*/
VARIABLES *Variables;         //VARIABLES must be declared before use, only needed for input


/*DEFINE SUBROUTINE Read_Input*/
void Read_Input(int argc, char **argv, POTENTIAL ***Potential,
                PREPROCESS *Preprocess, DISPERSION *Dispersion,
                HARMONIC *Harmonic, ANHARMONIC *Anharmonic,
                POSTPROCESS *Postprocess) {
  //DECLARE LOCAL VARIABLES
  int i, j, k;                //Counters
  string str;                 //Generic string
  ifstream Input;             //Input file
  string Read_Next(ifstream &Input);//Function to parse through input
  string PotDefine(ifstream &Input);//Function to define interatomic potential
  string GetPotential(ifstream &Input, POTENTIAL **Pot);


  //ALLOCATE DYNAMIC MEMORY
  Variables = new VARIABLES();//Linked list of variable names and values
  Material = new MATERIAL();  //Material names and masses
  Lattice = new LATTICE();    //Lattice definition
  Unit_Cell = new UNIT_CELL();//Unit cell definition
  Symmetry = new SYMMETRY();  //Symmetry operations


  //READ INPUT FILE
  char filename[1024];
  if(argc==1) {strcpy(filename, "LD_Input.txt");}
  else {strcpy(filename, argv[1]);}
  Input.open(filename);
  //Check for file open
  if(!Input.is_open()) {
    Log <<"Could not open input file '"<<filename<<"'."<<endl;exit(0);
  }
  else {
    Log  <<"\nParsing File '"<<filename<<"' For Input Data...\n"<<endl;
    cout <<"\nParsing File '"<<filename<<"' For Input Data...\n"<<endl;
  }
  //Parse through file
  str = Read_Next(Input);
  while(!Input.eof()) {
    //identify category keyword
    i = 0; while(str[i]) {str[i]=toupper(str[i]);i++;}
    if(str.compare(0, 3, "VARIABLES"  , 3)==0) {str=Variables->Define(Input);}        else
    if(str.compare(0, 3, "MATERIALS"  , 3)==0) {str=Material->Define(Input, Potential);} else
    if(str.compare(0, 3, "LATTICE"    , 3)==0) {str=Lattice->Define(Input);}          else
    if(str.compare(0, 3, "UNIT_CELL"  , 3)==0) {str=Unit_Cell->Define(Input);}        else
    if(str.compare(0, 3, "SYMMETRY"   , 3)==0) {str=Symmetry->ReadInput(Input);}      else
    if(str.compare(0, 3, "POTENTIAL"  , 3)==0) {str=GetPotential(Input, *Potential);} else
    if(str.compare(0, 3, "PREPROCESS" , 3)==0) {str=Preprocess->ReadInput(Input);}    else
    if(str.compare(0, 3, "DISPERSION" , 3)==0) {str=Dispersion->Define(Input);}       else
    if(str.compare(0, 3, "HARMONIC"   , 3)==0) {str=Harmonic->ReadInput(Input);}      else
    if(str.compare(0, 3, "SHIFT_WIDTH", 3)==0) {Log <<"SHIFT_WIDTH should be changed to ANHARMONIC"<<endl;str=Anharmonic->ReadInput(Input);} else
    if(str.compare(0, 3, "ANHARMONIC" , 3)==0) {str=Anharmonic->ReadInput(Input);}    else
    if(str.compare(0, 3, "TEMPERATURE", 3)==0) {Log <<"TEMPERATURE should be listed under ANHARMONIC"<<endl;Anharmonic->Temperature=atof(Read_Next(Input).c_str());str=Read_Next(Input);} else
    if(str.compare(0, 3, "POSTPROCESS", 3)==0) {str=Postprocess->ReadInput(Input);}   else
    {Log <<"Unexpected expression: "<<str<<endl;str = Read_Next(Input);}
  }
  Log  <<"\nDone reading input file.  Outputing stored data to log file for error checking..."<<endl;
  cout <<"\nDone reading input file.  Outputing stored data to log file for error checking..."<<endl;


  //OUTPUT DATA FOR ERROR CHECKING
  Log <<"\nVARIABLES"<<endl;
  delete Variables;  Variables = NULL;
  Lattice->Output();
  Material->Output();
  Unit_Cell->Output();
  Symmetry->Output();
  Preprocess->Output();
  Dispersion->Output();
  Harmonic->Output();
  Anharmonic->Output();
  Postprocess->Output();

  return;
}


/*SUBROUTINE Read_Next: USED TO READ FILE ONE WORD AT A TIME*/
string Read_Next(ifstream &Input) {
  string Math_Operation(string str);
  const int n = 7;
  int delim[n] = {9,10,13,32,44,61,124};/*{'\t','\n','?CR?',' ',',','=','|'}*/
  int i;
  char c;
  string str;
  str.clear();

  while (Input.get(c)) {
    //Check for comments
    if((int)c=='%') {
      Input.ignore(10001, '\n');
      if(!str.empty()) {return(Math_Operation(str));}
      else {continue;}
    }
    //Check for deliminators
    for(i=0;i<n;i++) {
      if((int)c==delim[i]) {
        if(str.empty()) {return(Read_Next(Input));}
        return(Math_Operation(str));
      }
    }
    str += c;
  }

  return(str);
}


/*SUBROUTINE Math_Operation: PERFORMS MATHEMATICAL OPERATIONS DEFINED IN A STRING*/
string Math_Operation(string str) {
  if( (str.find_first_of("$()^*/+-", 0))==string::npos ) {return(str);}
  //CHECK FOR VARIABLES
  string Op_Variable(string str);
  str = Op_Variable(str);
  //CHECK FOR PARENTHESES
  string Op_Parentheses(string str);
  str = Op_Parentheses(str);
  //CHECK FOR AND PERFORM REGULAR MATH OPERATIONS (^*/+-)
  string Op_Mathematics(string str);
  str = Op_Mathematics(str);

  return(str);
}


/*SUBROUTINE Op_Variable: SEARCHES FOR AND REPLACES VARIABLES IN str*/
string Op_Variable(string str) {
  //DECLARE LOCAL VARIABLES
  int i, j, n_char = str.size();
  //SEARCH FOR VARIABLES IN str
  for(i=0;i<n_char;i++) {
    if(int(str[i])=='$') {
      //Replace variable with value
      j=str.find_first_of("$()^*/+-", i+1);
      double value = Variables->GetValue(str.substr(i,j-i));
      string Op_Sign(string str, double &value, int &beg);
      str.replace(i, j-i, Op_Sign(str, value, i));
      i = 0;
      n_char = str.size();
    }
  }
  return(str);
}


/*SUBROUTINE Op_Parentheses: SEARCHES FOR EXPRESSIONS BOUNDED BY PARENTHESES*/
string Op_Parentheses(string str) {
  //DECLARE LOCAL VARIABLES
  int i, j, c, str_size = str.size();
  int n_open, n_close;
  //LOOK FOR INITIAL OPEN PARENTHESIS
  for(i=0;i<str_size;i++) {
    c = str[i];
    if(c==')') {Log <<"Bad parentheses: "<<str<<endl;exit(0);} else
    if(c=='(') {              //Initial open parenthesis
      j = i;
      n_open = 1;
      n_close = 0;
      while(n_open>n_close) { //Find matching close parenthesis
        if( (j+=1)>=str_size ) {Log <<"Bad parentheses: "<<str<<endl;exit(0);}
        c = str[j];
        if(c=='(') {n_open += 1;} else
        if(c==')') {n_close += 1;}
      }
      double value = strtod(Math_Operation(str.substr(i+1, j-i-1)).c_str(), NULL);
      string Op_Sign(string str, double &value, int &beg);
      str.replace(i, j-i+1, Op_Sign(str, value, i));
      str_size = str.size();
    }
  }
  return(str);
}


/*Op_Sign: DETERMINES IF THERE IS AN OPTIONAL SIGN IN FRONT OF AN EVALUATED VALUE*/
string Op_Sign(string str, double &value, int &beg) {
  if( (value<0.0)&&(beg!=0) ) {
    int op = int(str[beg-1]);
    if( (op=='+')||(op=='-') ) {
      beg-=1;
      if(beg!=0) {
        op = int(str[beg-1]);
        if( (op=='(')||(op=='*')||(op=='^')||(op=='*')||(op=='/')||(op=='+')||(op=='-') ) {
          if(int(str[beg])=='-') {value = -value;}
        }
      }
      else {if(int(str[beg])=='-') {value = -value;}}
    }
  }
  char c_value[101];
  sprintf(c_value, "%g", value);
  return(string(c_value));
}


/*SUBROUTINE Op_Mathematics: SEARCHES FOR MATHEMATICAL OPERATIONS*/
string Op_Mathematics(string str) {
  MATH_EXPRESSION num_op(str.c_str());//Constructor tokenizes expression
  if(num_op.op==0) {return(str);}
  num_op.Evaluate(int('^'), -1);
  num_op.Evaluate(int('*'), int('/'));
  num_op.Evaluate(int('+'), int('-'));
  char ch[101];
  sprintf(ch, "%g", num_op.number);
  str = ch;
  return(str);
}


/*SUBROUTINE Coordinate: READS AND TRANSFORMS VECTORS TO CARTESIAN COORDINATES*/
void Coordinate(int coord, double *x, ifstream &Input) {
  int i, j;
  double **T=NULL, value;
  string Read_Next(ifstream &Input);
  switch(coord) {
   case 0:
    for(i=0;i<DIM;i++) {x[i] = atof(Read_Next(Input).c_str());}
   return;
   case 3:
    T = Lattice->a;
   break;
   case 4:
    T = Lattice->b;
   break;
  }
  for(i=0;i<DIM;i++) {x[i] = 0.0;}
  for(i=0;i<DIM;i++) {
    value = atof(Read_Next(Input).c_str());
    for(j=0;j<DIM;j++) {x[j] += value*T[i][j];}
  }
  return;
}
