/*                           Variables.cpp                          */
/*                            12/02/2008                            */
/*********************************************************************
*    Source file for functions in the VARIABLE and MATH_EXPRESSION   *
*  classes.                                                          *
*********************************************************************/

/*DECLARE HEADERS*/
#include "LDCode.h"
#include "Variables.h"
#include <cmath>


/*DEFINE PREPROCESSOR VARIABLES*/


/*VARIABLES CONSTRUCTOR: INITIALIZES VALUES*/
VARIABLES::VARIABLES() {
  symb.clear();
  value = 0.0;
  next = this;
  return;
}


/*VARIABLES FUNCTION Define: READS AND STORES VARIABLES*/
string VARIABLES::Define(ifstream &Input) {
  //DECLARE LOCAL VARIABLES
  string Read_Next(ifstream &Input);
  string str;
  int i;
  char c;
  const int n = 7;
  int delim[n] = {9,10,13,32,44,61,124};/*{'\t','\n','?CR?',' ',',','=','|'}*/


  //GET VARIABLE NAME (WITHOUT LETTING Read_Next TRY TO REPLACE)
  str.clear();
  while (Input.get(c)) {
    //Check for comments
    if(int(c)=='%') {
      Input.ignore(10001, '\n');
      continue;
    }
    //Check for deliminators
    for(i=0;i<n;i++) {
      if(int(c)==delim[i]) {
        i=-1;
        break;
      }
    }
    if(i!=-1) {str += c;break;}
  }
  if(int(c)!='$') {
    //Not a valid variable name (return)
    Input.putback(c);
    return(Read_Next(Input));
  }
  str = Read_Next(Input);//NOTE: this allows deliminators between '$' and varaible name
  str.insert(0, "$");


  //FIND VARIABLE TO REDEFINE IT OR FIND UNUSED SLOT
  VARIABLES *var = this;
  while(var!=var->next) {
    if(str.compare(var->symb)==0) {break;}
    var = var->next;
  }


  //STORE SYMBOL AND VALUE
  var->symb = str;
  str = Read_Next(Input);
  i = int(str[0]);
  if( (i<43)||(i>57)||(i==44)||(i==47) ) {//Check if a number is present
    Log <<str<<" is not a valid expression for variable "<<var->symb<<endl;
    exit(0);
  }
  var->value = atof(str.c_str());
  if(var==var->next) {var->next = new VARIABLES();}//Allocate memory for next variable

  return(Define(Input));
}


double VARIABLES::GetValue(string str) {
  if(str.compare(symb)==0) {return(value);}
  else {
    if(next==this) {Log <<"Undefined variable "<<str<<endl;exit(0);}
    return(next->GetValue(str));
  }
}


/*VARIABLES DESTRUCTOR: PRINTS VARIABLES AND DEALLOCATES MEMORY*/
VARIABLES::~VARIABLES() {
  if(next!=this) {
    Log <<"  "<<symb<<" = "<<value<<endl;
    symb.clear();
    delete next;  next=NULL;
  }
  return;
}


/*MATH_EXPRESSION CONSTRUCTOR: CREATES LINKED LIST OF NUMBERS AND OPERATORS*/
MATH_EXPRESSION::MATH_EXPRESSION(const char *c_str) {
  //DECLARE LOCAL VARIABLES AND INITIALIZE
  char *str_next;
  number=double(op=0);
  next=NULL;
  //FIND NEXT NUMBER AND OPERATOR (IF PRESENT)
  if((op=int(c_str[0]))==0) {Log <<"Expression ends with operator."<<endl;exit(0);}
  number = strtod(c_str, &str_next);
  if(c_str==str_next) {Log <<"Could not interpert expression: "<<c_str<<endl;exit(0);}
  if((op=int(str_next[0]))!=0) {
    if( (op!='^')&&(op!='*')&&(op!='/')&&(op!='+')&&(op!='-') ) {
      Log <<"Could not interpert expression: "<<c_str<<endl;exit(0);
    }
    next = new MATH_EXPRESSION(str_next+1);
  }
  return;
}


/*MATH_EXPRESSION FUNCTION Evaluate: EVALUATES EVERY OCCURANCE OF THE OPERATION*/
void MATH_EXPRESSION::Evaluate(int op1, int op2) {
  if( (op1==op)||(op2==op) ) {
    //Evaluate mathematical operation
    if(op=='^') {number = pow(number, next->number);} else
    if(op=='*') {number = number * next->number;} else
    if(op=='/') {number = number / next->number;} else
    if(op=='+') {number = number + next->number;} else
    if(op=='-') {number = number - next->number;} else
    {Log <<"Operation: "<<op<<" is undefined."<<endl;}
    //Reduce linked list by copying numbers and operators down one
    op = next->op;
    if(next->next==NULL) {
      delete next;
      next=NULL;
      return;
    }
    next->Copy();
    Evaluate(op1, op2);
  }
  if(next!=NULL) {next->Evaluate(op1, op2);}
  return;


  /*if(_op==op) {
    //Evaluate mathematical operation
    if(op=='^') {number = pow(number, next->number);} else
    if(op=='*') {number = number * next->number;} else
    if(op=='/') {number = number / next->number;} else
    if(op=='+') {number = number + next->number;} else
    if(op=='-') {number = number - next->number;} else
    {Log <<"Operation: "<<op<<" is undefined."<<endl;}
    //Reduce linked list by copying numbers and operators down one
    op = next->op;
    if(next->next==NULL) {
      delete next;
      next=NULL;
      return;
    }
    next->Copy();
    Evaluate(_op);
  }
  if(next!=NULL) {next->Evaluate(_op);}
  return;*/
}


/*MATH_EXPRESSION FUNCTION Copy: COPIES NEXT ELEMENT IN LINKED LIST TO PREVIOUS ELEMENT AND DELETES EMPTY ELEMENTS*/
void MATH_EXPRESSION::Copy(void) {
  if(next==NULL) {return;}
  number = next->number;
  op = next->op;
  if(next->next==NULL) {delete next;next=NULL;}
  else {next->Copy();}
  return;
}


/*MATH_EXPRESSION DESTRUCTOR: DEALLOCATES HIGHER ELEMENTS IN LIST*/
MATH_EXPRESSION::~MATH_EXPRESSION() {
  if(next!=NULL) {delete next;  next=NULL;}
  return;
}
