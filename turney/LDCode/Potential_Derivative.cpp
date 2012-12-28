/*                     Potential_Derivative.cpp                     */
/*                            06/08/2009                            */
/*********************************************************************
*    Subroutines for the POT_DER classes.                            *
*********************************************************************/

/*DECLARE HEADERS*/
#include "LDCode.h"
#include "Potential_Derivative.h"
#include <cmath>


/*DEFINE PREPROCESSOR COMMANDS*/


/********************* Class POT_DER_LINK_LIST **********************/
/*POT_DER_LINK_LIST Constructor: Allocates memory*/
POT_DER_LINK_LIST::POT_DER_LINK_LIST(int O_der) {
  next = NULL;
  int i, j;
  int n_force = 1;
  for(i=0;i<O_der;i++) {n_force *= DIM;}
  fc = new double[n_force];
  for(i=0;i<n_force;i++) {fc[i] = 0.0;}
  return;
}


/*POT_DER_LINK_LIST Destructor: Deallocates memory*/
POT_DER_LINK_LIST::~POT_DER_LINK_LIST() {
  if(next!=NULL) {delete next;  next=NULL;}
  delete[] fc;  fc=NULL;
  return;
}


/*POT_DER_LINK_LIST Compute: Computes the force constants, returns true if all are zero*/
void POT_DER_LINK_LIST::Compute(int O_der, double **X, double **DX, POTENTIAL *pot) {
  //DECLARE LOCAL VARIABLES
  int i, j;                   //Counters
  int n_force;                //Total number of force constants
  double **dx;                //New derivative pointer to work with
  int *dx_counter;            //Tracks which dimension we are on for each atom


  //ALLOCATE MEMORY AND INITIALIZE
  dx_counter = new int[O_der];
  dx = new double*[O_der];
  for(j=0;j<O_der;j++) {
    dx[j] = DX[j];
    dx_counter[j] = 0;
  }
  n_force = 1;
  for(i=0;i<O_der;i++) {n_force *= DIM;}


  //COMPUTE FORCE CONSTANTS TO Order^th ORDER
  for(i=0;i<n_force;i++) {
    //Compute and store energy derivative
    fc[i] += pot->Num_Der(O_der, X, dx);
    //Manage pointers:ordering= x...xx, x...xy, x...xz, x...yx, ..., z...zz
    for(j=O_der-1;j>=0;j--) {// 0...00, 0...01, 0...02, 0...10, ..., 2...22
      if((dx_counter[j]+=1)==DIM) {//If counter resets continue to next atom
        dx_counter[j] = 0;
        dx[j] = DX[j];
      }
      else {                  //If counter does not reset, increment and break loop
        dx[j] += 1;
        break;
      }
    }
  }


  //DEALLOCATE MEMORY AND RETURN
  delete[] dx;  dx=NULL;
  delete[] dx_counter;  dx_counter=NULL;
  return;
}


/*POT_DER_LINK_LIST Zero: Zeros small force constants returns true if all zero*/
bool POT_DER_LINK_LIST::Scale(int O_der, double factor) {
  //DECLARE LOCAL VARIABLES
  int i, j=1;                 //Counters
  bool flag = true;           //true if all force constants are zero
  for(i=0;i<O_der;i++) {j *= DIM;}
  for(i=0;i<j;i++) {
    fc[i] *= factor;
    if(fabs(fc[i])<1.0e-13) {fc[i] = 0.0;}//ZERO SMALL VALUES
    else {flag = false;}
  }
  return(flag);
}


/*POT_DER_LINK_LIST Print: OUTPUTS FORCE CONSTANTS*/
void POT_DER_LINK_LIST::Print(int O_der) {
  int i, j = 1, k = 0;
  for(i=0;i<O_der;i++) {j*= DIM;}
  for(i=0;i<j;i++) {
    Log  <<" "<<fc[i];
    cout <<" "<<fc[i];
  }
  Log  <<endl;
  cout <<endl;
  return;
}


/********************* Class POT_DER_IDENTIFIER *********************/
/*POT_DER_IDENTIFIER Constructor: ALLOCATES MEMORY AND ASSIGNS DATA IF GIVEN*/
POT_DER_IDENTIFIER::POT_DER_IDENTIFIER(int N_atom, int *B=NULL, int **L=NULL) {
  next = NULL;
  N_atom -= 1;
  if(N_atom<=0) {
    b = NULL;
    l = NULL;
  }
  else {
    int i, j;
    b = new int[N_atom];
    l = new int*[N_atom];
    l[0] = new int[N_atom*DIM];
    for(i=0;i<N_atom;i++) {
      if(B!=NULL) {b[i] = B[i];}
      if(i+1<N_atom) {l[i+1] = l[i] + DIM;}
      if(L!=NULL) {for(j=0;j<DIM;j++) {l[i][j] = L[i][j];}}
    }
  }
  return;
}


/*POT_DER_IDENTIFIER Constructor: ALLOCATES MEMORY AND ASSIGNS DATA*/
POT_DER_IDENTIFIER::~POT_DER_IDENTIFIER() {
  delete next; next=NULL;
  if(l!=NULL) {delete[] l[0]; l[0]=NULL;}
  delete[] l; l=NULL;
  delete[] b; b=NULL;
  return;
}


/*POT_DER_IDENTIFIER Position: DETERMINIES THE POSITION TO WRITE DATA
function returns a position in the list to allocate memory (alloc==true) and/or write data
N_atom should be the number of uniqe atoms in the derivative excluding b0*/
int POT_DER_IDENTIFIER::Position(int N_atom, int *B, int **L, bool &alloc) {
  int i, j;                   //Counters
  bool eq = true;             //True if B and L equal b and l
  for(i=0;i<N_atom;i++) {
    if(B[i]>b[i]) {alloc=true; return(0);}
    else if(B[i]==b[i]) {
      for(j=0;j<DIM;j++) {
        if(L[i][j]>l[i][j]) {alloc=true; return(0);}
        else if(L[i][j]!=l[i][j]) {eq=false; break;}
      }
    }
    else {eq = false;}
    if(!eq) {break;}
  }
  if(eq) {alloc=false; return(0);}
  int pos = 1;
  if(next==NULL) {alloc=true;}
  else {pos += next->Position(N_atom, B, L, alloc);}
  if( (pos==1)&&(alloc) ) {   //Allocate memory and store b and l
    POT_DER_IDENTIFIER *ptr = next;
    next = new POT_DER_IDENTIFIER(N_atom+1, B, L);
    next->next = ptr;
  }
  return(pos);
}

/************************* Class POT_DER_X **************************/
/*POT_DER_X Constructor: Does nothing*/
POT_DER_X::POT_DER_X(int b0_, int n) {
  b0=b0_;
  N = new int[n+1];
  id = new POT_DER_IDENTIFIER*[n+1];
  for(int i=0;i<=n;i++) {
    N[i] = 0;
    id[i] = NULL;
  }
  return;
}


/*POT_DER_X Destructor: Deallocates memory*/
POT_DER_X::~POT_DER_X() {
  delete[] N;  N =NULL;
  delete[] id; id=NULL;
  return;
}


/*POT_DER_X Position: RETURNS THE POSITION TO INSERT NEW ENTRIES
returns alloc=true to indicate allocation needed*/
int POT_DER_X::Position(int N_atom, int *B, int **L, bool &alloc) {
  //CHECK FOR EMPTY LIST
  if(id[N_atom]==NULL) {
    alloc = true;
    N[N_atom] = 1;
    id[N_atom] = new POT_DER_IDENTIFIER(N_atom, B, L);
    return(0);
  }
  //FIND POSITION IN LIST
  int pos;
  pos = id[N_atom]->Position(N_atom-1, B, L, alloc);
  if(alloc) {                 //Memory has been or will be allocated
    N[N_atom] += 1;           //Increment number of entries
    if(pos==0) {              //Allocate memory at the begining of the list
      POT_DER_IDENTIFIER *ptr = id[N_atom];
      id[N_atom] = new POT_DER_IDENTIFIER(N_atom, B, L);
      id[N_atom]->next = ptr;
    }
  }
  return(pos);
}


/*POT_DER_X Allocate: ALLOCATES MEMORY FOR LINKED LIST IF NEEDED AND RETURNS POSITION*/
POT_DER_LINK_LIST *POT_DER_X::Allocate(bool alloc, int pos, int O_der, POT_DER_LINK_LIST **pdll) {
  POT_DER_LINK_LIST *next = (*pdll);
  if(pos==0) {                //HANDLE pos==0 SEPARATELY
    if(alloc) {
      (*pdll) = new POT_DER_LINK_LIST(O_der);
      pdll[0]->next = next;
      next = (*pdll);
    }
  }
  else {                      //NORMAL CASE
    POT_DER_LINK_LIST *prev = (*pdll);
    for(int i=1;i<pos;i++) {prev = prev->next;}
    next = prev->next;
    if(alloc) {
      POT_DER_LINK_LIST *ptr = new POT_DER_LINK_LIST(O_der);
      prev->next = ptr;
      ptr->next = next;
      next = ptr;
    }
  }

  return(next);
}


/*POT_DER_X Remove: Removes entry in list leaving others unchanged*/
void POT_DER_X::Remove(POT_DER_LINK_LIST *pdll) {
  cout <<"Enter Remove"<<endl;
/*  if(pdll==NULL) {return;}
  POT_DER_LINK_LIST *prev, *next;
  prev = pdll->prev;
  next = pdll->next;
  if(prev!=NULL) {prev->next = next;}
  if(next!=NULL) {next->prev = prev;}
  pdll->next = NULL;
  delete pdll;
  pdll=next;*/
  return;
}

/********************** Class PD_ENERGY_FORCE ***********************/
/*PD_ENERGY_FORCE Constructor*/
PD_ENERGY_FORCE::PD_ENERGY_FORCE(int b0_) : POT_DER_X(b0_, 1) {
  F0 = NULL;
  F1_0 = NULL;
  return;
}


/*PD_ENERGY_FORCE Destructor*/
PD_ENERGY_FORCE::~PD_ENERGY_FORCE() {
  delete F0;   F0=NULL;
  delete F1_0; F1_0=NULL;
  return;
}


/*PD_ENERGY_FORCE Compute: Computes the potential energy derivatives*/
void PD_ENERGY_FORCE::Compute(int d0, int n, int *B, int **L, double **X, POTENTIAL *pot) {
  //cout <<"Enter Compute"<<endl;
  //DECLARE LOCAL VARIABLES
  bool alloc;                 //Allocates memory if true
  int pos;                    //Position of entry in linked list
  int *b_list=NULL;           //List of atom indicies
  int **l_list=NULL;          //List of unit cell indicies
  double **DX=new double*[1]; //List positions for taking the derivative
  POT_DER_LINK_LIST *ptr;     //Temporary pointer to linked list location


  //COMPUTE HARMONIC FORCE CONSTANTS
  //F0 (Energy), F1_0 (-Force)
  pos = Position(1, b_list, l_list, alloc);
  DX[0] = X[d0];
  ptr = Allocate(alloc, pos, 0, &F0);
  ptr -> Compute(0, X, DX, pot);
  ptr = Allocate(alloc, pos, 1, &F1_0);
  ptr -> Compute(1, X, DX, pot);


  //DEALLOCATE MEMORY AND RETURN
  delete[] b_list; b_list=NULL;
  delete[] l_list; l_list=NULL;
  delete[] DX;     DX    =NULL;
  return;
}


/*PD_ENERGY_FORCE Send: Divides by mass and sends to other PEs*/
void PD_ENERGY_FORCE::Send() {
  //DECLARE LOCAL VARIABLES
  int rPE, i;                 //Counters
  //Pointers to energy and forces
  POT_DER_LINK_LIST *P0   = F0;
  POT_DER_LINK_LIST *P1_0 = F1_0;


  //NEGATE F1_0 TO MAKE FORCES
  for(i=0;i<N[1];i++) {
    P1_0->Scale(1, -1.0);
    P0   = P0  ->next;
    P1_0 = P1_0->next;
  }


  //SEND DATA TO OTHER PEs
#if defined PARALLEL
  for(rPE=0;rPE<nPE;rPE++) {
    if(rPE==PE) {continue;}
    MPI_Ssend(N, 2, MPI_INT, rPE, 1, MPI_COMM_WORLD);
  }
  for(rPE=0;rPE<nPE;rPE++) {
    if(rPE==PE) {continue;}
    P0   = F0;
    P1_0 = F1_0;
    for(i=0;i<N[1];i++) {
      MPI_Ssend(P0  ->fc, 1, MPI_DOUBLE, rPE, 0, MPI_COMM_WORLD);
      MPI_Ssend(P1_0->fc, DIM, MPI_DOUBLE, rPE, 10, MPI_COMM_WORLD);
      P0   = P0  ->next;
      P1_0 = P1_0->next;
    }
  }
#endif

  return;
}


/*PD_ENERGY_FORCE Recv: Receives data from other PEs*/
void PD_ENERGY_FORCE::Recv(int sPE) {
#if defined PARALLEL
  //DECLARE LOCAL VARIABLES
  int i;                      //Counters
  MPI_Status status;
  //Pointers to ID, energy, and forces
  F0   = new POT_DER_LINK_LIST(0);
  F1_0 = new POT_DER_LINK_LIST(1);
  POT_DER_LINK_LIST *P0   = F0;
  POT_DER_LINK_LIST *P1_0 = F1_0;


  //RECEIVE DATA
  MPI_Recv(N, 2, MPI_INT, sPE, 1, MPI_COMM_WORLD, &status);
  for(i=0;i<N[1];i++) {
    MPI_Recv(P0  ->fc, 1, MPI_DOUBLE, sPE, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(P1_0->fc, DIM, MPI_DOUBLE, sPE, 10, MPI_COMM_WORLD, &status);
    if(i<N[1]-1) {
      P0   = new POT_DER_LINK_LIST(0);
      P1_0 = new POT_DER_LINK_LIST(1);
      P0   = P0  ->next;
      P1_0 = P1_0->next;
    }
  }
#endif
  return;
}


/*PD_ENERGY_FORCE Print: OUTPUTS DATA TO STDOUT*/
void PD_ENERGY_FORCE::Print(int flag) {
  Log  <<"Energy and forces on atom "<<b0<<" = "<<F0->fc[0]*Parameter->energy<<" J,";
  cout <<"Energy and forces on atom "<<b0<<" = "<<F0->fc[0]*Parameter->energy<<" J,";
  for(int i=0;i<DIM;i++) {
    Log  <<" "<<F1_0->fc[i];
    cout <<" "<<F1_0->fc[i];
  }
  Log  <<endl;
  cout <<endl;
  if(flag) {
    Log  <<"N1 = "<<N[1]<<endl;
    cout <<"N1 = "<<N[1]<<endl;
  }
  return;
}

/************************ Class PD_HARMONIC *************************/
/*PD_HARMONIC Constructor*/
PD_HARMONIC::PD_HARMONIC(int b0_) : POT_DER_X(b0_, 2) {
  F2_00 = NULL;
  F2_01 = NULL;
  return;
}


/*PD_HARMONIC Destructor*/
PD_HARMONIC::~PD_HARMONIC() {
  delete F2_00; F2_00=NULL;
  delete F2_01; F2_01=NULL;
  return;
}


/*PD_HARMONIC Compute: Computes the potential energy derivatives*/
void PD_HARMONIC::Compute(int d0, int n, int *B, int **L, double **X, POTENTIAL *pot) {
  //DECLARE LOCAL VARIABLES
  bool alloc;                 //Allocates memory if true
  int pos;                    //Position of entry in linked list
  int d1;                     //Derivative index 1, 2, and 3
  int *b_list=new int[1];     //List of atom indicies
  int **l_list=new int*[1];   //List of unit cell indicies
  double **DX=new double*[2]; //List positions for taking the derivative
  POT_DER_LINK_LIST *ptr;     //Temporary pointer to linked list location


  //COMPUTE HARMONIC FORCE CONSTANTS
  //F2_00 (Self terms)
  pos = Position(1, b_list, l_list, alloc);
  DX[0] = DX[1] = X[d0];
  ptr = Allocate(alloc, pos, 2, &F2_00);
  ptr -> Compute(2, X, DX, pot);
  for(d1=0;d1<n;d1++) {
    if(d1==d0) {continue;}
    b_list[0] = B[d1];
    l_list[0] = L[d1];
    pos = Position(2, b_list, l_list, alloc);
    //F2_01
    DX[1] = X[d1];
    ptr = Allocate(alloc, pos, 2, &F2_01);
    ptr -> Compute(2, X, DX, pot);
  }


  //DEALLOCATE MEMORY AND RETURN
  delete[] b_list; b_list=NULL;
  delete[] l_list; l_list=NULL;
  delete[] DX;     DX    =NULL;
  return;
}


/*PD_HARMONIC Send: Divides by mass and sends to other PEs*/
void PD_HARMONIC::Send() {
  //DECLARE LOCAL VARIABLES
  int rPE, i;                 //Counters
  bool z;
  int *mat = Unit_Cell->mat; //Pointer to material
  double *M = Material->mass;//Pointer to mass
  double M0, M1;
  M0 = 1.0/sqrt(M[mat[b0]]);
  //Pointers to harmonic force constants
  POT_DER_IDENTIFIER *id2 = id[2];
  POT_DER_LINK_LIST *P2_00 = F2_00;
  POT_DER_LINK_LIST *P2_01 = F2_01;


  //DIVIDE BY MASSES
  //Lists containing 1 atom
  for(i=0;i<N[1];i++) {
    P2_00->Scale(2, M0*M0);
    P2_00 = P2_00->next;
  }
  //Lists containing 2 atoms
  for(i=0;i<N[2];i++) {
    z  = true;
    M1 = 1.0/sqrt(M[mat[id2->b[0]]]);
    if(!P2_01->Scale(2, M0*M1)) {z = false;}
/*    if(z) {
      N[2] -= 1;
      i  -= 1;
      Remove(P2_01);
    }
    else {*/
      id2 = id2->next;
      P2_01 = P2_01->next;
//    }
  }


  //SEND DATA
#if defined PARALLEL
  for(rPE=0;rPE<nPE;rPE++) {
    if(rPE==PE) {continue;}
    MPI_Ssend(N, 3, MPI_INT, rPE, 1, MPI_COMM_WORLD);
  }
  for(rPE=0;rPE<nPE;rPE++) {
    //Send lists containing 1 atom
    if(rPE==PE) {continue;}
    id2 = id[2];
    P2_00 = F2_00;
    P2_01 = F2_01;
    for(i=0;i<N[1];i++) {
      MPI_Ssend(P2_00->fc, DIM*DIM, MPI_DOUBLE, rPE, 200, MPI_COMM_WORLD);
      P2_00 = P2_00->next;
    }
    //Send lists containing 2 atoms
    for(i=0;i<N[2];i++) {
      MPI_Ssend(id2->b, 1, MPI_INT, rPE, 2, MPI_COMM_WORLD);
      MPI_Ssend(id2->l[0], DIM, MPI_INT, rPE, 3, MPI_COMM_WORLD);
      id2 = id2->next;
      MPI_Ssend(P2_01->fc, DIM*DIM, MPI_DOUBLE, rPE, 201, MPI_COMM_WORLD);
      P2_01 = P2_01->next;
    }
  }
#endif
  return;
}


/*PD_HARMONIC Recv: Receives data from other PEs*/
void PD_HARMONIC::Recv(int sPE) {
#if defined PARALLEL
  //RECEIVE NUMBER OF ENTRIES IN LISTS
  int i;                      //Counters
  MPI_Status status;
  MPI_Recv(N, 3, MPI_INT, sPE, 1, MPI_COMM_WORLD, &status);
  //RECEIVE LISTS CONTAINING 1 ATOM
  F2_00 = new POT_DER_LINK_LIST(2);
  POT_DER_LINK_LIST *P2_00 = F2_00;
  for(i=0;i<N[1];i++) {
    MPI_Recv(P2_00->fc, DIM*DIM, MPI_DOUBLE, sPE, 200, MPI_COMM_WORLD, &status);
    if(i<N[1]-1) {
      P2_00->next = new POT_DER_LINK_LIST(2);
      P2_00 = P2_00->next;
    }
  }
  //RECEIVE LISTS CONTAINING 2 ATOMS
  id[2] = new POT_DER_IDENTIFIER(2);
  POT_DER_IDENTIFIER *id2 = id[2];
  F2_01 = new POT_DER_LINK_LIST(2);
  POT_DER_LINK_LIST *P2_01 = F2_01;
  for(i=0;i<N[2];i++) {
    MPI_Recv(id2->b, 1, MPI_INT, sPE, 2, MPI_COMM_WORLD, &status);
    MPI_Recv(id2->l[0], DIM, MPI_INT, sPE, 3, MPI_COMM_WORLD, &status);
    MPI_Recv(P2_01->fc, DIM*DIM, MPI_DOUBLE, sPE, 201, MPI_COMM_WORLD, &status);
    if(i<N[2]-1) {
      id2->next = new POT_DER_IDENTIFIER(2);
      P2_01->next = new POT_DER_LINK_LIST(2);
      id2 = id2->next;
      P2_01 = P2_01->next;
    }
  }
#endif
  return;
}


/*PD_HARMONIC Print: OUTPUTS DATA TO STDOUT*/
void PD_HARMONIC::Print(int flag) {
  Log  <<"Harmonic force constants for atom "<<b0<<": N1, N2 = "<<N[1]<<", "<<N[2]<<endl;
  cout <<"Harmonic force constants for atom "<<b0<<": N1, N2 = "<<N[1]<<", "<<N[2]<<endl;
  if(flag) {
    //Print actual force constants
    ;
  }
  return;
}

/*********************** Class PD_ANHARMONIC ************************/
/*PD_ANHARMONIC Constructor*/
PD_ANHARMONIC::PD_ANHARMONIC(int b0_) : POT_DER_X(b0_, 4) {
  F3_000  = NULL;
  F4_0000 = NULL;
  F3_001  = NULL;
  F3_011  = NULL;
  F4_0001 = NULL;
  F4_0011 = NULL;
  F4_0111 = NULL;
  F3_012  = NULL;
  F4_0012 = NULL;
  F4_0112 = NULL;
  F4_0122 = NULL;
  F4_0123 = NULL;
  return;
}


/*PD_ANHARMONIC Destructor*/
PD_ANHARMONIC::~PD_ANHARMONIC() {
  delete F3_000;  F3_000 =NULL;
  delete F4_0000; F4_0000=NULL;
  delete F3_001;  F3_001 =NULL;
  delete F3_011;  F3_011 =NULL;
  delete F4_0001; F4_0001=NULL;
  delete F4_0011; F4_0011=NULL;
  delete F4_0111; F4_0111=NULL;
  delete F3_012;  F3_012 =NULL;
  delete F4_0012; F4_0012=NULL;
  delete F4_0112; F4_0112=NULL;
  delete F4_0122; F4_0122=NULL;
  delete F4_0123; F4_0123=NULL;
  return;
}


/*PD_ANHARMONIC Compute: Computes the potential energy derivatives*/
void PD_ANHARMONIC::Compute(int d0, int n, int *B, int **L, double **X, POTENTIAL *pot) {
  //DECLARE LOCAL VARIABLES
  bool alloc;                 //Allocates memory if true
  int pos;                    //Position of entry in linked list
  int d1, d2, d3;             //Derivative index 1, 2, and 3
  int *b_list=new int[3];     //List of atom indicies
  int **l_list=new int*[3];   //List of unit cell indicies
  double **DX=new double*[4]; //List positions for taking the derivative
  POT_DER_LINK_LIST *ptr;     //Temporary pointer to linked list location


  //COMPUTE THIRD AND FOURTH ORDER FORCE CONSTANTS
  //F3_000, F4_0000 (Self terms)
  pos = Position(1, b_list, l_list, alloc);
  DX[0] = DX[1] = DX[2] = DX[3] = X[d0];
  ptr = Allocate(alloc, pos, 3, &F3_000);
  ptr -> Compute(3, X, DX, pot);
  ptr = Allocate(alloc, pos, 4, &F4_0000);
#if defined FREQ_SHIFT
  ptr -> Compute(4, X, DX, pot);
#endif
  //All other 3rd and 4th order force constants
  for(d1=0;d1<n;d1++) {
    if(d1==d0) {continue;}
    b_list[0] = B[d1];
    l_list[0] = L[d1];
    pos = Position(2, b_list, l_list, alloc);
    //F4_0001 (F4_0010, F4_0100)
    DX[1] = DX[2] = X[d0];
    DX[3] = X[d1];
    ptr = Allocate(alloc, pos, 4, &F4_0001);
#if defined FREQ_SHIFT
    ptr -> Compute(4, X, DX, pot);
#endif
    //F3_001, F4_0011 (F3_010, F4_0101, F4_0110)
    DX[2] = X[d1];
    ptr = Allocate(alloc, pos, 3, &F3_001);
    ptr -> Compute(3, X, DX, pot);
    ptr = Allocate(alloc, pos, 4, &F4_0011);
#if defined FREQ_SHIFT
    ptr -> Compute(4, X, DX, pot);
#endif
    //F3_011, F4_0111
    DX[1] = X[d1];
    ptr = Allocate(alloc, pos, 3, &F3_011);
    ptr -> Compute(3, X, DX, pot);
    ptr = Allocate(alloc, pos, 4, &F4_0111);
#if defined FREQ_SHIFT
    ptr -> Compute(4, X, DX, pot);
#endif
    for(d2=d1+1;d2<n;d2++) {
      if(d2==d0) {continue;}
      b_list[1] = B[d2];
      l_list[1] = L[d2];
      pos = Position(3, b_list, l_list, alloc);
      //F4_0012 (F4_0102, F4_0120, F4_0021, F4_0201, F4_0210)
      DX[1] = X[d0];
      DX[2] = X[d1];
      DX[3] = X[d2];
      ptr = Allocate(alloc, pos, 4, &F4_0012);
#if defined FREQ_SHIFT
      ptr -> Compute(4, X, DX, pot);
#endif
      //F4_0112 (F4_0121, F4_0211)
      DX[1] = X[d1];
      ptr = Allocate(alloc, pos, 4, &F4_0112);
#if defined FREQ_SHIFT
      ptr -> Compute(4, X, DX, pot);
#endif
      //F3_012, F4_0122 (F3_021, F4_0212, F4_0221)
      DX[2] = X[d2];
      ptr = Allocate(alloc, pos, 3, &F3_012);
      ptr -> Compute(3, X, DX, pot);
      ptr = Allocate(alloc, pos, 4, &F4_0122);
#if defined FREQ_SHIFT
      ptr -> Compute(4, X, DX, pot);
#endif
      for(d3=d2+1;d3<n;d3++) {
        if(d3==d0) {continue;}
        b_list[2] = B[d3];
        l_list[2] = L[d3];
        pos = Position(4, b_list, l_list, alloc);
        //F4_0123 (F4_0231, F4_0312, F4_0321, F4_0213, F4_0132)
        DX[3] = X[d3];
        ptr = Allocate(alloc, pos, 4, &F4_0123);
#if defined FREQ_SHIFT
        ptr -> Compute(4, X, DX, pot);
#endif
      }
    }
  }


  //DEALLOCATE MEMORY AND RETURN
  delete[] b_list; b_list=NULL;
  delete[] l_list; l_list=NULL;
  delete[] DX;     DX    =NULL;
  return;
}


/*PD_ANHARMONIC Send: Divides by mass and sends to other PEs*/
void PD_ANHARMONIC::Send() {
  //DECLARE LOCAL VARIABLES
  int rPE, i;                 //Counters
  bool z;                     //True if all force constant entries are zero
  int *mat = Unit_Cell->mat;  //Pointer to material
  double *M = Material->mass; //Pointer to mass
  double M0, M1, M2, M3;
  M0 = 1.0/sqrt(M[mat[b0]]);
  //Pointers to anharmonic force constants
  POT_DER_LINK_LIST *P3_000  = F3_000;
  POT_DER_LINK_LIST *P4_0000 = F4_0000;
  POT_DER_IDENTIFIER *id2 = id[2];
  POT_DER_LINK_LIST *P3_001  = F3_001;
  POT_DER_LINK_LIST *P3_011  = F3_011;
  POT_DER_LINK_LIST *P4_0001 = F4_0001;
  POT_DER_LINK_LIST *P4_0011 = F4_0011;
  POT_DER_LINK_LIST *P4_0111 = F4_0111;
  POT_DER_IDENTIFIER *id3 = id[3];
  POT_DER_LINK_LIST *P3_012  = F3_012;
  POT_DER_LINK_LIST *P4_0012 = F4_0012;
  POT_DER_LINK_LIST *P4_0112 = F4_0112;
  POT_DER_LINK_LIST *P4_0122 = F4_0122;
  POT_DER_IDENTIFIER *id4 = id[4];
  POT_DER_LINK_LIST *P4_0123 = F4_0123;


  //DIVIDE BY MASSES AND ZEROS SMALL VALUES
  //Lists containing 1 atom
  for(i=0;i<N[1];i++) {
    P3_000 ->Scale(3, M0*M0*M0);
    P4_0000->Scale(4, M0*M0*M0*M0);
    P3_000  = P3_000 ->next;
    P4_0000 = P4_0000->next;
  }
  //Lists containing 2 atoms
  for(i=0;i<N[2];i++) {
    z = true;
    M1 = 1.0/sqrt(M[mat[id2->b[0]]]);
    if(!P3_001 ->Scale(3, M0*M0*M1)) {z = false;}
    if(!P3_011 ->Scale(3, M0*M1*M1)) {z = false;}
    if(!P4_0001->Scale(4, M0*M0*M0*M1)) {z = false;}
    if(!P4_0011->Scale(4, M0*M0*M1*M1)) {z = false;}
    if(!P4_0111->Scale(4, M0*M1*M1*M1)) {z = false;}
/*      if(z) {
      N[2] -= 1;
      i  -= 1;
      Remove(P3_001);
      Remove(P3_011);
      Remove(P4_0001);
      Remove(P4_0011);
      Remove(P4_0111);
    }
    else {*/
      id2 = id2->next;
      P3_001  = P3_001 ->next;
      P3_011  = P3_011 ->next;
      P4_0001 = P4_0001->next;
      P4_0011 = P4_0011->next;
      P4_0111 = P4_0111->next;
//      }
  }
  //Lists containing 3 atoms
  for(i=0;i<N[3];i++) {
    z  = true;
    M1 = 1.0/sqrt(M[mat[id3->b[0]]]);
    M2 = 1.0/sqrt(M[mat[id3->b[1]]]);
    if(!P3_012 ->Scale(3, M0*M1*M2))    {z = false;}
    if(!P4_0012->Scale(4, M0*M0*M1*M2)) {z = false;}
    if(!P4_0112->Scale(4, M0*M1*M1*M2)) {z = false;}
    if(!P4_0122->Scale(4, M0*M1*M2*M2)) {z = false;}
/*      if(z) {
      N[3] -= 1;
      i  -= 1;
      Remove(P3_012);
      Remove(P4_0012);
      Remove(P4_0112);
      Remove(P4_0122);
    }
    else {*/
      id3 = id3->next;
      P3_012  = P3_012 ->next;
      P4_0012 = P4_0012->next;
      P4_0112 = P4_0112->next;
      P4_0122 = P4_0122->next;
//      }
  }
  //Lists containing 4 atoms
  for(i=0;i<N[4];i++) {
    z = true;
    M1 = 1.0/sqrt(M[mat[id4->b[0]]]);
    M2 = 1.0/sqrt(M[mat[id4->b[1]]]);
    M3 = 1.0/sqrt(M[mat[id4->b[2]]]);
    if(!P4_0123->Scale(4, M0*M1*M2*M3)) {z = false;}
/*      if(z) {
      N[4] -= 1;
      i  -= 1;
      Remove(P4_0123);
    }
    else {*/
      id4 = id4->next;
      P4_0123 = P4_0123->next;
//      }
  }


#if defined PARALLEL
  //SEND DATA
  for(rPE=0;rPE<nPE;rPE++) {
    if(rPE==PE) {continue;}
    MPI_Ssend(N, 5, MPI_INT, rPE, 1, MPI_COMM_WORLD);
  }
  for(rPE=0;rPE<nPE;rPE++) {
    if(rPE==PE) {continue;}
    id2 = id[2];
    id3 = id[3];
    id4 = id[4];
    P3_000  = F3_000;
    P4_0000 = F4_0000;
    P3_001  = F3_001;
    P3_011  = F3_011;
    P4_0001 = F4_0001;
    P4_0011 = F4_0011;
    P4_0111 = F4_0111;
    P3_012  = F3_012;
    P4_0012 = F4_0012;
    P4_0112 = F4_0112;
    P4_0122 = F4_0122;
    P4_0123 = F4_0123;
    //Send lists containing 1 atom
    for(i=0;i<N[1];i++) {
      MPI_Ssend(P3_000 ->fc, DIM*DIM*DIM, MPI_DOUBLE, rPE, 3000, MPI_COMM_WORLD);
      MPI_Ssend(P4_0000->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, rPE, 40000, MPI_COMM_WORLD);
      P3_000  = P3_000 ->next;
      P4_0000 = P4_0000->next;
    }
    //Send lists containing 2 atoms
    for(i=0;i<N[2];i++) {
      MPI_Ssend(id2->b, 1, MPI_INT, rPE, 2, MPI_COMM_WORLD);
      MPI_Ssend(id2->l[0], DIM, MPI_INT, rPE, 3, MPI_COMM_WORLD);
      id2 = id2->next;
      MPI_Ssend(P3_001 ->fc, DIM*DIM*DIM, MPI_DOUBLE, rPE, 3001, MPI_COMM_WORLD);
      MPI_Ssend(P3_011 ->fc, DIM*DIM*DIM, MPI_DOUBLE, rPE, 3011, MPI_COMM_WORLD);
      MPI_Ssend(P4_0001->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, rPE, 40001, MPI_COMM_WORLD);
      MPI_Ssend(P4_0011->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, rPE, 40011, MPI_COMM_WORLD);
      MPI_Ssend(P4_0111->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, rPE, 40111, MPI_COMM_WORLD);
      P3_001  = P3_001 ->next;
      P3_011  = P3_011 ->next;
      P4_0001 = P4_0001->next;
      P4_0011 = P4_0011->next;
      P4_0111 = P4_0111->next;
    }
    //Send lists containing 3 atoms
    for(i=0;i<N[3];i++) {
      MPI_Ssend(id3->b, 2, MPI_INT, rPE, 2, MPI_COMM_WORLD);
      MPI_Ssend(id3->l[0], 2*DIM, MPI_INT, rPE, 3, MPI_COMM_WORLD);
      id3 = id3->next;
      MPI_Ssend(P3_012 ->fc, DIM*DIM*DIM, MPI_DOUBLE, rPE, 3012, MPI_COMM_WORLD);
      MPI_Ssend(P4_0012->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, rPE, 40012, MPI_COMM_WORLD);
      MPI_Ssend(P4_0112->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, rPE, 40112, MPI_COMM_WORLD);
      MPI_Ssend(P4_0122->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, rPE, 40122, MPI_COMM_WORLD);
      P3_012  = P3_012 ->next;
      P4_0012 = P4_0012->next;
      P4_0112 = P4_0112->next;
      P4_0122 = P4_0122->next;
    }
    //Send lists containing 4 atoms
    for(i=0;i<N[4];i++) {
      MPI_Ssend(id4->b, 3, MPI_INT, rPE, 2, MPI_COMM_WORLD);
      MPI_Ssend(id4->l[0], 3*DIM, MPI_INT, rPE, 3, MPI_COMM_WORLD);
      id4 = id4->next;
      MPI_Ssend(P4_0123->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, rPE, 40123, MPI_COMM_WORLD);
      P4_0123 = P4_0123->next;
    }
  }
#endif

  return;
}


/*PD_ANHARMONIC Recv: Receives data from other PEs*/
void PD_ANHARMONIC::Recv(int sPE) {
#if defined PARALLEL
  //RECEIVE NUMBER OF ENTRIES IN LISTS
  int i;                      //Counters
  MPI_Status status;
  MPI_Recv(N, 5, MPI_INT, sPE, 1, MPI_COMM_WORLD, &status);
  //RECEIVE LISTS CONTAINING 1 ATOM
  F3_000  = new POT_DER_LINK_LIST(3);
  F4_0000 = new POT_DER_LINK_LIST(4);
  POT_DER_LINK_LIST *P3_000  = F3_000;
  POT_DER_LINK_LIST *P4_0000 = F4_0000;
  for(i=0;i<N[1];i++) {
    MPI_Recv(P3_000 ->fc, DIM*DIM*DIM, MPI_DOUBLE, sPE, 3000, MPI_COMM_WORLD, &status);
    MPI_Recv(P4_0000->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, sPE, 40000, MPI_COMM_WORLD, &status);
    if(i<N[1]-1) {
      P3_000 ->next = new POT_DER_LINK_LIST(3);
      P4_0000->next = new POT_DER_LINK_LIST(4);
      P3_000  = P3_000 ->next;
      P4_0000 = P4_0000->next;
    }
  }
  //RECEIVE LISTS CONTAINING 2 ATOMS
  id[2] = new POT_DER_IDENTIFIER(2);
  POT_DER_IDENTIFIER *id2 = id[2];
  F3_001  = new POT_DER_LINK_LIST(3);
  F3_011  = new POT_DER_LINK_LIST(3);
  F4_0001 = new POT_DER_LINK_LIST(4);
  F4_0011 = new POT_DER_LINK_LIST(4);
  F4_0111 = new POT_DER_LINK_LIST(4);
  POT_DER_LINK_LIST *P3_001  = F3_001;
  POT_DER_LINK_LIST *P3_011  = F3_011;
  POT_DER_LINK_LIST *P4_0001 = F4_0001;
  POT_DER_LINK_LIST *P4_0011 = F4_0011;
  POT_DER_LINK_LIST *P4_0111 = F4_0111;
  for(i=0;i<N[2];i++) {
    MPI_Recv(id2->b, 1, MPI_INT, sPE, 2, MPI_COMM_WORLD, &status);
    MPI_Recv(id2->l[0], DIM, MPI_INT, sPE, 3, MPI_COMM_WORLD, &status);
    MPI_Recv(P3_001 ->fc, DIM*DIM*DIM, MPI_DOUBLE, sPE, 3001, MPI_COMM_WORLD, &status);
    MPI_Recv(P3_011 ->fc, DIM*DIM*DIM, MPI_DOUBLE, sPE, 3011, MPI_COMM_WORLD, &status);
    MPI_Recv(P4_0001->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, sPE, 40001, MPI_COMM_WORLD, &status);
    MPI_Recv(P4_0011->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, sPE, 40011, MPI_COMM_WORLD, &status);
    MPI_Recv(P4_0111->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, sPE, 40111, MPI_COMM_WORLD, &status);
    if(i<N[2]-1) {
      id2->next = new POT_DER_IDENTIFIER(2);
      P3_001 ->next = new POT_DER_LINK_LIST(3);
      P3_011 ->next = new POT_DER_LINK_LIST(3);
      P4_0001->next = new POT_DER_LINK_LIST(4);
      P4_0011->next = new POT_DER_LINK_LIST(4);
      P4_0111->next = new POT_DER_LINK_LIST(4);
      id2 = id2->next;
      P3_001  = P3_001 ->next;
      P3_011  = P3_011 ->next;
      P4_0001 = P4_0001->next;
      P4_0011 = P4_0011->next;
      P4_0111 = P4_0111->next;
    }
  }
  //RECEIVE LISTS CONTAINING 3 ATOMS
  id[3] = new POT_DER_IDENTIFIER(3);
  POT_DER_IDENTIFIER *id3 = id[3];
  F3_012  = new POT_DER_LINK_LIST(3);
  F4_0012 = new POT_DER_LINK_LIST(4);
  F4_0112 = new POT_DER_LINK_LIST(4);
  F4_0122 = new POT_DER_LINK_LIST(4);
  POT_DER_LINK_LIST *P3_012  = F3_012;
  POT_DER_LINK_LIST *P4_0012 = F4_0012;
  POT_DER_LINK_LIST *P4_0112 = F4_0112;
  POT_DER_LINK_LIST *P4_0122 = F4_0122;
  for(i=0;i<N[3];i++) {
    MPI_Recv(id3->b, 2, MPI_INT, sPE, 2, MPI_COMM_WORLD, &status);
    MPI_Recv(id3->l[0], 2*DIM, MPI_INT, sPE, 3, MPI_COMM_WORLD, &status);
    MPI_Recv(P3_012 ->fc, DIM*DIM*DIM, MPI_DOUBLE, sPE, 3012, MPI_COMM_WORLD, &status);
    MPI_Recv(P4_0012->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, sPE, 40012, MPI_COMM_WORLD, &status);
    MPI_Recv(P4_0112->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, sPE, 40112, MPI_COMM_WORLD, &status);
    MPI_Recv(P4_0122->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, sPE, 40122, MPI_COMM_WORLD, &status);
    if(i<N[3]-1) {
      id3->next = new POT_DER_IDENTIFIER(3);
      P3_012 ->next = new POT_DER_LINK_LIST(3);
      P4_0012->next = new POT_DER_LINK_LIST(4);
      P4_0112->next = new POT_DER_LINK_LIST(4);
      P4_0122->next = new POT_DER_LINK_LIST(4);
      id3 = id3->next;
      P3_012  = P3_012 ->next;
      P4_0012 = P4_0012->next;
      P4_0112 = P4_0112->next;
      P4_0122 = P4_0122->next;
    }
  }
  //RECEIVE LISTS CONTAINING 4 ATOMS
  id[4] = new POT_DER_IDENTIFIER(4);
  POT_DER_IDENTIFIER *id4 = id[4];
  F4_0123 = new POT_DER_LINK_LIST(4);
  POT_DER_LINK_LIST *P4_0123 = F4_0123;
  for(i=0;i<N[4];i++) {
    MPI_Recv(id4->b, 3, MPI_INT, sPE, 2, MPI_COMM_WORLD, &status);
    MPI_Recv(id4->l[0], 3*DIM, MPI_INT, sPE, 3, MPI_COMM_WORLD, &status);
    MPI_Recv(P4_0123->fc, DIM*DIM*DIM*DIM, MPI_DOUBLE, sPE, 40123, MPI_COMM_WORLD, &status);
    if(i<N[4]-1) {
      id4->next = new POT_DER_IDENTIFIER(4);
      P4_0123->next = new POT_DER_LINK_LIST(4);
      id4 = id4->next;
      P4_0123 = P4_0123->next;
    }
  }
#endif
  return;
}


/*PD_ANHARMONIC Print: OUTPUTS DATA TO STDOUT*/
void PD_ANHARMONIC::Print(int flag) {
  Log  <<"Anharmonic force constants for atom "<<b0<<": N1, N2, N3, N4 = "<<N[1]<<", "<<N[2]<<", "<<N[3]<<", "<<N[4]<<endl;
  cout <<"Anharmonic force constants for atom "<<b0<<": N1, N2, N3, N4 = "<<N[1]<<", "<<N[2]<<", "<<N[3]<<", "<<N[4]<<endl;
  if(flag) {
    //Print actual force constants
    int i, d;
    POT_DER_LINK_LIST *P3_000  = F3_000;
    POT_DER_LINK_LIST *P4_0000 = F4_0000;
    POT_DER_IDENTIFIER *id2 = id[2];
    POT_DER_LINK_LIST *P3_001  = F3_001;
    POT_DER_LINK_LIST *P3_011  = F3_011;
    POT_DER_LINK_LIST *P4_0001 = F4_0001;
    POT_DER_LINK_LIST *P4_0011 = F4_0011;
    POT_DER_LINK_LIST *P4_0111 = F4_0111;
    POT_DER_IDENTIFIER *id3 = id[3];
    POT_DER_LINK_LIST *P3_012  = F3_012;
    POT_DER_LINK_LIST *P4_0012 = F4_0012;
    POT_DER_LINK_LIST *P4_0112 = F4_0112;
    POT_DER_LINK_LIST *P4_0122 = F4_0122;
    POT_DER_IDENTIFIER *id4 = id[4];
    POT_DER_LINK_LIST *P4_0123 = F4_0123;
    for(i=0;i<N[1];i++) {
//      Log <<id[1]->b[0];
//      for(d=0;d<DIM;d++) {Log <<" "<<id[1]->l[0][d];}
//      Log <<endl;
      Log <<"P4_0000 = ";
      P4_0000->Print(4);
      P4_0000 = P4_0000->next;
    }
    for(i=0;i<N[2];i++) {
      Log <<id2->b[0];
      for(d=0;d<DIM;d++) {Log <<" "<<id2->l[0][d];}
      Log <<endl;
      Log <<"P4_0001 = ";
      P4_0001->Print(4);
      Log <<"P4_0011 = ";
      P4_0011->Print(4);
      Log <<"P4_0111 = ";
      P4_0111->Print(4);
      id2 = id2->next;
      P4_0001 = P4_0001->next;
      P4_0011 = P4_0011->next;
      P4_0111 = P4_0111->next;
    }
  }
  return;
}
