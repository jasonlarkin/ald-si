/*								           Symmetry.cpp		               					*/
/*                            12/06/2008                            */
/*********************************************************************
*    Function definitions for the SYMMETRY class.                    *
*********************************************************************/

/*DEFINE HEADERS*/
#include "LDCode.h"
#include <cmath>


/*SYMMETRY CLASS CONSTRUCTOR: USES LATTICE AND UNIT CELL TO DETERMINE CRYSTAL SYMMETRY*/
SYMMETRY::SYMMETRY() {
  k_test = new int[3];
  dir = new int[3];
  for(int i=0;i<3;i++) {dir[i] = -1;}
  Sym_ptr = NULL;
  strcpy(Symb, "1");
  wv_out = false;
  return;
}


/*DESTRUCTOR FOR THE SYMMETRY CLASS*/
SYMMETRY::~SYMMETRY() {
  delete[] k_test;      k_test     =NULL;
  delete[] dir;         dir        =NULL;
  delete[] unique[0];   unique[0]  =NULL;
  delete[] unique;      unique     =NULL;
  delete[] F_map[0][0]; F_map[0][0]=NULL;
  delete[] F_map[0];    F_map[0]   =NULL;
  delete[] F_map;       F_map      =NULL;
  delete[] E_map[0][0]; E_map[0][0]=NULL;
  delete[] E_map[0];    E_map[0]   =NULL;
  delete[] E_map;       E_map      =NULL;
  return;
}


/*FUNCTION Define: READS INPUT DATA*/
string SYMMETRY::ReadInput(ifstream &Input) {
  int i;
  string str;
  string Read_Next(ifstream &Input);
  strcpy(Symb, Read_Next(Input).c_str());
  while(!Input.eof()) {
    //identify keyword in HARMONIC category
    str = Read_Next(Input);
    i = 0;
    while(str[i]) {str[i]=tolower(str[i]);i++;}
    if(str.compare(0, 3, "wave_vector", 3)==0) {
      str = Read_Next(Input);
      if     (str.compare(0, 1, "true",  1)==0) {wv_out = true;}
      else if(str.compare(0, 1, "false", 1)==0) {wv_out = false;}
      else {Log <<"'wave_vector' accepts a flag of true or false"<<endl;exit(0);}
    }
    else{break;}
  }
  return(str);
}


/*SYMMETRY FUNCTION Output: WRITES DATA TO LOG FILE*/
void SYMMETRY::Output() {
  Log <<"\nSYMMETRY = "<<Symb<<endl;
  Log <<"  wave_vector = "<<(wv_out ? "true" : "false")<<endl;
  return;
}


/*SYMMETRY FUNCTION Initialize: INITIALIZES VARIALES*/
void SYMMETRY::Initialize() {
  //DECLARE LOCAL VARIABLES
  int i, j;
  int p, m, swap;
  char c_swap;
  int *UC  = Lattice->N;
  int *low = Lattice->lower;
  int *up  = Lattice->upper;


  //INTERPERT SYMMETRY
  sym = atoi(Symb);
  if(sym<1) {strcpy(Symb, "1");  sym = 1;}
  p = 1;
  for(i=0;i<DIM;i++) {p *= 10;}
  for(i=0;i<DIM;i++) {
    p = p/10;
    m = sym%p;
    dir[i] = (sym-m)/p - 1;
    if(dir[i]>=DIM) {
      Log <<"Symmetry direction greater than dimension."<<endl;exit(0);
    }
    sym = m;
    //Order numbers
    for(j=0;j<i;j++) {
      if(dir[i]<dir[j]) {
        swap = dir[i];
        dir[i] = dir[j];
        dir[j] = swap;
        if(swap!=-1) {
          c_swap = Symb[i];
          Symb[i] = Symb[j];
          Symb[j] = c_swap;
        }
      }
    }
  }
  sym = DIM;
  m = 0;
  for(i=0;i<DIM;i++) {m -= dir[i];}
  if(m==DIM) {sym=0;for(i=0;i<DIM;i++) {dir[i]=i;}}//All entries -1
  else {
    while(dir[0]==-1) {
      for(i=0;i<DIM-1;i++) {dir[i] = dir[i+1];}
      dir[DIM-1] = -1;
      if( (sym-=1)<=0 ) {sym=0;break;}
    }
    m = sym;
    if(dir[0]==-1) {for(i=0;i<DIM;i++) {dir[m] = i; m++;}}
    else {for(i=0;i<dir[0];i++) {dir[m] = i; m++;}}
    for(i=0;i<sym-1;i++) {
      if(m>=DIM) {break;}
      if(dir[i+1]==-1) {for(j=1;j<DIM-dir[i];j++) {dir[m] = dir[i] + j; m++;}}
      else {for(j=1;j<dir[i+1]-dir[i];j++) {dir[m] = dir[i] + j; m++;}}
    }
    for(i=m;i<DIM;i++) {dir[i] = i;}
  }


  //VERIFY THAT SYMMETRY IS POSSIBLE
  double *a_mag = new double[DIM];
  for(i=0;i<DIM;i++) {
    a_mag[i] = 0.0;
    for(j=0;j<DIM;j++) {a_mag[i] += Lattice->a[i][j]*Lattice->a[i][j];}
    a_mag[i] = sqrt(a_mag[i]);
  }
  if(sym==SYM3D) {
    if( (UC[0]!=UC[1])||(UC[1]!=UC[2])||(UC[0]!=UC[2]) ) {
      Log <<"Symmetry not possible (N)."<<endl;exit(0);
    }
    if( (fabs(a_mag[0]-a_mag[1])>1.0e-6)||(fabs(a_mag[1]-a_mag[2])>1.0e-6)
        ||(fabs(a_mag[0]-a_mag[2])>1.0e-6) ) {
      Log <<"Symmetry not possible (a)."<<endl;exit(0);
    }
  } else
  if(sym==SYM2D) {
    if(UC[dir[0]]!=UC[dir[1]]) {Log <<"Symmetry not possible (N)."<<endl;exit(0);}
    if(fabs(a_mag[dir[0]]-a_mag[dir[1]])>1.0e-6) {
      Log <<"Symmetry not possible (a)."<<endl;exit(0);
    }
  }
  delete[] a_mag;  a_mag = NULL;


  //INITIALIZE VARIABLES
	TotalDOF = UC_DOF;
	for(i=0;i<DIM;i++) {
		TotalDOF *= UC[i];
		k_test[i] = (UC[i]+1)/2;
	}
	Sym_E = TotalDOF;
  i = UC[1] - ((UC[1]+1)%2);
  j = UC[2] - ((UC[2]+1)%2);
  Sym_E += UC_DOF*(low[2]+low[1]*j+low[0]*i*j);
	if(sym==SYM3D) {
	  Sym_F = UC_DOF*( (up[0]+1)*(up[0]+2)*(up[0]+3) )/6;
	  Sym_ptr = &SYMMETRY::Sym_3D;
	} else
	if(sym==SYM2D) {
	  Sym_F = UC_DOF*(up[dir[2]]+1)*( (up[dir[0]]+1)*(up[dir[0]]+2) )/2;
	  Sym_ptr = &SYMMETRY::Sym_2D;
	}
	else {
	  Sym_F = Sym_E;
	  Sym_ptr = &SYMMETRY::Sym_None;
  }


  //CREATE LIST OF UNIQUE PHONONS AND MAPPING LIST
  unique = new int*[Sym_F];
  unique[0] = new int[DIM*Sym_F];
  F_map = new int** [UC[0]]; //Directs wavevector to proper frequency
	F_map[0] = new int* [UC[0]*UC[1]];
	E_map = new int** [UC[0]]; //Directs wavevector to proper eigenvector
	E_map[0] = new int* [UC[0]*UC[1]];
	int *wv = new int[DIM];
	int *n = new int[DIM];
	int inner_F = 0;
	int outer_E = 0;
	m = 0;
	for(wv[0]=up[0];wv[0]>=low[0];wv[0]--) {
	  n[0] = up[0] - wv[0];
	  if(wv[0]>low[0]) {
      F_map[n[0]+1] = F_map[n[0]] + UC[1];
      E_map[n[0]+1] = E_map[n[0]] + UC[1];
	  }
	  F_map[n[0]][0] = new int [UC[1]*UC[2]];
	  E_map[n[0]][0] = new int [UC[1]*UC[2]];
	  for(wv[1]=up[1];wv[1]>=low[1];wv[1]--) {
      n[1] = up[1] - wv[1];
      if(wv[1]>low[1]) {
        F_map[n[0]][n[1]+1] = F_map[n[0]][n[1]] + UC[2];
        E_map[n[0]][n[1]+1] = E_map[n[0]][n[1]] + UC[2];
      }
      for(wv[2]=up[2];wv[2]>=low[2];wv[2]--) {
        n[2] = up[2] - wv[2];
        //Find location of eigenvector in list
        if(Symmetry->Sym_Negative(wv)<0.0) {
          E_map[n[0]][n[1]][n[2]] = E_map[up[0]+wv[0]][up[1]+wv[1]][up[2]+wv[2]];
        }
        else {
          E_map[n[0]][n[1]][n[2]] = outer_E;
          outer_E += UC_DOF;
        }
        //Find location of frequency in list
        if(!Symmetry->Sym_Freq(wv, &i)) {
          if(m<Sym_F-1) {unique[m+1] = unique[m] + DIM;}
          for(j=0;j<DIM;j++) {unique[m][j] = wv[j];}
          m += 1;
          inner_F += UC_DOF;
        }
        F_map[n[0]][n[1]][n[2]] = i;
      }
	  }
	}
	if(Sym_F!=m*UC_DOF) {Log <<"SYMMETRY ERROR"<<endl;exit(0);}


	if(wv_out) {
	  double x;
	  ofstream output("Symmetry.txt");
	  output <<"List of "<<Sym_F/UC_DOF<<" unique wave vectors out of a total of "
           <<TotalDOF/UC_DOF<<" in the order output by the HARMONIC and "
           <<"ANHARMONIC calculations.  The first "<<DIM<<" values are "
           <<"the 'l's in the expression l*b/N, where b is the reciprocal "
           <<"lattice vector and N is the number unit cells.  The next "
           <<DIM<<" numbers are the Cartesian values of the recipocal "
           <<"wave vectors in 1/nm."<<endl;
    for(m=0;m<Sym_F/UC_DOF;m++) {
      for(j=0;j<DIM;j++) {output <<unique[m][j]<<'\t';}
      for(j=0;j<DIM;j++) {
        x = 0.0;
        for(i=0;i<DIM;i++) {x += Lattice->b[i][j]*unique[m][i];}
        output <<'\t'<<x*1.0e-9/Parameter->length;
      }
      output <<endl;
    }
    output.close();
	}

	return;
}


/*SYMMETRY FUNCTION Sym_Negative: TESTS IF NEGATIVE OF WAVE VECTOR IS
  CONSIDERED AND RETURNS THE FACTOR USED IN 4-PHONON INTERACTIONS*/
double SYMMETRY::Sym_Negative(int *k) {
  //COMPUTE FACTOR USED IN 4-PHONON INTERACTIONS
  //RETURNS -1.0 IF THE WAVE VECTOR IS CONSIDERED AT ANOTHER TIME
	if(k[0]==k_test[0]) {return(1.0);}
	if(k[1]==k_test[1]) {return(1.0);}
	if(k[2]==k_test[2]) {return(1.0);}
	if(k[0]<0) {return(-1.0);}  //Symmetry condition k=-k
	if(k[0]==0) {
		if(k[1]<0) {return(-1.0);}//Symmetry condition k=-k
		if(k[1]==0) {
		  if(k[2]<0) {return(-1.0);}//Symmetry condition k=-k
		  if(k[2]==0) {return(1.0);}
    }
	}
	return(2.0);
}


/*SYMMETRY FUNCTION Sym_3Phonon: COMPUTES FACTOR USED IN 3-PHONON
  INTERACTIONS OR RETURNS -1.0 IF CONSIDERED AT ANOTHER TIME*/
double SYMMETRY::Sym_3Phonon(int *k1, int *k2) {
  //COMPUTE FACTOR USED IN 3-PHONON INTERACTIONS
  if(k1[0]>k2[0]) {return(2.0);}
  else if(k1[0]==k2[0]) {
    if(k1[1]>k2[1]) {return(2.0);}
    else if(k1[1]==k2[1]) {
      if(k1[2]>k2[2]) {return(2.0);}
      else if(k1[2]==k2[2]) {return(1.0);}
    }
  }
  return(-1.0);
}


/*SYMMETRY FUNCTION Sym_Freq: CALLS THE PROPER FREQUENCY SYMMETRY ROUTINE*/
bool SYMMETRY::Sym_Freq(int *k, int *l, int *s) {
  if(s!=NULL) {
    for(int i=0;i<DIM;i++) {s[i]=i; s[DIM+i]=1;}
  }
  return((*this.*Sym_ptr)(k, l, s));
}


/*SYMMETRY FUNCTION Sym_NONE: CHECKS ALL SYMMETRY CONDITIONS ON FREQUENCIES*/
bool SYMMETRY::Sym_None(int *k, int *l, int *s) {
  //DECLARE LOCAL VARIABLES
  int flag;
  int ll = 0;
  int n_x = k[dir[0]];
  int n_y = k[dir[1]];
  int n_z = k[dir[2]];


  //PERFORM SYMMETRY OPERATION
  if(n_x==0) {
    if(n_y==0) {
      if(n_z<0) {flag = -1;}
      else if(n_z==0) {flag = 1;}
      else if(n_z==k_test[2]) {flag = 1;}
      else {flag = 2;}
    }
    else if(n_y<0) {
      if(n_z==k_test[dir[2]]) {ll -= n_y;  flag = 1;}
      else {flag = -1;}
    }
    else if(n_y==k_test[1]) {flag = 1;}
    else {flag = 2;}
    if( (l)&&(ll) ) {
      ll += (Lattice->upper[dir[0]]*Lattice->N[dir[1]] + Lattice->upper[dir[1]] )
             *Lattice->N[dir[2]] + Lattice->upper[dir[2]] - Sym_E/UC_DOF;
    }
  }
  else if(n_x<0) {
    if(n_y==k_test[dir[1]]) {
      ll += Lattice->lower[dir[2]]-n_z - 1;
      if(Lattice->upper[dir[2]]==k_test[dir[2]]) {
        ll += Lattice->lower[dir[1]]-n_y;
        ll += (Lattice->lower[dir[0]]-n_x)*(Lattice->N[dir[1]]+Lattice->N[dir[2]]-1);
      }
      else {ll += (Lattice->lower[dir[0]]-n_x)*Lattice->N[dir[2]];}
      flag = 1;
    }
    else if(n_z==k_test[dir[2]]) {
      ll += Lattice->lower[dir[1]]-n_y - 1;
      if(Lattice->upper[dir[1]]==k_test[dir[1]]) {
        ll += (Lattice->lower[dir[0]]-n_x)*(Lattice->N[dir[1]]+Lattice->N[dir[2]]-1);
      }
      else {ll += (Lattice->lower[dir[0]]-n_x)*Lattice->N[dir[1]];}
      flag = 1;
    }
    else {flag = -1;}
  }
  else if(n_x==k_test[dir[0]]) {flag = 1;}
  else if(n_y==k_test[dir[1]]) {flag = 1;}
  else if(n_z==k_test[dir[2]]) {flag = 1;}
  else {flag = 2;}

  if(flag==-1) {
    n_x = -k[dir[0]];
    n_y = -k[dir[1]];
    n_z = -k[dir[2]];
  }
  if(l!=NULL) {
    if(ll) {(*l) = Sym_E/UC_DOF + ll;}
    else {
      (*l) = ( (Lattice->upper[dir[0]]-n_x) * Lattice->N[dir[1]]
               + (Lattice->upper[dir[1]]-n_y) )*Lattice->N[dir[2]]
               + (Lattice->upper[dir[2]]-n_z);
    }
    (*l) *= UC_DOF;
  }

  if(flag<0) {
    if(s!=NULL) {for(ll=0;ll<DIM;ll++) {s[DIM+ll] = -1;}}
    return(true);
  }
  else {return(false);}
}


/*SYMMETRY FUNCTION Sym_2D: CHECKS ALL SYMMETRY CONDITIONS ON FREQUENCIES*/
bool SYMMETRY::Sym_2D(int *k, int *l, int *s) {
  //DECLARE LOCAL VARIABLES
  int j;
  bool flag = false;
  int n[3];
  n[0] = k[0];
  n[1] = k[1];
  n[2] = k[2];
  //PERFORM SYMMETRY OPERATION
  if(n[dir[0]]<0) {n[dir[0]]=-n[dir[0]]; flag=true; if(s) {s[DIM+dir[0]]=-1;}}
  if(n[dir[1]]<0) {n[dir[1]]=-n[dir[1]]; flag=true; if(s) {s[DIM+dir[1]]=-1;}}
  if(n[dir[2]]<0) {n[dir[2]]=-n[dir[2]]; flag=true; if(s) {s[DIM+dir[2]]=-1;}}
  if(n[dir[1]]>n[dir[0]]) {j=n[dir[1]]; n[dir[1]]=n[dir[0]]; n[dir[0]]=j; flag=true; if(s) {j=s[dir[0]]; s[dir[0]]=s[dir[1]]; s[dir[1]]=j;}}
  if(l!=NULL) {
    (*l) = Sym_F - UC_DOF*(1+n[2]);
    if( (dir[0]==0)&&(dir[1]==1) ) {
      (*l) -= UC_DOF*(Lattice->upper[2]+1)*( n[0]*(n[0]+1) )/2;
      (*l) -= UC_DOF*(Lattice->upper[2]+1)*n[1];
    }
    if( (dir[0]==0)&&(dir[1]==2) ) {
      (*l) -= UC_DOF*(Lattice->upper[1]+1)*( n[0]*(n[0]+1) )/2;
      (*l) -= UC_DOF*(n[0]+1)*n[1];
    }
    if( (dir[0]==1)&&(dir[1]==2) ) {
      (*l) -= UC_DOF*( (((Lattice->upper[2]+1)*Lattice->upper[2])/2+Lattice->upper[2]+1)*n[0] );
      (*l) -= UC_DOF*( (n[1])*(n[1]+1) )/2;
    }
  }
  return(flag);
}


/*SYMMETRY FUNCTION Sym_3D: CHECKS ALL SYMMETRY CONDITIONS ON FREQUENCIES*/
bool SYMMETRY::Sym_3D(int *k, int *l, int *s) {
  //DECLARE LOCAL VARIABLES
  int j;
  int *t = new int[DIM];
  for(j=0;j<DIM;j++) {t[j] = j;}
  bool flag = false;
  int x = k[dir[0]];
  int y = k[dir[1]];
  int z = k[dir[2]];
  //PERFORM SYMMETRY OPERATION
  if(x<0) {x=-x; flag=true; if(s) {s[DIM+dir[0]]=-1;}}
  if(y<0) {y=-y; flag=true; if(s) {s[DIM+dir[1]]=-1;}}
  if(z<0) {z=-z; flag=true; if(s) {s[DIM+dir[2]]=-1;}}
  if(y>x) {j=y; y=x; x=j; flag=true; j=t[dir[0]]; t[dir[0]]=t[dir[1]]; t[dir[1]]=j;}
  if(z>y) {j=y; y=z; z=j; flag=true; j=t[dir[1]]; t[dir[1]]=t[dir[2]]; t[dir[2]]=j;}
  if(y>x) {j=y; y=x; x=j; flag=true; j=t[dir[0]]; t[dir[0]]=t[dir[1]]; t[dir[1]]=j;}
  if(l!=NULL) {
    (*l) = Sym_F - UC_DOF;
    (*l) -= UC_DOF*( (x)*(x+1)*(x+2) )/6;
    (*l) -= UC_DOF*( (y)*(y+1) )/2;
    (*l) -= UC_DOF*  (z);
  }
  if(s!=NULL) {
    s[t[dir[0]]] = 0;
    s[t[dir[1]]] = 1;
    s[t[dir[2]]] = 2;
  }
  delete[] t; t=NULL;
  return(flag);
}
