/*                           Harmonic.cpp                           */
/*                            06/05/2009                            */
/*********************************************************************
*    Contains the function definitions for the HARMONIC class.       *
*********************************************************************/


/*DECLARE HEADERS*/
#include "LDCode.h"
#include "Harmonic.h"
#include <cmath>


/*CONSTRUCTOR*/
HARMONIC::HARMONIC() {
  master_flag = true;
  F_flag   = true;
  V_flag   = true;
  E_flag   = true;
  DOS_flag = true;
  thermo   = false;
  return;
}


/*DESTRUCTOR*/
HARMONIC::~HARMONIC() {return;}


/*DEFINE SUBROUTINE Read_Input*/
string HARMONIC::ReadInput(ifstream &Input) {
  //DECLARE LOCAL VARIABLES
  int i;
  string str;
  string Read_Next(ifstream &Input);


  //PARSE THROUGH FILE
  str = Read_Next(Input);
  i = 0;
  while(str[i]) {str[i]=tolower(str[i]);i++;}
  if     (str.compare(0, 1, "true",  1)==0) {master_flag = true;}
  else if(str.compare(0, 1, "false", 1)==0) {master_flag = false;}
  else {Log <<"'HARMONIC' accepts a flag of true or false"<<endl;exit(0);}
  while(!Input.eof()) {
    //identify keyword in HARMONIC category
    str = Read_Next(Input);
    i = 0;
    while(str[i]) {str[i]=tolower(str[i]);i++;}
    if(str.compare(0, 1, "frequency"    , 1)==0) {F_flag   = Boolean(Read_Next(Input));} else
    if(str.compare(0, 1, "velocity"     , 1)==0) {V_flag   = Boolean(Read_Next(Input));} else
    if(str.compare(0, 1, "evect"        , 1)==0) {E_flag   = Boolean(Read_Next(Input));} else
    if(str.compare(0, 3, "dos"          , 3)==0) {DOS_flag = Boolean(Read_Next(Input));} else
    if(str.compare(0, 5, "thermodynamic", 5)==0) {thermo   = Boolean(Read_Next(Input));}
    else{break;}
  }

  return(str);
}


/* */
bool HARMONIC::Boolean(string str) {
  if     (str.compare(0, 1, "true",  1)==0) {return(true);}
  else if(str.compare(0, 1, "false", 1)==0) {return(false);}
  else {Log <<"True or false expected in place of '"<<str<<"'."<<endl;exit(0);}
  return(false);
}


/*HARMONIC FUNCTION Output: WRITES DATA TO LOG FILE*/
void HARMONIC::Output() {
  if(master_flag) {
    Log <<"\nHARMONIC = true"<<endl;
    if( (V_flag)||(E_flag)||(DOS_flag) ) {
      if(!F_flag) {
        Log  <<"  % Flag 'frequency' set to true to fulfill dependencies."<<endl;
        cout <<"  % Flag 'frequency' set to true to fulfill dependencies."<<endl;
        F_flag = true;
      }
    }
    Log <<"  frequency = "<<(F_flag   ? "true" : "false")<<endl;
    Log <<"  velocity  = "<<(V_flag   ? "true" : "false")<<endl;
    Log <<"  evect     = "<<(E_flag   ? "true" : "false")<<endl;
    Log <<"  DOS       = "<<(DOS_flag ? "true" : "false")<<endl;
//    Log <<"  thermodynamic = "<<(thermo ? "true" : "false")<<endl;
  }
  else {Log <<"\nHARMONIC = false"<<endl;}
  return;
}


/*LATTICE FUNCTION Initialize: INITALIZES PARAMETERS AND NON-DIMENSIONALIZES*/
bool HARMONIC::Initialize() {return(master_flag);}


/* */
void HARMONIC::Compute(int &anh_iter, PD_HARMONIC **FC) {
  //EARLY RETURN IF POSSIBLE
  if(!master_flag) {return;}
  if(!F_flag) {return;}
	int dim;				                //Degree of freedom counter
	for(dim=0;dim<DIM;dim++) {
	  if(Lattice->N[dim]<=0) {
		  Log  <<"\n\nSkipping Harmonic Calculation (N"<<dim<<" = "<<Lattice->N[dim]<<").\n"<<endl;
		  cout <<"\n\nSkipping Harmonic Calculation (N"<<dim<<" = "<<Lattice->N[dim]<<").\n"<<endl;
		  return;
	  }
	}


	//DECLARE LOCAL VARIABLES
#if defined PARALLEL
	MPI_Status status;
#endif
	//Misc variables
	int i, j;		                      //Generic counters
	int mntr;				                  //Monitor for output information
	int *UC = Lattice->N;             //Rename number of unit cells in lattice
	int Sym;                          //Number of modes to output
  double dim_time;		              //Time dimensionalizing factor (1/ps)
  double dim_vel;                   //Velocity dimensionalizing factor (m/s)
  double factor;                    //Generic factor
	bool flag;                        //Generic boolean
	//Wave vector
	int *wv_N=new int[DIM];           //Counter over wave vectors
	int *lower=Lattice->lower;        //Lower bound for wave vector
	int *upper=Lattice->upper;        //Upper bound for wave vector
	double *wv=new double[DIM];       //Dimensionless wave vector
	//Values to compute
	double ***D;			                //Dynamical matrix or eigenvectors
	double *Freq;			                //Frequencies
	double **V;				                //Velocities
	//Output files
	ofstream Frequency;	              //Output file for frequencies (THz)
  ofstream Eigenvector;             //Output file for mode shapes ()
  ofstream Vel_out;		              //Output file for velocity (m/s)
	//Subroutines
  void Dynamical_Matrix(PD_HARMONIC **FC, double *wv, double ***D, int dim=-1);//Computes dynamical matrix
	void Eigen1(int UC_DOF, double *Freq, double ***D);//Computes eigenvalues only of D
	void Eigen2(int UC_DOF, double *Freq, double ***D);//Computes eigenvalues/vectors of D
	void Velocity(PD_HARMONIC **FC, double *wv, double *Freq, double **V,
    double ***E);//Computes group velocity (dw/dk) ()
  int LU(int n, double **A, int *permute);//LU decomposition
  double LU_Determinant(int n, double **A, double sign);//Computes determinant of A


	//INITIALIZE VARIABLES AND ALLOCATE MEMORY
	Sym = E_flag ? Symmetry->Sym_E : Symmetry->Sym_F;
	dim_time = 1.0e-12/Parameter->time;//divide by 2Pi to get THz
	dim_vel = Parameter->length/Parameter->time;
	Freq = new double[UC_DOF];	//Frequency of vibration = sqrt(Eigenvalues)
	V = new double*[UC_DOF];		//Velocities
	V[0] = new double[UC_DOF*DIM];
	D = new double** [UC_DOF];	//Dynamical Matrix
	double **A2D = new double*[UC_DOF*UC_DOF];
	double  *A1D = new double [UC_DOF*UC_DOF*2];
	double *A_ptr = A1D;
	for(i=0;i<UC_DOF;i++) {
	  if(i!=UC_DOF-1) {V[i+1] = V[i] + DIM;}
    D[i] = A2D;
    A2D += UC_DOF;
	  for(dim=0;dim<UC_DOF;dim++) {D[i][dim] = A_ptr; A_ptr += 2;}
	}


	//OPEN OUTPUT FILES AND PRINT HEADERS
	if(PE==0) {
    Log  <<"\n\nBeginning Harmonic Calculation of ";
    cout <<"\n\nBeginning Harmonic Calculation of ";
    if(F_flag) {
      Log  <<"Frequency";
      cout <<"Frequency";
      Frequency.open("Frequency.txt");
      Frequency <<Symmetry->Sym_F<<" Out of "<<Symmetry->TotalDOF
        <<" Normal Mode Frequencies in 1/ps for a "<<UC[0]<<"x"<<UC[1]<<"x"<<UC[2];
      for(i=0;i<Material->n_mat;i++) {Frequency <<" "<<Material->symb[i];}
      Frequency <<" Lattice."<<endl;
      Frequency.precision(6);
      fixed(Frequency);
    }
    if(V_flag) {
      if(F_flag) {Log <<", ";  cout <<", ";}
      Log  <<"Velocity";
      cout <<"Velocity";
      Vel_out.open("Velocity.txt");
      Vel_out <<Symmetry->Sym_F<<" Out of "<<Symmetry->TotalDOF
        <<" Velocities in m/s for a "<<UC[0]<<"x"<<UC[1]<<"x"<<UC[2];
      for(i=0;i<Material->n_mat;i++) {Vel_out <<" "<<Material->symb[i];}
      Vel_out <<" Lattice.\n";
      Vel_out.precision(4);
      fixed(Vel_out);
    }
    if(E_flag) {
      if(F_flag||V_flag) {Log <<", ";  cout <<", ";}
      Log  <<"Eigenvector";
      cout <<"Eigenvector";
      Eigenvector.open("Eigenvector.txt");
      Eigenvector <<Symmetry->Sym_E<<" Out of "<<Symmetry->TotalDOF
        <<" Normal Modes for a "<<UC[0]<<"x"<<UC[1]<<"x"<<UC[2];
      for(i=0;i<Material->n_mat;i++) {Eigenvector <<" "<<Material->symb[i];}
      Eigenvector <<" Lattice.\n";
      Eigenvector.precision(10);
      fixed(Eigenvector);
    }
    Log  <<".\n"<<endl;
    cout <<".\n"<<endl;
	}


	//CALCULATE DENSITY OF STATES
	i = 0;
	mntr = 1;
  int loop = 0;
	for(wv_N[0]=upper[0];wv_N[0]>=lower[0];wv_N[0]--) {
    wv[0] = double(wv_N[0])/double(UC[0]);
  for(wv_N[1]=upper[1];wv_N[1]>=lower[1];wv_N[1]--) {
    wv[1] = double(wv_N[1])/double(UC[1]);
  for(wv_N[2]=upper[2];wv_N[2]>=lower[2];wv_N[2]--) {
    wv[2] = double(wv_N[2])/double(UC[2]);

    //Check symmetry conditions
    if(Symmetry->Sym_Negative(wv_N)<0.0) {continue;}//Continue if k=-k
    flag = Symmetry->Sym_Freq(wv_N);//Defined symmetry operations
    if( (!E_flag)&&flag ) {continue;}

    //Output data for each phonon mode//
    //Compute frequencies and eigenvectors
    if(loop!=PE) {
      if(PE!=0) {
        loop += 1;
        continue;
      }
    }
    else {
      loop -= nPE;
      Dynamical_Matrix(FC, wv, D);//Compute Dynamical Matrix
    }

    if(E_flag) {
      if( (PE!=0)||(loop==-nPE) ) {Eigen2(UC_DOF, Freq, D);}//Compute Eigenvectors and Eigenvalues
#if defined PARALLEL
      if(PE!=0) {MPI_Ssend(D, UC_DOF*UC_DOF*2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);}
      else if(loop!=-nPE) {MPI_Recv(D, UC_DOF*UC_DOF*2, MPI_DOUBLE, loop+nPE, 0, MPI_COMM_WORLD, &status);}
#endif
      if(PE==0) {
        for(j=0;j<UC_DOF;j++) {
          for(dim=0;dim<UC_DOF;dim++) {
            Eigenvector.width(15);
            Eigenvector <<D[dim][j][0]<<" ";
            Eigenvector.width(15);
            Eigenvector <<D[dim][j][1]<<" ";
          }
          Eigenvector <<endl;
        }
      }
    }
    else {if( (PE!=0)||(loop==-nPE) ) {Eigen1(UC_DOF, Freq, D);}}
    //Output frequencies
    if(F_flag&&!flag) {
      //Check for zero or imaginary frequencies
      if( (wv_N[0]==0)&&(wv_N[1]==0)&&(wv_N[2]==0) ) {
        for(j=UC_DOF-1;j>=UC_DOF-DIM-1;j--) {if(fabs(Freq[j])<1.0e-7) {Freq[j]=0.0;}}//Check last four for possible torsion mode
      }
#if defined PARALLEL
      if(PE!=0) {MPI_Ssend(Freq, UC_DOF, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);}
      else if(loop!=-nPE) {MPI_Recv(Freq, UC_DOF, MPI_DOUBLE, loop+nPE, 1, MPI_COMM_WORLD, &status);}
#endif
      if(PE==0) {
        //Update List
        for(j=0;j<UC_DOF;j++) {
          if(Freq[j]<0.0) {
            Log  <<"Imaginary frequencies detected.  No anharmonic calculation possible."<<endl;
            cout <<"Imaginary frequencies detected.  No anharmonic calculation possible."<<endl;
            anh_iter = 0;
          }
          else {Freq[j] = sqrt(Freq[j])*dim_time;}
          Frequency.width(12);
          Frequency <<Freq[j];
        }
        Frequency <<endl;
      }
    }
    //Output velocities
    if(V_flag&&!flag) {
      Velocity(FC, wv, Freq, V, D);//Compute velocities
#if defined PARALLEL
      if(PE!=0) {MPI_Ssend(V, UC_DOF, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);}
      else if(loop!=-nPE) {MPI_Recv(V, UC_DOF, MPI_DOUBLE, loop+nPE, 2, MPI_COMM_WORLD, &status);}
#endif
      if(PE==0) {
        for(j=0;j<UC_DOF;j++) {
          for(dim=0;dim<DIM;dim++) {Vel_out.width(12);Vel_out <<V[j][dim]*dim_vel;}
        }
        Vel_out <<endl;
      }
    }
    //Output progress
    if(PE==0) {
      i += UC_DOF;
      if(Sym*mntr<=100*i) {
        Log <<"Harmonic Calculation is "<<int(double(i)/Sym*100)<<" % Complete.";
        Log <<"  "<<i<<" of "<<Sym<<" waves calculated."<<endl;
        cout <<int(double(i)/Sym*100)<<"% done.";
        cout <<"  "<<i<<" of "<<Sym<<" waves calculated."<<endl;
        while(Sym*mntr<=100*i) {mntr += 1;}
      }
    }

    loop += 1;

  }  //for(wv_n[2])
  }  //for(wv_n[1])
	}  //for(wv_n[0])

#if defined PARALLEL
  if(PE==0) {
    for(i=1;i<nPE;i++) {
      MPI_Ssend(&anh_iter, 1, MPI_INT, i, 10, MPI_COMM_WORLD);
    }
#endif
    if(F_flag) {Frequency.close();}
    if(V_flag) {Vel_out.close();}
    if(E_flag) {Eigenvector.close();}
#if defined PARALLEL
	}
	else {
	  MPI_Recv(&anh_iter, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &status);
	}
#endif


	//DEALLOCATE MEMORY
	delete[] wv_N;    wv_N   =NULL;
	delete[] wv;      wv     =NULL;
	delete[] Freq;    Freq   =NULL;
	delete[] V;       V      =NULL;
	delete[] D[0][0]; D[0][0]=NULL;
	delete[] D[0];    D[0]   =NULL;
	delete[] D;       D      =NULL;

	return;
}


/*COMPUTE THE PHONON DENSITY OF STATES*/
void HARMONIC::DOS(void) {
  if(PE!=0) {return;}
  if(!DOS_flag) {return;}
  Log  <<"\nComputing phonon density of states."<<endl;
  cout <<"\nComputing phonon density of states."<<endl;
	//DELCARE LOCAL VARIABLES
	int i, j, jj, mntr;	  //Counters
	int *wv_N=new int[DIM];			//Wave vector counters
	int MaxSym = Symmetry->Sym_F;//Phonons at each wave vector, total number of phonons
	int *lower = Lattice->lower, *upper = Lattice->upper;	//Lower and upper bounds of 1st Brillioun zone
	int bins = 1000;			      //Number of bins to use for histogram
	double factor;				      //Multiplicity factor for symmetry conditions
	double F_max;				        //Maximum frequency or bin size
	double *DOS;		            //Density of states


	//INITIALIZE FILES AND DATA
	ifstream Frequency("Frequency.txt");
	if(!Frequency.is_open()) {Log <<"Could not open file '"<<"Frequency.txt"<<"' for DOS calculation."<<endl;exit(0);}
	Frequency.ignore(10000, '\n');
	double *QH_F = new double[MaxSym];
	F_max = 0;
	for(i=0;i<MaxSym;i++) {
	  Frequency >>QH_F[i];
	  if(F_max<QH_F[i]) {F_max = QH_F[i];}
	}
	Frequency.close();
	bins = Symmetry->TotalDOF/(DIM*2000);//6000 phonons per bin
	if(bins<50) {bins = 50;}	    //Minimum of 50 bins
	F_max = F_max/double(bins-1); //Bin size
	DOS = new double [bins+1];
	for(i=0;i<=bins;i++) {DOS[i] = 0.0;}


	//MAKE HISTOGRAM FROM DATA
	i = 0;
	mntr = 1;
	for(wv_N[0]=upper[0];wv_N[0]>=lower[0];wv_N[0]--) {
  for(wv_N[1]=upper[1];wv_N[1]>=lower[1];wv_N[1]--) {
  for(wv_N[2]=upper[2];wv_N[2]>=lower[2];wv_N[2]--) {

    //Perform symmetry operations
    if((factor=Symmetry->Sym_Negative(wv_N))<0.0) {continue;}
    Symmetry->Sym_Freq(wv_N, &i);

    //Accumulate The DOS
    for(j=0;j<UC_DOF;j++) {
      jj = 0;
      while(jj <= bins) {
        if( QH_F[i]<=(F_max*(jj+0.5)) ) {
          DOS[jj] += factor;
          break;
        }
        jj += 1;
      }
      i += 1;

      //Output Progress To Screen
      if(MaxSym*mntr<=10*i) {
        cout <<int(double(i)/MaxSym*10)*10<<"% done.";
        cout <<"  "<<i<<" of "<<MaxSym<<" waves calculated."<<endl;
        while(MaxSym*mntr<=10*i) {mntr += 1;}
      }

    }  //for(j)
  }  //for(wv_N[2])
  }  //for(wv_N[1])
	}  //for(wv_N[0])


	//PRINT RESULTS AND CLOSE FILES
	DOS[0] = DOS[0] - 3; //Remove 3 translational modes (k = [0 0 0])
//	factor = 1.0e-18*(1.0e-12*1.0e30/F_max/V);  //scale to a reasonable range (1e12 Phonons/m^3)
	factor = 1.0/(Lattice->V*DIM*F_max);  //scale to integrate to 10^12
	for(i=0;i<DIM;i++) {factor /= Lattice->N[i]*Parameter->length*1.0e9;}
	ofstream Output("DOS.xls");
	Output <<"Frequency (1/ps)"<<'\t'<<"Number Density per Volume (Phonons*ps/nm^3)\tDOS integrates to 3*n/V_uc"<<endl;
	for(i=0;i<=bins;i++) {
		Output <<F_max*(i+0.5)<<'\t'<<DOS[i]*factor<<endl;
	}
	Output.close();


	//DEALLOCATE MEMORY AND RETURN
	delete[] wv_N;       wv_N      =NULL;
	delete[] QH_F;       QH_F      =NULL;
	delete[] DOS;        DOS       =NULL;

	return;
}
