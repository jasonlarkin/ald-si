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
  thermo = false;
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
    if(str.compare(0, 5, "thermodynamic", 5)==0) {thermo = Boolean(str);}
    else{break;}
  }

  return(str);
}


/* */
bool HARMONIC::Boolean(string str) {
  if     (str.compare(0, 1, "true",  1)==0) {return(true);}
  else if(str.compare(0, 1, "false", 1)==0) {return(false);}
  else {Log <<"True or false expected."<<endl;exit(0);}
  return(false);
}


/*HARMONIC FUNCTION Output: WRITES DATA TO LOG FILE*/
void HARMONIC::Output() {
  if(master_flag) {
    Log <<"\nHARMONIC = true"<<endl;
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
	int dim;				                //Degree of freedom counter
	for(dim=0;dim<DIM;dim++) {
	  if(Lattice->N[dim]<=0) {
		  Log  <<"\n\nSkipping Density of States Calculation (N"<<dim<<" = "<<Lattice->N[dim]<<").\n"<<endl;
		  cout <<"\n\nSkipping Density of States Calculation (N"<<dim<<" = "<<Lattice->N[dim]<<").\n"<<endl;
		  return;
	  }
	}


	//DECLARE LOCAL VARIABLES
	//Misc variables
	int i, j;		                      //Generic counters
	int mntr;				                  //Monitor for output information
	int *UC = Lattice->N;             //Rename number of unit cells in lattice
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
	//Density of states
	int bins;				                  //Number of bins for DOS histogram
	double MaxFreq;			              //Maximum frequency in spectrum
	double *List;			                //List of all the frequencies for DOS
	double **Density;		              //Density of states histogram
	//Subroutines
  void Dynamical_Matrix(PD_HARMONIC **FC, double *wv, double ***D, int dim=-1);//Computes dynamical matrix
	void Eigen1(int UC_DOF, double *Freq, double ***D);//Computes eigenvalues only of D
	void Eigen2(int UC_DOF, double *Freq, double ***D);//Computes eigenvalues/vectors of D
	void Velocity(PD_HARMONIC **FC, double *wv, double *Freq, double **V,
    double ***E);//Computes group velocity (dw/dk) ()
  int LU(int n, double **A, int *permute);//LU decomposition
  double LU_Determinant(int n, double **A, double sign);//Computes determinant of A


	//INITIALIZE VARIABLES
	dim_time = 1.0e-12/Parameter->time;//divide by 2Pi to get THz
	dim_vel = Parameter->length/Parameter->time;
	bins = Symmetry->TotalDOF/(DIM*6000);//6000 waves per bin
	if(bins<50) {bins = 50;}	           //Minimum of 50 bins


	//ALLOCATE MEMORY
	Freq = new double[UC_DOF];	//Frequency of vibration = sqrt(Eigenvalues)
	List = new double[Symmetry->Sym_E];	  //All frequencies
	V = new double*[UC_DOF];		//Velocities
	V[0] = new double[UC_DOF*DIM];
	D = new double** [UC_DOF];	//Dynamical Matrix
	D[0] = new double* [UC_DOF*UC_DOF];
	for(i=0;i<UC_DOF;i++) {
	  if(i!=UC_DOF-1) {
		V[i+1] = V[i] + DIM;
		D[i+1] = D[i] + UC_DOF;
	  }
	  D[i][0] = new double [UC_DOF*2];
	  for(dim=0;dim<UC_DOF-1;dim++) {D[i][dim+1] = D[i][dim] + 2;}
	}


	//OPEN OUTPUT FILES AND PRINT HEADERS
  Log <<"\n\nBeginning Density of States and Velocity Calculation. (Using "<<bins<<" bins)\n\n";
  cout <<"\n\nBeginning Density of States and Velocity Calculation. (Using "<<bins<<" bins)\n\n";
  ofstream Frequency("Frequency.txt");	//Output file for frequencies (THz)
  ofstream Eigenvector("Eigenvector.txt");//Output file for mode shapes ()
  ofstream Vel_out("Velocity.txt");		//Output file for velocity (m/s)
  Frequency <<Symmetry->Sym_F<<" Out of "<<Symmetry->TotalDOF
    <<" Normal Mode Frequencies in 1/ps for a "<<UC[0]<<"x"<<UC[1]<<"x"<<UC[2];
  for(i=0;i<Material->n_mat;i++) {Frequency <<" "<<Material->symb[i];}
  Frequency <<" Lattice.\n";
  Eigenvector <<Symmetry->Sym_E<<" Out of "<<Symmetry->TotalDOF
    <<" Normal Modes for a "<<UC[0]<<"x"<<UC[1]<<"x"<<UC[2];
  for(i=0;i<Material->n_mat;i++) {Eigenvector <<" "<<Material->symb[i];}
  Eigenvector <<" Lattice.\n";
  Vel_out <<Symmetry->Sym_F<<" Out of "<<Symmetry->TotalDOF
    <<" Velocities in m/s for a "<<UC[0]<<"x"<<UC[1]<<"x"<<UC[2];
  for(i=0;i<Material->n_mat;i++) {Vel_out <<" "<<Material->symb[i];}
  Vel_out <<" Lattice.\n";
  Frequency.precision(6);
  fixed(Frequency);
  Vel_out.precision(4);
  fixed(Vel_out);
  Eigenvector.precision(10);
  fixed(Eigenvector);


	//CALCULATE DENSITY OF STATES
	MaxFreq = 0.0;
	i = 0;
	mntr = 1;
	for(wv_N[0]=upper[0];wv_N[0]>=lower[0];wv_N[0]--) {
	  wv[0] = double(wv_N[0])/double(UC[0]);
	  for(wv_N[1]=upper[1];wv_N[1]>=lower[1];wv_N[1]--) {
      wv[1] = double(wv_N[1])/double(UC[1]);
      for(wv_N[2]=upper[2];wv_N[2]>=lower[2];wv_N[2]--) {
        wv[2] = double(wv_N[2])/double(UC[2]);

        //Check symmetry conditoins
        if(Symmetry->Sym_Negative(wv_N)<0.0) {continue;}//Continue if k=-k
        flag = Symmetry->Sym_Freq(wv_N);//Defined symmetry operations

        //Output data for each phonon mode//
        //Compute frequencies and eigenvectors
        Dynamical_Matrix(FC, wv, D);//Compute Dynamical Matrix

        Eigen2(UC_DOF, Freq, D);//Compute Eigenvectors and Eigenvalues
        //Check for zero or imaginary frequencies
        if( (wv_N[0]==0)&&(wv_N[1]==0)&&(wv_N[2]==0) ) {
          for(j=UC_DOF-1;j>=UC_DOF-DIM-1;j--) {if(fabs(Freq[j])<1.0e-7) {Freq[j]=0.0;}}//Check last four for possible torsion mode
        }
        //Update List and output eigenvectors
        for(j=0;j<UC_DOF;j++,i++) {
          if(Freq[j]<0.0) {
            Log  <<"Imaginary frequencies detected.  No anharmonic calculation possible."<<endl;
            cout <<"Imaginary frequencies detected.  No anharmonic calculation possible."<<endl;
            anh_iter = 0;
          }
          else {Freq[j] = sqrt(Freq[j])*dim_time;}
          List[i] = Freq[j];
          for(dim=0;dim<UC_DOF;dim++) {
            Eigenvector.width(15);
            Eigenvector <<D[dim][j][0]<<" ";
            Eigenvector.width(15);
            Eigenvector <<D[dim][j][1]<<" ";
          }
          Eigenvector <<endl;
        }
        //Output progress
        if(Symmetry->Sym_E*mntr<=100*i) {
          Log <<"DOS Calculation is "<<int(double(i)/Symmetry->Sym_E*100)<<" % Complete.";
          Log <<"  "<<i<<" of "<<Symmetry->Sym_E<<" waves calculated."<<endl;
          cout <<int(double(i)/Symmetry->Sym_E*100)<<"% done.";
          cout <<"  "<<i<<" of "<<Symmetry->Sym_E<<" waves calculated."<<endl;
          while(Symmetry->Sym_E*mntr<=100*i) {mntr += 1;}
        }
        //Skip frequency and velocity if phonon is not unique
        if(flag) {continue;}
        //Output frequencies
        for(j=0;j<UC_DOF;j++) {
          Frequency.width(12);
          Frequency <<Freq[j];
          if(MaxFreq<Freq[j]) {MaxFreq = Freq[j];}
        }
        Frequency <<endl;
        //Output velocities
        Velocity(FC, wv, Freq, V, D);//Compute velocities
        for(j=0;j<UC_DOF;j++) {
          for(dim=0;dim<DIM;dim++) {Vel_out.width(12);Vel_out <<V[j][dim]*dim_vel;}
        }
        Vel_out <<endl;

      }  //for(wv_n[2])
	  }  //for(wv_n[1])
	}  //for(wv_n[0])
	Frequency.close();
	Eigenvector.close();
	Vel_out.close();


	//BUILD HISTOGRAM
	factor = double(MaxFreq)/double(bins);    //Generic scaling factor
	Density = new double* [bins+1];
	Density[0] = new double [(bins+1)*2];
	for(i=0;i<=bins;i++) {
		if(i!=bins) {Density[i+1] = Density[i] + 2;}
		Density[i][0] = factor * (i + 0.5);
		Density[i][1] = 0.0;
	}
	i=0;
	for(wv_N[0]=upper[0];wv_N[0]>=lower[0];wv_N[0]--) {
	  wv[0] = double(wv_N[0])/double(UC[0]);
	  for(wv_N[1]=upper[1];wv_N[1]>=lower[1];wv_N[1]--) {
      wv[1] = double(wv_N[1])/double(UC[1]);
      for(wv_N[2]=upper[2];wv_N[2]>=lower[2];wv_N[2]--) {
        wv[2] = double(wv_N[2])/double(UC[2]);

        //Check Conditions
        if( (mntr=int(Symmetry->Sym_Negative(wv_N)+0.001))<=0 ) {continue;}

        for(dim=0;dim<UC_DOF;dim++) {
          j = 0;
          while(j <= bins) {
            if( List[i+dim]<=(factor*(j+1.0)) ) {
              Density[j][1] += mntr;
              break;
            }
            j += 1;
          }
        }
        i += UC_DOF;

      }
	  }
	}


	//SCALE DATA - int(n(w)dw) = 3(N-1)/V = 12/a^3 * (number of atoms per FCC site)
	Density[0][1] = Density[0][1] - 3; //Remove 3 translational modes (k = [0 0 0])
	factor = 1.0e-12/(2.0*Density[0][0]);
	factor /= Lattice->V;
	for(i=0;i<DIM;i++) {factor /= Lattice->N[i]*Parameter->length;}
	factor *= 1.0e-18;  //scale to a reasonable range (1e18 Phonons/m^3)


	//OUTPUT DOS HISTOGRAM
	ofstream Out_DOS("DOS.xls");
	Out_DOS <<"Frequency (1/ps)"<<'\t'<<"Number Density per Volume (1e18 Phonons/m^3)"<<endl;
	for(i=0;i<=bins;i++) {
		Out_DOS <<Density[i][0]<<'\t'<<Density[i][1]*factor<<endl;
	}
	Out_DOS.close();


	//DEALLOCATE MEMORY
	delete[] wv_N;  wv_N = NULL;
	delete[] wv;  wv = NULL;
	delete[] Freq;  Freq = NULL;
	delete[] V;  V = NULL;
	for(i=0;i<UC_DOF;i++) {delete[] D[i][0];  D[i][0] = NULL;}
	delete[] D[0];  D[0] = NULL;
	delete[] D;  D = NULL;
	delete[] List;  List = NULL;
	delete[] Density[0];  Density[0] = NULL;
	delete[] Density;  Density = NULL;

	return;
}
