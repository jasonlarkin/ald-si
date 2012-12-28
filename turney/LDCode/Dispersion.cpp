/*                          Dispersion.cpp                          */
/*                            12/02/2008                            */
/*********************************************************************
*    Source code file for functions in the DISPERSION class.         *
*********************************************************************/

/*DECLARE HEADERS*/
#include "LDCode.h"
#include "Dispersion.h"
#include <cmath>


/*DEFINE PREPROCESSOR VARIABLES*/


DISPERSION::DISPERSION() {
  points = -1;                //Default do not compute dispersion
  coordinate = CARTESIAN;     //Default cartesian coordinate system
  begin = end = NULL;
  return;
}


DISPERSION::~DISPERSION() {
  delete[] begin;
  begin = NULL;
  delete[] end;
  end = NULL;
}


string DISPERSION::Define(ifstream &Input) {
  //DECLARE LOCAL VARIABLES
  int i;
  string str;
  string Read_Next(ifstream &Input);
  void Coordinate(int coord, double *x, ifstream &Input);

  //PARSE THROUGH FILE
  str = Read_Next(Input);
  points = int(atof(str.c_str()));
  if(begin==NULL) {
    begin = new double[DIM];
    end = new double[DIM];
  }
  while(!Input.eof()) {
    //identify keyword in DISPERSION category
    str = Read_Next(Input);
    i = 0;while(str[i]) {str[i]=tolower(str[i]);i++;}
    if(str.compare(0, 5, "coordinate", 5)==0) {coordinate = KeywordToInt(Input);}
    else if(str.compare(0, 3, "begin", 3)==0) {Coordinate(coordinate, begin, Input);}
    else if(str.compare(0, 3, "end",   3)==0) {Coordinate(coordinate, end, Input);}
    else {break;}
  }

  return(str);
}


/*DISPERSION FUNCTION KeywordToInt: IDENTIFIES A KEYWORD USED IN THE DISPERSION*/
int DISPERSION::KeywordToInt(ifstream &Input) {
  int i=0;
  string str;
  string Read_Next(ifstream &Input);
  str = Read_Next(Input);
  while(str[i]) {str[i]=tolower(str[i]);i++;}
  if(str.compare(0, 3, "cartesian",  3)==0) {return(CARTESIAN);}//x, y, z
  if(str.compare(0, 3, "cylindrical",3)==0) {return(CYLINDRICAL);}//r, theta, z
  if(str.compare(0, 3, "spherical",  3)==0) {return(SPHERICAL);}//r, theta (from x axis), phi (from z axis)
  if(str.compare(0, 3, "direct",     3)==0) {return(DIRECT);}//a1, a2, a3
  if(str.compare(0, 3, "reciprocal", 3)==0) {return(RECIPROCAL);}//b1, b2, b3
  Log <<"Error. Unknown keyword for DISPERSION: "<<str<<endl;
  exit(0);
}


/*SUBROUTINE Compute: COMPUTES DISPERSION RELATION*/
void DISPERSION::Compute(PD_HARMONIC **FC) {
	//EARLY RETURN IF POSSIBLE
	if(points<0) {return;}


	//DECLARE LOCAL VARIABLES
	int i, dim, n, mntr;		//Degree of freedom and wave vector counters
	double wv_mag;					//Wave vector magnitude
	double dim_time = 1.0e-12/(Parameter->time);
	ofstream Out_Dispersion;
	void Dynamical_Matrix(PD_HARMONIC **FC, double *wv, double ***D, int dim=-1);
	void Eigen1(int UC_DOF, double *Freq, double ***D);//Computes eigenvalues only of D
	void Eigen2(int UC_DOF, double *Freq, double ***D);//Computes eigenvalues/vectors of D


	//ALLOCATE MEMORY
	double *wv = new double[DIM];				//Psuedo-Dimensional wave vector
	double *Freq = new double[UC_DOF];    //Frequency of vibration = sqrt(Eigenvalues)
	double ***D = new double** [UC_DOF];  //Dynamical matrix
	D[0] = new double* [UC_DOF*UC_DOF];
	for(i=0;i<UC_DOF;i++) {
		if(i!=UC_DOF-1) {D[i+1] = D[i] + UC_DOF;}
		D[i][0] = new double [UC_DOF*2];
		for(dim=0;dim<UC_DOF-1;dim++) {D[i][dim+1] = D[i][dim] + 2;}
	}


	//OUTPUT INITIAL DATA TO 'Out_Dispersion'
#if defined PARALLEL
	if(PE==0) {
#endif
    Out_Dispersion.open("Dispersion.xls");
    Out_Dispersion <<"Wave Vector Magnitude ()"<<'\t';
    for(i=0;i<UC_DOF;i++) {Out_Dispersion <<"Frequency (branch "<<i<<") (1/ps)"<<'\t';}
    Out_Dispersion <<endl;
    Log <<"\n\nBegining Dispersion Relation Calculation.\n"<<endl;
    cout <<"\n\nBegining Dispersion Relation Calculation.\n"<<endl;
#if defined PARALLEL
	}
#endif


	//TRANSFORM begin AND end TO RECIPROCAL LATTICE COORDINATES
	double *Beg = new double[DIM];
	double *End = new double[DIM];
	for(i=0;i<DIM;i++) {
	  Beg[i] = 0.0;
	  End[i] = 0.0;
	  for(dim=0;dim<DIM;dim++) {
	    Beg[i] += Lattice->a[i][dim]*begin[dim];
	    End[i] += Lattice->a[i][dim]*end[dim];
	  }
	  if(fabs(Beg[i])<1.0e-8) {Beg[i] = 0.0;}
	  if(fabs(End[i])<1.0e-8) {End[i] = 0.0;}
	}


	//CYCLE THROUGH WAVE VECTORS
	mntr=1;
	for(n=PE;n<=points;n+=nPE) {
		wv_mag = 0.0;
		for(dim=0;dim<DIM;dim++) {
		  if(points) {
        wv[dim] = Beg[dim]+(End[dim]-Beg[dim])*double(n)/double(points);
		  }
		  else {wv[dim] = Beg[dim];}
			wv_mag += wv[dim]*wv[dim];
		}
		wv_mag = sqrt(wv_mag);

		Dynamical_Matrix(FC, wv, D);//Compute Dynamical Matrix

/*		Log <<"DYNAMICAL MATRIX"<<endl;
		for(i=0;i<UC_DOF;i++) {	//Output Dynamical Matrix
		  Log <<"D"<<i<<" = [";
			for(dim=0;dim<UC_DOF;dim++) {Log <<D[i][dim][0]<<" + "<<D[i][dim][1]<<"*i,   ";}
			Log <<"]"<<endl;
		}// */

		Eigen1(UC_DOF, Freq, D);	//Compute Eigenvalues (1) and optionally Eigenvectors (2)

/*		Log <<"EIGENVALUES"<<endl;//Be sure to change to Eigen2
		for(i=0;i<UC_DOF;i++) {Log <<"lambda"<<i<<" = "<<Freq[i]<<endl;}//Output Frequencies and Eigenvectors
		Log <<"EIGENVECTORS"<<endl;
		for(i=0;i<UC_DOF;i++) {
		  Log <<"E"<<i<<" = [";
			for(dim=0;dim<UC_DOF;dim++) {Log <<D[dim][i][0]<<" + "<<D[dim][i][1]<<"*i,   ";}
			Log <<"]"<<endl;
		}// */

#if defined PARALLEL
		if(PE==0) {
#endif
      //Compute and output frequencies
      Out_Dispersion <<wv_mag<<'\t';
      for(i=0;i<UC_DOF;i++) {
        if(Freq[i]>=0.0) {Out_Dispersion <<sqrt(Freq[i])*dim_time<<'\t';}
        else {Out_Dispersion <<-sqrt(fabs(Freq[i]))*dim_time<<'\t';}
      }
      Out_Dispersion <<endl;


#if defined PARALLEL
      //Gather and output data
      MPI_Status status;
      for(dim=1;dim<nPE;dim++) {
        if(n+dim>points) {break;}
        MPI_Recv(&wv_mag,   1, MPI_DOUBLE, dim, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(wv,      DIM, MPI_DOUBLE, dim, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(Freq, UC_DOF, MPI_DOUBLE, dim, 2, MPI_COMM_WORLD, &status);
        Out_Dispersion <<wv_mag<<'\t';
        for(i=0;i<UC_DOF;i++) {
          if(Freq[i]>=0.0) {Out_Dispersion <<sqrt(Freq[i])*dim_time<<'\t';}
          else {Out_Dispersion <<-sqrt(fabs(Freq[i]))*dim_time<<'\t';}
        }
        Out_Dispersion <<endl;
      }
#endif

      //Output progress
      i = n + nPE - 1;
      if(i>points) {i=points;}
      if(points==0) {
        Log <<"Dispersion is 100% complete. k* = < ";
        for(dim=0;dim<DIM;dim++) {Log <<wv[dim]<<" ";}
        Log <<">, |k*| = "<<wv_mag<<endl;
        cout <<"Dispersion is 100% complete. k* = < ";
        for(dim=0;dim<DIM;dim++) {cout <<wv[dim]<<" ";}
        cout <<">, |k*| = "<<wv_mag<<endl;
      } else
      if(points*mntr<=100*i) {
        Log <<"Dispersion is "<<int(double(i)/points*100)<<"% complete. k* = < ";
        for(dim=0;dim<DIM;dim++) {Log <<wv[dim]<<" ";}
        Log <<">, |k*| = "<<wv_mag<<endl;
        cout <<"Dispersion is "<<int(double(i)/points*100)<<"% complete. k* = < ";
        for(dim=0;dim<DIM;dim++) {cout <<wv[dim]<<" ";}
        cout <<">, |k*| = "<<wv_mag<<endl;
        while(points*mntr<=100*i) {mntr += 1;}
      }
#if defined PARALLEL
    }
    else {
      //Send data to master PE
      MPI_Ssend(&wv_mag,   1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Ssend(wv,      DIM, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      MPI_Ssend(Freq, UC_DOF, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }
#endif
	}



	//DEALLOCATE MEMORY
  if(PE==0) {Out_Dispersion.close();}
	delete[] wv;  wv = NULL;
	delete[] Freq;  Freq = NULL;
	for(i=0;i<UC_DOF;i++) {delete[] D[i][0];  D[i][0] = NULL;}
	delete[] D[0];  D[0] = NULL;
	delete[] D;  D = NULL;
	delete[] Beg;  Beg = NULL;
	delete[] End;  End = NULL;

	return;
}


/*DISPERSION FUNCTION Output: WRITES DATA TO LOG FILE*/
void DISPERSION::Output() {
  int i;
  Log <<"\nDISPERSION = "<<points<<endl;
  Log <<"  Begin =";
  for(i=0;i<DIM;i++) {Log <<" "<<begin[i];}
  Log <<"\n  End =";
  for(i=0;i<DIM;i++) {Log <<" "<<end[i];}
  Log <<endl;
  return;
}


/*DISPERSION FUNCTION Initialize: INITALIZES PARAMETERS AND NON-DIMENSIONALIZES*/
bool DISPERSION::Initialize() {
  if(points<0) {return(false);}
  int i;
  for(i=0;i<DIM;i++) {begin[i] *= Parameter->length;}
  for(i=0;i<DIM;i++) {end[i] *= Parameter->length;}
  return(true);
}
