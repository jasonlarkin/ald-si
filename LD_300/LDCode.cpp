/*						                LDCode.cpp		                  			*/
/*							              12/30/2008	              						*/
/*********************************************************************
*    This program performs lattice dynamics related calculations on  *
*  a crystal structure and outputs information as directed.	    	   *
*********************************************************************/

/*DEFINE HEADERS*/
#include "LDCode.h"
#include "PreProcess.h"
#include "Dispersion.h"
#include "Harmonic.h"
#include "Anharmonic.h"
#include "PostProcess.h"
#include "Neighbor_Force.h"
#include <ctime>


/*DECLARE GLOBAL VARIABLES*/
int PE, nPE;
ofstream Log1;                //Log file
//int DIM;                      //Dimension of system
int UC_DOF;                   //Degrees of freedom = DIM*Unit_Cell->n_atom
int BODY = 0;                 //Highest order of atomic interactions
MATERIAL  *Material;          //MATERIAL does not depend upon other classes
LATTICE   *Lattice;           //LATTICE does not depend upon other classes
UNIT_CELL *Unit_Cell;         //UNIT_CELL must be after LATTICE, MATERIAL
SYMMETRY  *Symmetry;          //SYMMETRY must be after LATTICE
PARAMETER *Parameter;				  //(Non-)dimensionalizing parameters


/*MAIN ROUTINE*/
int main(int argc, char **argv) {
  //DECLARE LOCAL VARIABLES
  int Begin_time = time(0);	  //Begin time
  int rtrn = 0;
  PREPROCESS  *Preprocess  = new PREPROCESS();
  POTENTIAL   **Potential  = NULL;
  DISPERSION  *Dispersion  = new DISPERSION();
  HARMONIC    *Harmonic    = new HARMONIC();
  ANHARMONIC  *Anharmonic  = new ANHARMONIC();
  POSTPROCESS *Postprocess = new POSTPROCESS();
#if defined PARALLEL
	//BEGIN PARALLEL CODE
	MPI_Init(&argc, &argv);
	MPI_Comm_size( MPI_COMM_WORLD, &nPE );
	MPI_Comm_rank( MPI_COMM_WORLD, &PE );
	if(PE==0) {
	  MPI_Status status;
	  Log1.open("LDCode.log");
	  Log1 <<"Parallel job using "<<nPE<<" processing elements."<<endl;
	  Log1 <<"Processing element "<<PE<<" (Master) reports."<<endl;
	  cout <<"Parallel job using "<<nPE<<" processing elements."<<endl;
	  cout <<"Processing element "<<PE<<" (Master) reports."<<endl;
	  int node;
	  for(rtrn=1;rtrn<nPE;rtrn++) {
		MPI_Recv(&node, 1, MPI_INT, rtrn, 10, MPI_COMM_WORLD, &status);
	    Log1 <<"Processing element "<<node<<" (Slave) reports."<<endl;
	    cout <<"Processing element "<<node<<" (Slave) reports."<<endl;
	  }
	}
	else {MPI_Ssend(&PE, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);}
	MPI_Barrier(MPI_COMM_WORLD);	//Ensure nodes do not get too far out of sync
#else
	Log1.open("LDCode.log");
	PE = 0;
	nPE = 1;
#endif
  Log  <<"INITIALIZED AND RUNNING '"<<argv[0]<<"'.\n"<<endl;
  cout <<"INITIALIZED AND RUNNING '"<<argv[0]<<"'.\n"<<endl;


	//READ INPUT DATA AND NON-DIMENSIONALIZE
	void Read_Input(int, char**, POTENTIAL***, PREPROCESS*,
      DISPERSION*, HARMONIC*, ANHARMONIC*, POSTPROCESS*);
	Read_Input(argc, argv, &Potential, Preprocess, Dispersion,
             Harmonic, Anharmonic, Postprocess);
#if defined PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);	//Ensure nodes do not get too far out of sync
#endif


  //INITIALIZE DATA AND COPY POTENTIALS
  Parameter = new PARAMETER(Potential);
  Lattice->Initialize();
  Material->Initialize();
  Unit_Cell->Initialize();
  Symmetry->Initialize();
  Anharmonic->Temperature /= Parameter->temp;
  void InitializePotential(POTENTIAL**);
  InitializePotential(Potential);


  //PRE-PROCESSING (call individual functions so you can see the structure)
  //Output the structure (unit cell/total, initial then optimized if applicable)
	//Optimize the structure of the unit cell (none, steepest decent, MD)
	if(Preprocess->Initialize()) {
    Preprocess->Optimize(Potential);
    Preprocess->Structure();
    Preprocess->Energy_Force(Potential);
	}
	delete Preprocess; Preprocess=NULL;



	//PERFORM HARMONIC LATTICE DYNAMICS CALCULATIONS (add thermodynamics calculations here and in post-processing)
	if( Dispersion->Initialize() || Harmonic->Initialize() ) {
    Log  <<"\n\nComputing Harmonic Force Constants...\n"<<endl;
    cout <<"\n\nComputing Harmonic Force Constants...\n"<<endl;
    PD_HARMONIC **HFC = NULL;
    HFC = Neighbor_Force<PD_HARMONIC>(HFC, Potential);
    if(!PE) {for(rtrn=0;rtrn<Unit_Cell->natom;rtrn++) {HFC[rtrn]->Print(0);}}
    Dispersion->Compute(HFC);
    if(PE==0) {Harmonic->Compute(Anharmonic->iterations, HFC);}
#if defined PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);	//Ensure nodes do not get too far out of sync
#endif
    for(rtrn=0;rtrn<Unit_Cell->natom;rtrn++) {delete HFC[rtrn]; HFC[rtrn]=NULL;}
    delete[] HFC; HFC=NULL;
	}


	//PERFORM ANHARMONIC LATTICE DYNAMICS CALCULATIONS
//	int Anharmonic_LDCode(int argc, char **argv);
//	rtrn = Anharmonic_LDCode(argc, argv);
  if(Anharmonic->iterations>0) {
    Log  <<"\n\nComputing Anharmonic Force Constants...\n"<<endl;
    cout <<"\n\nComputing Anharmonic Force Constants...\n"<<endl;
    PD_ANHARMONIC **AFC = NULL;
    AFC = Neighbor_Force(AFC, Potential);
    if(!PE) {for(rtrn=0;rtrn<Unit_Cell->natom;rtrn++) {AFC[rtrn]->Print(0);}}
    int Conductivity(int iter, double Temperature);
    //CALCULATE FREQUENCY SHIFT AND LINEWIDTH (INITIAL GUESS IS EMBEDDED IN Initialize)
    Anharmonic->Initialize(AFC);
    Log  <<endl;
    cout <<endl;
    for(int i=Anharmonic->continue_+1;i<=Anharmonic->iterations+Anharmonic->continue_;i++) {
      Log  <<"\nComputing frequency shift and linewidth (iter "<<i<<")...\n"<<endl;
      cout <<"\nComputing frequency shift and linewidth (iter "<<i<<")...\n"<<endl;
      Anharmonic->Compute(i);
      if(PE==0) {
        Log  <<"\nComputing thermal conductivity (iter "<<i<<")...\n"<<endl;
        cout <<"\nComputing thermal conductivity (iter "<<i<<")...\n"<<endl;
        Conductivity(i, Anharmonic->Temperature);
      }
    }
    for(rtrn=0;rtrn<Unit_Cell->natom;rtrn++) {delete AFC[rtrn]; AFC[rtrn]=NULL;}
    delete[] AFC; AFC=NULL;
  }


  //POST-PROCESSING (extract data, conductivity, thermodynamics properties, QC for now)
  if( (Postprocess->master_flag)&&(PE==0) ) {
    Postprocess->Initialize();
    if(Postprocess->TC_beg>0) {
      for(int i=Postprocess->TC_beg;i<=Postprocess->TC_end;i++) {
        Log  <<"\nComputing thermal conductivity (iter "<<i<<")...\n"<<endl;
        cout <<"\nComputing thermal conductivity (iter "<<i<<")...\n"<<endl;
        int Conductivity(int iter, double Temperature);
        Conductivity(i, Anharmonic->Temperature);
      }
    }
    Postprocess->Data_Dir();
  }
  //OUTPUT DATA IN A CONVENIENT, READABLE FORM
	//PERFORM THERMAL CONDUCTIVITY CALCULATION
	//int Conductivity(int iter, double Temperature);
	//rtrn = Conductivity(5, Anharmonic->Temperature);
  //PERFORM QUANTUM CORRECTION CALCULATIONS
	//int QC_Temperature(int argc, char **argv);
	//rtrn = QC_Temperature(argc, argv);


	//CLOSE LOG FILE AND DEALLOCATE MEMORY
	delete Parameter;   Parameter   = NULL;
	delete Anharmonic;  Anharmonic  = NULL;
	delete Dispersion;  Dispersion  = NULL;
	delete Symmetry;    Symmetry    = NULL;
	delete Unit_Cell;   Unit_Cell   = NULL;
	for(argc=0;argc<Material->n_mat;argc++) {
	  delete Potential[argc]; Potential[argc]=NULL;
	}
	delete[] Potential; Potential   = NULL;
	delete Lattice;     Lattice     = NULL;
	delete Material;    Material    = NULL;

	if(PE==0) {
    Begin_time = time(0) - Begin_time;
    rtrn = Begin_time/86400;
    Log  <<"\n# PE = "<<nPE<<"\nTotal Run Time = ";
    cout <<"\n# PE = "<<nPE<<"\nTotal Run Time = ";
    Log  <<rtrn<<" d  ";
    cout <<rtrn<<" d  ";
    Begin_time = Begin_time%86400;
    rtrn = Begin_time/3600;
    Log  <<rtrn<<" h  ";
    cout <<rtrn<<" h  ";
    Begin_time = Begin_time%3600;
    rtrn = Begin_time/60;
    Log  <<rtrn<<" m  ";
    cout <<rtrn<<" m  ";
    Begin_time = Begin_time%60;
    rtrn = Begin_time;
    Log  <<rtrn<<" s."<<endl;
    cout <<rtrn<<" s."<<endl;
    Log1.close();
  }


#if defined PARALLEL
	//END PARALLEL CODE
	MPI_Finalize();
#endif

	return(0);
}
