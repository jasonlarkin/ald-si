/*                          PreProcess.cpp                          */
/*                            06/30/2009                            */
/*********************************************************************
*    Subroutines for the PREPROCESS class.                           *
*********************************************************************/

/*DECLARE HEADERS*/
#include "LDCode.h"
#include "PreProcess.h"
#include "Neighbor_Force.h"
#include <iomanip>
#include <cmath>


/*DEFINE PREPROCESSOR COMMANDS*/


/*POSTPROCESS CLASS Constructor*/
PREPROCESS::PREPROCESS() {
  master_flag = false;        //Master flag for turning PREPROCESS on/off
  opt_itr = 0;                //Max number of iterations for optimization
  struct_flag = false;        //Flag for outputing structure
  ef_flag = false;            //Flag for computing and outputing energies and forces
  return;
}


/*POSTPROCESS CLASS Destructor*/
PREPROCESS::~PREPROCESS() {return;}


/*SUBROUTINE ReadInput: USED TO IDENTIFY DATA TO STORE*/
string PREPROCESS::ReadInput(ifstream &Input) {
  //DECLARE LOCAL VARIABLES
  int i;
  string str;
  string Read_Next(ifstream &Input);


  //PARSE THROUGH FILE
  str = Read_Next(Input);
  i = 0; while(str[i]) {str[i]=tolower(str[i]);i++;}
  if     (str.compare(0, 1, "true",  1)==0) {master_flag = true;}
  else if(str.compare(0, 1, "false", 1)==0) {master_flag = false;}
  else {Log <<"'master_flag' accepts true or false"<<endl;exit(0);}
  while(!Input.eof()) {
    //identify keyword in POSTPROCESS category
    str = Read_Next(Input);
    i = 0; while(str[i]) {str[i]=tolower(str[i]);i++;}
    if     (str.compare(0, 3, "optimization", 3)==0) {opt_itr = int(atof(Read_Next(Input).c_str())+0.001);}
    else if(str.compare(0, 3, "structure"   , 3)==0) {
      str = Read_Next(Input);
      i = 0; while(str[i]) {str[i]=tolower(str[i]);i++;}
      if     (str.compare(0, 1, "true",  1)==0) {struct_flag = true;}
      else if(str.compare(0, 1, "false", 1)==0) {struct_flag = false;}
      else {Log <<"'structure' accepts true or false"<<endl;exit(0);}
    }
    else if( (str.compare(0, 3, "energy", 3)==0)||(str.compare(0, 3, "force", 3)==0) ) {
      str = Read_Next(Input);
      i = 0; while(str[i]) {str[i]=tolower(str[i]);i++;}
      if     (str.compare(0, 1, "true",  1)==0) {ef_flag = true;}
      else if(str.compare(0, 1, "false", 1)==0) {ef_flag = false;}
      else {Log <<"'energy'/'force' accepts true or false"<<endl;exit(0);}
    }
    else{break;}
  }
  return(str);
}


/*PREPROCESS FUNCTION Output: WRITES DATA TO LOG FILE*/
void PREPROCESS::Output() {
  if(master_flag) {
    Log <<"\nPREPROCESS = true"<<endl;
    Log <<"  optimization = "<<opt_itr<<endl;
    Log <<"  structure    = "<<(struct_flag ? "true" : "false")<<endl;
  }
  else {Log <<"\nPREPROCESS = false"<<endl;}
  return;
}


/*PREPROCESS FUNCTION Initialize: INITALIZES PARAMETERS AND NON-DIMENSIONALIZES*/
bool PREPROCESS::Initialize() {return(master_flag);}


/*PREPROCESS SUBROUTINE Optimize: OPTIMIZES THE POSITIONS OF THE ATOMS BY SIMPLE GRADIENT METHOD*/
void PREPROCESS::Optimize(POTENTIAL **Pot) {
  if(opt_itr<=0) {return;}
  Log  <<"\nOptimizing Unit Cell Structure...\n"<<endl;
	cout <<"\nOptimizing Unit Cell Structure...\n"<<endl;
  //DECLARE LOCAL VARIABLES
  double TEST = 1.0e-10;      //Convergence criterion
  double STEP = 20.000;       //Initial step size
  int natom=Unit_Cell->natom; //Number of atoms in the unit cell
  int iter, n, d;             //Counters
  int n_new, d_new;           //New atom and direction to move
  int n_old, d_old;           //Old atom and direction to move
  double X_step;              //Atomic positions at the last and current iteration
  double norm_old, norm_new;  //Vector norms of the energy gradient
  double *mass;               //Temporarily stores the mass of the atoms
  double gd;                  //Generic double
  double *max_step;           //Max step to take in x, y, and z


  //ALLOCATE MEMORY AND INITIALIZE
  //Temporarily set all masses to 1.0 (reset when deallocating)
  mass = new double[Material->n_mat];
  for(n=0;n<Material->n_mat;n++) {
    mass[n] = Material->mass[n];
    Material->mass[n] = 1.0;
  }
  //Set maximum step size to 1/10th the unit cell dimesnions
  max_step = new double[DIM];
  for(d=0;d<DIM;d++) {
    max_step[d] = gd = 0.0;
    for(n=0;n<DIM;n++) {
      if(Lattice->a[n][d]>max_step[d]) {max_step[d] = Lattice->a[n][d];}
      else if(Lattice->a[n][d]<gd) {gd = Lattice->a[n][d];}
    }
    max_step[d] -= gd;
    max_step[d] *= 0.1;
  }
  //Open temporary output file
  Log.close();
  Log.open("Opt_temp.txt");
  //Compute initial force constants
  Log  <<"Computing force constants for iteration 0."<<endl;
	cout <<"Computing force constants for iteration 0."<<endl;
	PD_ENERGY_FORCE **FC = NULL;
	FC = Neighbor_Force(FC, Pot);
	//Allocate memory and initialize vector norm
	norm_old = 0.0;
	n_old = 0;
	d_old = 0;
  for(n=0;n<natom;n++) {
    for(d=0;d<DIM;d++) {
      gd = FC[n]->F1_0->fc[d];
      if(fabs(gd)>fabs(FC[n_old]->F1_0->fc[d_old])) {n_old=n;  d_old=d;}
      norm_old += gd*gd;
    }
  }
  norm_new = norm_old = sqrt(norm_old/natom);
  n_new = n_old;
  d_new = d_old;
  cout <<"Vector norm of forces after 0 iterations = "<<norm_old<<endl;


  //LOOP
  iter = 0;
  while(norm_old>TEST) {
    iter += 1;
    if(iter>opt_itr) {
      Log  <<"Maximum number of iterations ("<<opt_itr<<") reached."<<endl;
      cout <<"Maximum number of iterations ("<<opt_itr<<") reached."<<endl;
      break;
    }


/* Simple gradient method */
    //COMPUTE STEP
    X_step = FC[n_old]->F1_0->fc[d_old];//*/


    //UPDATE ATOMIC POSTIONS
    cout <<norm_old<<" "<<norm_new<<endl;
    cout <<STEP*X_step*Parameter->length*-1.0e10<<" "<<n_old<<" "<<d_old<<endl;
    Unit_Cell->X[n_old][d_old] -= STEP*X_step;


    //COMPUTE NEW FORCE CONSTANTS AND VECTOR NORM
    Log  <<"\nComputing force constants for iteration "<<iter<<"."<<endl;
    cout <<"\nComputing force constants for iteration "<<iter<<"."<<endl;
    for(n=0;n<natom;n++) {delete FC[n];}
    delete[] FC;
    FC = Neighbor_Force(FC, Pot);
    n_new = 0;
    d_new = 0;
    norm_new = 0.0;
    for(n=0;n<natom;n++) {
      for(d=0;d<DIM;d++) {
        gd = FC[n]->F1_0->fc[d];
        if(fabs(gd)>fabs(FC[n_new]->F1_0->fc[d_new])) {n_new=n;  d_new=d;}
        norm_new += gd*gd;
      }
    }
    norm_new = sqrt(norm_new/natom);
    if( (norm_new>norm_old)||(STEP*X_step>max_step[d_old]) ) {   //Reduce step size
      Unit_Cell->X[n_old][d_old] += STEP*X_step;
      for(n=0;n<natom;n++) {delete FC[n];}
      delete[] FC;
      FC = Neighbor_Force(FC, Pot);
      STEP /= 2.0;
      if(STEP<1.0e-10) {
        Log  <<"Minimum step size ("<<STEP<<") reached."<<endl;
        cout <<"Minimum step size ("<<STEP<<") reached."<<endl;
        break;
      }
      iter -= 1;
    }
    else {
      norm_old = norm_new;
      n_old = n_new;
      d_old = d_new;
      Log  <<"Vector norm (STEP) of forces after "<<iter
           <<" iterations = "<<norm_old<<" ("<<STEP<<")."<<endl;
      cout <<"Vector norm (STEP) of forces after "<<iter
           <<" iterations = "<<norm_old<<" ("<<STEP<<")."<<endl;
    }

  }


  //COPY NEW COORDINATES TO Unit_Cell DATA
  Log.close();
	remove("Opt_temp.txt");
  Log.open("LDCode.log", ios::app);
  Log  <<"Updating atomic positions."<<endl;
  cout <<"\nUpdating atomic positions."<<endl;
  ofstream Opt("Opt_Struct.xyz");
  Opt <<natom<<"\nVector norm (step) is "<<norm_old<<" ("<<STEP<<") after "<<iter
      <<" iterations.  Optimized atomic positions in angstroms:";
  for(n=0;n<natom;n++) {
    Opt <<"\n  "<<Material->symb[Unit_Cell->mat[n]];
    for(d=0;d<DIM;d++) {Opt <<" "<<Unit_Cell->X[n][d]*Parameter->length*1.0e10;}
  }
  Opt <<"\n\n\nOptimized atomic positions in meters:";
  Opt.precision(14);
	Opt.setf(ios::scientific);
  for(n=0;n<natom;n++) {
    Opt <<"\n  "<<Material->symb[Unit_Cell->mat[n]];
    for(d=0;d<DIM;d++) {Opt <<" "<<setw(21)<<Unit_Cell->X[n][d]*Parameter->length;}
  }
  Opt.close();


	//DEALLOCATE MEMORY AND RETURN
  for(n=0;n<Material->n_mat;n++) {Material->mass[n] = mass[n];}
  delete[] mass;      mass     =NULL;
  delete[] max_step;  max_step =NULL;
  for(n=0;n<natom;n++) {delete FC[n]; FC[n]=NULL;}
	delete[] FC;        FC       =NULL;
	return;
}


/*PREPROCESS FUNCTION Structure: OUTPUTS THE STRUCTURE OF THE CRYSTAL*/
void PREPROCESS::Structure(void) {
  if(PE) {return;}
  if(!struct_flag) {return;}
  Log  <<"\nOutputing crystal structure to file 'Structure.xyz'...\n"<<endl;
  cout <<"\nOutputing crystal structure to file 'Structure.xyz'...\n"<<endl;
  //DECLARE LOCAL VARIABLES
  int b, d0, d1, l[3];
  double x;
  ofstream Struct_Out("Structure.xyz");


  //OUTPUT STRUCTURE
  b = Unit_Cell->natom;
  for(d0=0;d0<DIM;d0++) {b *= Lattice->N[d0];}
  Struct_Out <<b<<"\nAtomic positions in angstroms"<<endl;
  for(l[0]=0;l[0]<Lattice->N[0];l[0]++) {
    for(l[1]=0;l[1]<Lattice->N[1];l[1]++) {
      for(l[2]=0;l[2]<Lattice->N[2];l[2]++) {
        for(b=0;b<Unit_Cell->natom;b++) {
          Struct_Out <<Material->symb[Unit_Cell->mat[b]];
          for(d0=0;d0<DIM;d0++) {
            x = Unit_Cell->X[b][d0];
            for(d1=0;d1<DIM;d1++) {x += l[d1]*Lattice->a[d1][d0];}
            Struct_Out <<'\t'<<x*Parameter->length*1.0e10;
          }
          Struct_Out <<endl;
        }
      }
    }
  }
  Struct_Out.close();
  return;
}


/*PREPROCESS FUNCTION Energy_Force: COMPUTES THE ENERGIES AND FORCES ON ATOMS*/
void PREPROCESS::Energy_Force(POTENTIAL **Pot) {
  if(!ef_flag) {return;}
  int i;
  Log  <<"\n\nComputing Energies and Forces...\n"<<endl;
  cout <<"\n\nComputing Energies and Forces...\n"<<endl;
  PD_ENERGY_FORCE **EF = NULL;
  EF = Neighbor_Force(EF, Pot);
  for(i=0;i<Unit_Cell->natom;i++) {
    if(!PE) {EF[i]->Print(0);}
    delete EF[i]; EF[i]=NULL;
  }
  delete[] EF; EF=NULL;
  return;
}
