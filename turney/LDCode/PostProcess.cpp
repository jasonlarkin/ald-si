/*                          PostProcess.cpp                         */
/*                            06/20/2009                            */
/*********************************************************************
*    Subroutines for the POSTPROCESS class.                          *
*********************************************************************/

/*DECLARE HEADERS*/
#include "LDCode.h"
#include "PostProcess.h"


/*DEFINE PREPROCESSOR COMMANDS*/


/*POSTPROCESS CLASS Constructor*/
POSTPROCESS::POSTPROCESS() {
  master_flag = false;        //Master flag for turning POSTPROCESS on/off
  TC_beg = 0;                 //Beginning iteration for thermal conductivity calculation
  TC_end = 0;                 //Ending iteration for thermal conductivity calculation
  data_itr = -1;              //Iteration number for extracting data
  data_beg = NULL;            //Beginning wave vector for extracting data along arbitrary direction
  data_inc = NULL;            //Ending wave vector for extracting data along arbitrary direction
  data_num = -1;              //Number of times to use increment
  data_fBZ = DATA_NONE;       //Flag for outputing the data for the full BZ
  return;
}


/*POSTPROCESS CLASS Destructor*/
POSTPROCESS::~POSTPROCESS() {
  if(data_beg!=NULL) {
    delete[] data_beg; data_beg=NULL;
    delete[] data_inc; data_inc=NULL;
  }
  return;
}


/*SUBROUTINE ReadInput: USED TO IDENTIFY DATA TO STORE*/
string POSTPROCESS::ReadInput(ifstream &Input) {
  //DECLARE LOCAL VARIABLES
  int i;
  string str;
  string Read_Next(ifstream &Input);
  data_beg = new int[DIM];
  data_inc = new int[DIM];
  for(i=0;i<DIM;i++) {data_beg[i] = data_inc[i] = 0;}


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
    if     (str.compare(0, 4, "tc_beg"  , 4)==0) {TC_beg = int(atof(Read_Next(Input).c_str()));}
    else if(str.compare(0, 4, "tc_end"  , 4)==0) {TC_end = int(atof(Read_Next(Input).c_str()));}
    else if(str.compare(0, 8, "data_itr", 8)==0) {data_itr = int(atof(Read_Next(Input).c_str()));}
    else if(str.compare(0, 8, "data_beg", 8)==0) {for(i=0;i<DIM;i++) {data_beg[i] = int(atof(Read_Next(Input).c_str()));}}
    else if(str.compare(0, 8, "data_inc", 8)==0) {for(i=0;i<DIM;i++) {data_inc[i] = int(atof(Read_Next(Input).c_str()));}}
    else if(str.compare(0, 8, "data_num", 8)==0) {data_num = int(atof(Read_Next(Input).c_str()));}
    else if(str.compare(0, 8, "data_fbz", 8)==0) {
      str = Read_Next(Input);
      i = 0; while(str[i]) {str[i]=tolower(str[i]);i++;}
      if(str.compare(0, 4, "none"       , 4)==0) {data_fBZ = DATA_NONE;} else
      if(str.compare(0, 4, "no_symmetry", 4)==0) {data_fBZ = DATA_NO_SYM;} else
      if(str.compare(0, 3, "symmetry"   , 3)==0) {data_fBZ = DATA_SYM;} else
      {Log <<"Expected 'none', 'no_symmetry', or 'symmetry' as argument to 'data_fBZ'."<<endl;exit(0);}
    }
    else{break;}
  }
  return(str);
}


/*LATTICE FUNCTION Output: WRITES DATA TO LOG FILE*/
void POSTPROCESS::Output() {
  if(master_flag) {
    int i;
    Log <<"\nPOSTPROCESS = true"<<endl;
    Log <<"  TC_beg   = "<<TC_beg<<endl;
    Log <<"  TC_end   = "<<TC_end<<endl;
    Log <<"  data_itr = "<<data_itr<<endl;
    Log <<"  data_beg =";
    for(i=0;i<DIM;i++) {Log <<" "<<data_beg[i];}
    Log <<endl;
    Log <<"  data_inc =";
    for(i=0;i<DIM;i++) {Log <<" "<<data_inc[i];}
    Log <<endl;
    Log <<"  data_num = "<<data_num<<endl;
    Log <<"  data_fBZ = ";
    if(data_fBZ==DATA_NONE) {Log <<"none"<<endl;}
    else if (data_fBZ==DATA_SYM) {Log <<"symmetry"<<endl;}
    else if (data_fBZ==DATA_NO_SYM) {Log <<"no_symmetry"<<endl;}
  }
  else {Log <<"\nPOSTPROCESS = false"<<endl;}
  return;
}


/*LATTICE FUNCTION Initialize: INITALIZES PARAMETERS AND NON-DIMENSIONALIZES*/
void POSTPROCESS::Initialize() {return;}


/*LATTICE FUNCTION Data_Dir: OUTPUTS DATA ALONG ARBITRARY DIRECTION*/
void POSTPROCESS::Data_Dir() {
  if( (data_num<0)&&(data_fBZ==DATA_NONE) ) {return;}
  Log  <<"\n\nPrinting Data Along Specified Direction and/or For Full BZ\n";
  cout <<"\n\nPrinting Data Along Specified Direction and/or For Full BZ\n";
  //DECLARE LOCAL VARIABLES
  int i, j, d;                //Counters
  int Sym_F = Symmetry->Sym_F;//Number of unique phonons
  bool flag_AHLD = true;      //Flag, true if anharmonic files found
  char name[128];             //Name of files to open
  int pos;                    //Position in list
  int *wv = new int[DIM];     //Current wave vector
  int *k  = new int[DIM];     //Holds wave vector in first BZ
  int *n  = new int[2*DIM];   //For use in mapping velocities
  double *QH_f;	              //Quasi-harmonic phonon frequencies
  double *AH_fsQ;	            //Quantum anharmonic frequencies
	double *AH_fsC;             //Classical anharmonic frequencies
	double *AH_lwQ;             //Quantum anharmonic linewidth
	double *AH_lwC;             //Classical anharmonic linewidth
	double **Vel;               //Velocities
	double x;                   //Generic double
  ofstream file_out;          //Output file
  ifstream QH_Freq;           //Frequency.txt
  ifstream Q_AHLD;            //Shift_Width.txt
  ifstream C_AHLD;            //Classical_Shift_Width.txt
  ifstream QH_Vel;            //Velocity.txt


  //OPEN FILES
  QH_Freq.open("Frequency.txt");
  if(!QH_Freq.is_open()) {Log <<"Could not open '"<<name<<"'."<<endl;exit(0);}
  QH_Freq.ignore(10000, '\n');
  QH_Vel.open("Velocity.txt");
  if(!QH_Vel.is_open())  {Log <<"Could not open '"<<name<<"'."<<endl;exit(0);}
  QH_Vel.ignore(10000, '\n');
  sprintf(name, "%s%i.txt", "Shift_Width", data_itr);
  Q_AHLD.open(name);
  sprintf(name, "%s%i.txt", "Classical_Shift_Width", data_itr);
  C_AHLD.open(name);
  if( (!Q_AHLD.is_open())||(!C_AHLD.is_open()) ) {
    flag_AHLD = false;
    Log  <<"Could not open quantum and/or classical AHLD file; suppressing all AHLD output."<<endl;
    cout <<"Could not open quantum and/or classical AHLD file; suppressing all AHLD output."<<endl;
  }
  else {
    Q_AHLD.ignore(10000, '\n');
    C_AHLD.ignore(10000, '\n');
  }


  //ALLOCATE MEMORY AND READ DATA INTO ARRAYS
  QH_f   = new double[Sym_F];	//Quasi-harmonic phonon frequencies
  if(flag_AHLD) {
    AH_fsQ = new double[Sym_F];	//Quantum anharmonic frequency shifts
    AH_fsC = new double[Sym_F]; //Classical anharmonic frequency shifts
    AH_lwQ = new double[Sym_F]; //Quantum anharmonic linewidth
    AH_lwC = new double[Sym_F]; //Classical anharmonic linewidth
  }
  else {AH_fsQ = AH_fsC = AH_lwQ = AH_lwC = NULL;}
	Vel    = new double*[Sym_F];//Velocities
	Vel[0] = new double[Sym_F*DIM];
	for(i=0;i<Sym_F;i++) {
	  if(i<Sym_F-1) {Vel[i+1] = Vel[i] + DIM;}
	  QH_Freq >>QH_f[i];
	  for(j=0;j<DIM;j++) {QH_Vel >>Vel[i][j];}
	  if(flag_AHLD) {
	    Q_AHLD >>x>>AH_fsQ[i]>>x>>AH_lwQ[i];
	    AH_fsQ[i] += x;
	    C_AHLD >>x>>AH_fsC[i]>>x>>AH_lwC[i];
	    AH_fsC[i] += x;
	  }
	}
	QH_Freq.close();
	Q_AHLD.close();
	C_AHLD.close();
	QH_Vel.close();


  //OUTPUT DATA ALONG ARBITRARY DIRECTION
  if(data_num>=0) {
    file_out.open("Data_Direction.xls");
    for(d=0;d<DIM;d++) {file_out <<"l"<<d<<'\t';}
    file_out <<"QH freq (rad/ps)\t";
    if(flag_AHLD) {
      file_out <<"Q freq shift (rad/ps)\tQ lifetime (ps)\t";
      file_out <<"C freq shift (rad/ps)\tC lifetime (ps)\t";
    }
    for(d=0;d<DIM;d++) {file_out <<"v"<<d<<" (m/s)\t";}
    file_out <<endl;
    for(d=0;d<DIM;d++) {k[d] = wv[d] = data_beg[d];}
    for(i=0;i<=data_num;i++) {
      for(d=0;d<DIM;d++) {
        if(k[d]<Lattice->lower[d]) {k[d] = Lattice->N[d]-((-k[d])%Lattice->N[d]);}
        else if(k[d]>Lattice->upper[d]) {k[d] = k[d]%Lattice->N[d]-Lattice->N[d];}
      }
      Symmetry->Sym_Freq(k, &pos, n);
      for(j=0;j<UC_DOF;j++) {
        for(d=0;d<DIM;d++) {file_out <<wv[d]<<'\t';}
        file_out <<QH_f[pos]<<'\t';
        if(flag_AHLD) {
          file_out <<AH_fsQ[pos]<<'\t'<<AH_lwQ[pos]<<'\t';
          file_out <<AH_fsC[pos]<<'\t'<<AH_lwC[pos]<<'\t';
        }
        for(d=0;d<DIM;d++) {
          file_out <<n[DIM+d]*Vel[pos][n[d]]<<'\t';
        }
        file_out <<endl;
        pos += 1;
      }
      for(d=0;d<DIM;d++) {
        k[d]  += data_inc[d];
        wv[d] += data_inc[d];
      }
    }
    file_out.close();
  }


  //OUTPUT DATA FOR FULL BRILLOUIN ZONE
  if(data_fBZ!=DATA_NONE) {
    file_out.open("Data_fullBZ.xls");
    for(d=0;d<DIM;d++) {file_out <<"l"<<d<<'\t';}
    file_out <<"QH freq (rad/ps)\t";
    if(flag_AHLD) {
      file_out <<"Q freq shift (rad/ps)\tQ lifetime (ps)\t";
      file_out <<"C freq shift (rad/ps)\tC lifetime (ps)\t";
    }
    for(d=0;d<DIM;d++) {file_out <<"v"<<d<<" (m/s)\t";}
    file_out <<endl;
    for(wv[0]=Lattice->upper[0];wv[0]>=Lattice->lower[0];wv[0]--) {
    for(wv[1]=Lattice->upper[1];wv[1]>=Lattice->lower[1];wv[1]--) {
    for(wv[2]=Lattice->upper[2];wv[2]>=Lattice->lower[2];wv[2]--) {
      if(Symmetry->Sym_Negative(wv)<0.0) {if(data_fBZ==DATA_SYM) {continue;}}
      if(Symmetry->Sym_Freq(wv, &pos, n)) {if(data_fBZ==DATA_SYM) {continue;}}
      for(j=0;j<UC_DOF;j++) {
        for(d=0;d<DIM;d++) {file_out <<wv[d]<<'\t';}
        file_out <<QH_f[pos]<<'\t';
        if(flag_AHLD) {
          file_out <<AH_fsQ[pos]<<'\t'<<AH_lwQ[pos]<<'\t';
          file_out <<AH_fsC[pos]<<'\t'<<AH_lwC[pos]<<'\t';
        }
        for(d=0;d<DIM;d++) {
          file_out <<n[DIM+d]*Vel[pos][n[d]]<<'\t';
        }
        file_out <<endl;
        pos += 1;
      }
    }
    }
    }
  }


  //DEALLOCATE MEMORY
  delete[] k;      k     =NULL;
  delete[] n;      n     =NULL;
  delete[] wv;     wv    =NULL;
  delete[] QH_f;   QH_f  =NULL;
  delete[] AH_fsQ; AH_fsQ=NULL;
  delete[] AH_fsC; AH_fsC=NULL;
  delete[] AH_lwQ; AH_lwQ=NULL;
  delete[] AH_lwC; AH_lwC=NULL;
  delete[] Vel[0]; Vel[0]=NULL;
  delete[] Vel;    Vel   =NULL;
  return;
}
