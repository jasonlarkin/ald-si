/*                            Symmetry.h                            */
/*                            12/06/2008                            */
/*********************************************************************
*    Header file for the SYMMETRY class.  This class reads and       *
*  stores the proper symmetry conditions to use.                     *
*********************************************************************/

#ifndef SYMMETRY_H
#define SYMMETRY_H


/*DECLARE HEADERS*/
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;


/*DEFINE PREPROCESSOR COMMANDS*/
#define SYM3D 3
#define SYM2D 2


/*DEFINE CLASSES*/
class SYMMETRY {
public:
	SYMMETRY();                 //Constructor
	~SYMMETRY();                //Destructor
	string ReadInput(ifstream &Input);//
	int KeywordToInt(ifstream &Input);//Defines keywords
	void Output();        //Writes the variables
	void Initialize();    //Initializes variables after they have been defined (non-dimensionalizes)
	double Sym_Negative(int *k);//
	double Sym_3Phonon(int *k1, int *k2);//
	bool Sym_3D(int *k, int *pnt, int *k_IBZ);//Three equivalent directions
	bool Sym_2D(int *k, int *pnt, int *k_IBZ);//Two equivalent directions
	bool Sym_None(int *k, int *pnt, int *k_IBZ);//No symmetry
	int *k_test;                //Condition used in Sym_Negative only
	char Symb[4];               //Character string of equivalent directions
	bool wv_out;                //Flag, true to output unique wave vectors
	int sym;                    //Symmetry type identifier
	int *dir;                   //Equivalent directions (for 2D symmetry)
	bool (SYMMETRY::*Sym_ptr)(int *k, int *l, int *s);//Function pointer to proper frequency symmetry
	bool Sym_Freq(int *k, int *l=NULL, int *s=NULL);
	int TotalDOF;               //Total # of degrees of freedom = DIM*n*N
	int Sym_E;                  //Total # of unique eigenvectors
	int Sym_F;                  //Total # of unique frequencies
	int **unique;               //Stores the unique wave vectors
	int ***F_map;               //Points to the position of the data
	int ***E_map;               //Points to the correct eigenvector for each wave vector in 1st BZ
};


#endif // SYMMETRY_H
