/*                           Structure.cpp                          */
/*                            12/02/2008                            */
/*********************************************************************
*    Subroutine that builds the structure defined by Lattice and     *
*  Unit_Cell and outputs it to "Structure.xyz".                      *
*********************************************************************/

/*DECLARE HEADERS*/
#include "LDCode.h"


/*DEFINE PREPROCESSOR VARIABLES*/


/*DEFINE SUBROUTINE Structure*/
void Structure(void) {
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
