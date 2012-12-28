/*						             Neighbor_List.cpp           	  					*/
/*							              02/15/2009	              						*/
/*********************************************************************
*    File that contains the N_LIST functions.   				          	 *
*********************************************************************/

/*DEFINE HEADERS*/
#include "LDCode.h"


/*N_LIST CONSTRUCTOR: INITIALIZES VARIABLES*/
N_LIST::N_LIST(void) {
  l = new int[DIM];
  x = new double[DIM];
  next = NULL;
  return;
}


/*N_LIST DESTRUCTOR: DEALLOCATES MEMORY*/
N_LIST::~N_LIST(void) {
  if(next!=NULL) {delete next; next = NULL;}
  delete[] l;  l = NULL;
  delete[] x;  x = NULL;
  return;
}
