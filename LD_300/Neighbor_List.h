/*                          Neighbor_List.h                         */
/*                            02/15/2009                            */
/*********************************************************************
*    Header file for the N_LIST class which is a linked list of      *
*  neighbors.                                                        *
*********************************************************************/

#ifndef NEIGHBOR_LIST_H_INCLUDED
#define NEIGHBOR_LIST_H_INCLUDED


/*DECLARE HEADERS*/


/*DECLARE GLOBAL VARIABLES*/


/*DEFINE PREPROCESSOR COMMANDS*/


/*N_LIST CLASS USED AS LINKED LIST IN NEIGHBOR LIST*/
class N_LIST {
public:
  int b;                      //Atom number in unit cell
	int *l;                     //Unit cell
	double *x;                  //Atomic position vector
	N_LIST *next;               //Pointer to next entry in list
  N_LIST(void);               //Constructor
  ~N_LIST(void);              //Destructor
};

#endif // NEIGHBOR_LIST_H_INCLUDED
