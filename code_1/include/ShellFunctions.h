/* The Cell structure is defined here. The functions are responsible for
scalar-vector-matrix allocation and this header includes the multifunctional
"min" and "max" functions too (comparison of two number). */
#ifndef SHELLFUNCTIONS_H
#define SHELLFUNCTIONS_H

#include "BlockSizeDefiner.h"
#include <upc_strict.h>

#include <stdbool.h>  // bool variables
#include <string.h>   // string handling


#if __UPC__ == 0


#define THREADS  8
#define MYTHREAD  (int)8*rand()

#define shared_block(var)
#define shared

#define upc_barrier
//#define upc_forall(a) ();

#define UPUT(a,b,c)
#define UGET(a,b,c)
#define UFREE(a)
#define UAFREE(a)

#else
#define shared_block(var) shared [var]
#define UPUT(a,b,c) upc_memput(a,b,c)
#define UGET(a,b,c) upc_memget(a,b,c)
#define UFREE(a) upc_free(a)
#define UAFREE(a) upc_all_free(a)
#endif

//#define main_thread(c) if (MYTHREAD == 0) { c }
#define main_thread if(MYTHREAD == 0)


////////////////////////////////////////////////////
///////////////////// DEFINE ///////////////////////
////////////////////////////////////////////////////

//#define UPC_MAX_BLOCK_SIZE 145000000


//#define BLOCKSIZE 802*((int)(83/THREADS+1))
//#define NN 802
//#define MM 83

////////////////////////////////////////////////////
////////////// SHARED DECLARATIONS /////////////////
////////////////////////////////////////////////////


typedef double MyReal;

/////////////////////////////////////////////////////////
/////////////////// STRUCT FOR CELLS ////////////////////
/////////////////////////////////////////////////////////

typedef struct
{
    //int Fluid;
    int    StreamLattice[19];
    int    ID;
    int    Boundary;
    int    BoundaryID;
    int    BC_ID[19];
    int    ThreadNumber;
    MyReal Q[19];
    MyReal CoordX;
    MyReal CoordY;
    MyReal CoordZ;
    MyReal U;
    MyReal V;
    MyReal W;
    MyReal Rho;
    MyReal Uo;
    MyReal Vo;
    MyReal Wo;
    MyReal DragXF;
    MyReal DragYF;
    MyReal LiftF;
    MyReal F[19];
    MyReal Feq[19];
    MyReal METAF[19];
    MyReal Fneighbours[19];
}CellProps ;

//shared [BLOCKSIZE+2*NN] CellProps  *SCells;
shared_block(5)             double sResiduals[5*THREADS]; // variable to store residuals
shared                      double shared_total_Residuals[5]; // variable to store residuals


//New things
shared_block(CELL_TOT_SIZE) CellProps  *WCells; // Writing Cells: cells to write data
shared_block(B_CELLS_SIZE)  CellProps  *BCells; // Boundary cells

shared  int    *NumNodes;       // This will store the number of lines of the read files
shared  int    *NumConn;        // This will store the number of lines of the read files
shared  int    *NumInletNodes;  // number of inlet nodes
shared  double *Delta;          // grid spacing
shared  double *MaxInletCoordY; // maximum inlet coordinate in y
shared  double *MinInletCoordY; // minimum inlet coordinate in y

shared int loop;





int    *Create1DArrayInt(int length);
float  *Create1DArrayFloat(int length);
double *Create1DArrayDouble(int length);

int    **Create2DArrayInt(int width, int height);
float  **Create2DArrayFloat(int width, int height);
double **Create2DArrayDouble(int width, int height);

int    ***Create3DArrayInt(int width, int height, int depth);
float  ***Create3DArrayFloat(int width, int height, int depth);
double ***Create3DArrayDouble(int width, int height, int depth);
bool   ***Create3DArrayBool(int width, int height, int depth);

void CreateDirectory(char* MainWorkDir);

void StringAddition(char* first, char* second, char* result);

int lID(int i, int j, int k);


#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

/*#define eq(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a == _b ? 1 : 0; })*/

#define eq(a,b) (a == b)

#endif
