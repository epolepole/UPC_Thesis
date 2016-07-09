/* The Cell structure is defined here. The functions are responsible for
scalar-vector-matrix allocation and this header includes the multifunctional
"min" and "max" functions too (comparison of two number). */
#ifndef SHELLFUNCTIONS_H
#define SHELLFUNCTIONS_H

#include <stdbool.h>  // bool variables
#include <string.h>   // string handling
#include <upc.h>



#if __UPC__ == 0


#define MYTHREAD  (int) rand()*8
#define THREADS  8

#define shared_block(var)
#define shared

#define upc_barrier
#define upc_forall(a) ();

#else
#define shared_block(var) shared [var]
#endif

//#define main_thread(c) if (MYTHREAD == 0) { c }
#define main_thread if(MYTHREAD == 0)


////////////////////////////////////////////////////
///////////////////// DEFINE ///////////////////////
////////////////////////////////////////////////////

//#define UPC_MAX_BLOCK_SIZE 145000000

#include "BlockSizeDefiner.h"

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
    int    Fluid;
    int    Corner;
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
shared_block(BLOCKSIZE)     CellProps  *WCells; // Writing Cells: cells to write data
shared_block(BLOCKSIZE + 2*LAYER)     CellProps  *WCells2; // Writing Cells: cells to write data
shared_block(2*LAYER)           CellProps  *BCells; // Boundary cells
shared_block(5)              double sResiduals[5*THREADS]; // variable to store residuals

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


#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })


#endif
