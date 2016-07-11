/* These functionc are responsible for the initialization and includes the collision model.
The microscopic and macroscopic variables are refreshed based on the last two functions.*/

#ifndef CELLFUNCTIONS_H
#define CELLFUNCTIONS_H

#include "ShellFunctions.h"
#include "BlockSizeDefiner.h"

// D3Q19 Variables of the lattice
double* w;              // weight values for the directions
int*   cx;              // x coordinate of the discrete lattice directions
int*   cy;              // y coordinate of the discrete lattice directions
int*   cz;              // z coordinate of the discrete lattice directions
int*  opp;              // opposite vector
int*    c;              // shift of lattice directions written in vector form
int** norm;
int** j_wall_unknown;

//void D3Q19Vars(double* w, int* cx, int* cy, int* cz, int* opp, int* c);
void D3Q19Vars();

void MRTInitializer(double** tm, double** stmiv, double Omega);


void CellIni(CellProps *Cells,
             float  **Nod,
             float  **Con,
             float  Uavg,
             float  Vavg,
             float  Wavg,
             int    InletProfile,
             int    CollisionModel,
             int*   opp,
             float  rho_ini);

void BGKW(int i, double Omega);

void TRT(CellProps *Cells, int i, double* w, int* cx, int* cy, int* opp, double Omega, double OmegaA);

void MRT(CellProps *Cells, int i, double** tm, double** stmiv);

void UpdateF(CellProps *Cells, int i);

void UpdateMacroscopic(CellProps *Cells, int i, int CalculateDragLift);

void CalculateDragLiftForces(CellProps *Cells, int j, int i, int CalculateDragLift, shared int* n, shared int* m);

int getIndex(const int x, const int y, const int z);
int getThread(int index);

#endif
