/* These functionc are responsible for the different boundary conditions.*/

#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#define multVec(a,b) a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


#include "ShellFunctions.h" // convenience



int tang[3];


void InletBC(CellProps *Cells, int i);

void WallBC(CellProps *Cells, int i, int* opp);

void EdgeBC(CellProps *Cells, int i);

void CornerBC(CellProps *Cells, int i);

void OutletBoundaries(CellProps *Cells, int j, int i);

void CurvedWallBoundaries(CellProps *Cells, int j, int i, int* opp);


void generalWall(CellProps* Cells, int i);
void getWall(CellProps *Cells, int i, int wall_ID, double* v);
void getTan(const int* c, const int* n, int* t);

#endif
