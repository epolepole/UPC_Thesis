/* These functionc are responsible for the different boundary conditions.*/

#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include "ShellFunctions.h" // convenience

void InletBC(CellProps *Cells, int i);

void WallBC(CellProps *Cells, int i, int* opp);

void EdgeBC(CellProps *Cells, int i);

void CornerBC(CellProps *Cells, int i);

void OutletBoundaries(CellProps *Cells, int j, int i);

void CurvedWallBoundaries(CellProps *Cells, int j, int i, int* opp);

#endif
