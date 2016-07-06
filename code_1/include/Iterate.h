/* This function is the SOLVER itself */

#ifndef ITERATE_H
#define ITERATE_H

#include "ShellFunctions.h" // convenience

void putCellsToShared(CellProps *Cells);

void getSharedToCells(CellProps *Cells);

void putCellsToWCells(CellProps *Cells);

void CollisionStep(CellProps *Cells, double* w, int* cx, int* cy, int* opp,
                   double Omega, double OmegaA, double **tm,
                   double **stmiv, int CollisionModel);

void StreamingStep(CellProps *Cells, int* c);

void HandleBoundariesStep(CellProps *Cells, int* cx, int* cy, int* c, int* opp, int OutletProfile, int CurvedBoundaries);

void UpdateMacroscopicStep(CellProps *Cells, int* cx, int* cy, int CalculateDragLift);

void CalculateDragLiftForcesStep(CellProps *Cells, shared int* m, shared int* n, int CalculateDragLift);

void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg,         float Vavg,  float Wavg,
               float rho_ini, float Viscosity,
               int InletProfile,   int CollisionModel, int CurvedBoundaries,
               int OutletProfile,  int Iterations,     int AutosaveAfter,
               int AutosaveEvery,  int postproc_prog,  int CalculateDragLift,
               float ConvergenceCritVeloc, float ConvergenceCritRho);

void print_cells_info(CellProps* Cells);
void print_boundary_type(CellProps* Cells);
void print_cell_line(FILE* file, const CellProps* Cell);
CellProps* cell_from_id(CellProps* Cells, int ID);

#endif
