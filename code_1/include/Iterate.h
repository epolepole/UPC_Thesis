/* This function is the SOLVER itself */

#ifndef ITERATE_H
#define ITERATE_H

#include "ShellFunctions.h" // convenience

#define init_measure_time tInstant1 = clock()
#define end_measure_time(t_var) tInstant2 = clock();\
                                t_var = (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC

#define func_between_time(func,t_container) \
    tInstant1 = clock(); // Start measuring time\
    func;\
    tInstant2 = clock();\
    t_container = t_container + (float)(tv.tInstant2-tv.tInstant1); / CLOCKS_PER_SEC;

////////////////////////////////////////////////////
/////////// LOCAL VARIABLES FOR EACH CPU ///////////
////////////////////////////////////////////////////

int i, j, k, iter;                      // variables for loops
FILE* resid_file;                       // file for residuals
FILE* log_file;                         // file for log
FILE* TimeMeasurementFile;              // file for time measurement results
char IterOutputFile[50];                // write results to this file after the iterations
char AutosaveOutputFile[50];            // autosave filename
char OutputFile[50];                    // initial data will be written to this file
char FinalOutputFile[50];               // final data will be written to this file
char logFile[50];                       // path of the .log file
int  AutosaveI;                         // autosave i variable, will be incremented after every autosave
int* ppp;                               // pointer of the postproc_prog variable

CellProps *Cells;                       // Struct for Cells

// Time measurement variables
float tInitialization;          // Time measurement of Initialization
float tCellsInitialization;     // Time measurement of Initialization
float tIteration;               // Time measurement of Iteration
float tCollision;               // Time measurement of Collision
float tUpdateF;                 // Time measurement of UpdateF
float tStreaming;               // Time measurement of Streaming
float tBoundaries;              // Time measurement of Boundaries
float tUpdateMacro;             // Time measurement of Update Macroscopic vars
float tResiduals;               // Time measurement of calculating residuals
float tWriting;                 // Time measurement of writing data
float tBCells;                  // Time measurement of handling boundaries
clock_t tInstant1, tInstant2;   // Time measurement points: universal
clock_t tIterStart, tIterEnd;   // Time measurement points: main loop

clock_t tStart;     // Time measurement: declaration, begin

// Variables for residuals
double *Residuals;
double *sumVel0;
double *sumVel1;
double *sumRho0;
double *sumRho1;

float **Nodes;          // matrices for the nodes
float **BCconn;         // matrices for connections

double Omega;
double OmegaA;          // collision frequency from the viscosity
double **tm;            // variable for the MRT collision model
double **stmiv;         // variable for the MRT collision model

// D3Q19 Variables of the lattice
double* w;              // weight values for the directions
int*   cx;              // x coordinate of the discrete lattice directions
int*   cy;              // y coordinate of the discrete lattice directions
int*   cz;              // z coordinate of the discrete lattice directions
int*  opp;              // opposite vector
int*    c;              // shift of lattice directions written in vector form



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
