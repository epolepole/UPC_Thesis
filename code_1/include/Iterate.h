/* This function is the SOLVER itself */

#ifndef ITERATE_H
#define ITERATE_H

#include <stdio.h>
#include <time.h>
#include "BoundaryConditions.h" // convenience

#define init_measure_time tInstant1 = clock()
#define end_measure_time(t_var) tInstant2 = clock();\
                                t_var = t_var + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC

#if __SAVE_ITER__ == 0
#define SAVE_ITERATION
#else //__SAVE_ITER__
#define SAVE_ITERATION save_iteration(postproc_prog, Iteration, AutosaveEvery)
#endif //__SAVE_ITER__


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
FILE* testing_file;              // file for time measurement results
char IterOutputFile[50];                // write results to this file after the iterations
char AutosaveOutputFile[50];            // autosave filename
char IterationOutputFile[50];            // autosave filename
char OutputFile[50];                    // initial data will be written to this file
char FinalOutputFile[50];               // final data will be written to this file
char fnMemCopyRes[50];
char testingFileName[50];
char logFile[50];                       // path of the .log file
int  AutosaveI;                         // autosave i variable, will be incremented after every autosave
int* ppp;                               // pointer of the postproc_prog variable
int iter_counter;

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



//Init/Alloc/Free functions

void init_vars(int *postproc_prog);
void time_meas_vars_init();

void alloc_cells();
void allocate_vars();
void allocate_residuals();
void allocate_lattice_vars();
void allocate_shared();


void free_vars();
void free_mesh_data_matrices();


//IO functions
void read_data(const char *NodeDataFile, const char *BCconnectorDataFile);

void print_init_info_to_log(float Uavg, float Vavg, float Wavg, float rho_ini, float Viscosity, int InletProfile,
                            int CollisionModel, int CurvedBoundaries, int OutletProfile, int Iterations,
                            int AutosaveAfter, int AutosaveEvery, int postproc_prog, int CalculateDragLift,
                            float ConvergenceCritVeloc, float ConvergenceCritRho);




void auto_save(int AutosaveAfter, int AutosaveEvery, int postproc_prog);
void save_iteration(int postproc_prog,int Iterations, int AutosaveEvery);


void write_boundary_cells_to_results(int postproc_prog);
void save_init_data(int postproc_prog);
void export_data(int postproc_prog, int AutosaveEvery);
void print_times(const char* fname, int AutosaveEvery);

void print_cells_info(CellProps* Cells);
void print_boundary_type(CellProps* Cells);
void print_cell_line(FILE* file, const CellProps* Cell);


void putCellsToShared();
void getSharedToCells();
void putCellsToWCells();



//Step functions
void CollisionStep(int CollisionModel);
void StreamingStep();
void HandleBoundariesStep(int OutletProfile, int CurvedBoundaries);
void UpdateMacroscopicStep(int CalculateDragLift);
void CalculateDragLiftForcesStep(int CalculateDragLift);



//Other functions
CellProps* cell_from_id(CellProps* Cells, int ID);
void calc_collision_freq(float Viscosity);


void main_while_loop(int CollisionModel, int CurvedBoundaries, int OutletProfile, int *Iterations, int AutosaveAfter,
                     int AutosaveEvery, int postproc_prog, int CalculateDragLift, float ConvergenceCritVeloc,
                     float ConvergenceCritRho);

void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg,         float Vavg,  float Wavg,
               float rho_ini, float Viscosity,
               int InletProfile,   int CollisionModel, int CurvedBoundaries,
               int OutletProfile,  int Iterations,     int AutosaveAfter,
               int AutosaveEvery,  int postproc_prog,  int CalculateDragLift,
               float ConvergenceCritVeloc, float ConvergenceCritRho);

void print_Cell(char*,CellProps*,int);

#endif
