#include <stdio.h>                      // printf();
#include <math.h>                       // need to compile with -lm
#include <stdlib.h>                     // for calloc();
#include <stdbool.h>                    // Include for bool type variables!
#include <string.h>                     // String operations
#include <time.h>                       // time functions
#include <upc_relaxed.h>                // Required for UPC 
#include <CellFunctions.h>
#include <ComputeResiduals.h>

////////////////////////////////////////////////////
////////////////// OWN HEADERS /////////////////////
////////////////////////////////////////////////////

#include "ShellFunctions.h"     // For convenience
#include "Iterate.h"            // Iteration takes place
#include "FilesReading.h"       // For reading files
#include "CellFunctions.h"      // For cell modifications


////////////////////////////////////////////////////
/////////////////// ITERATION //////////////////////
////////////////////////////////////////////////////

void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg,         float Vavg,         float Wavg,
               float rho_ini,      float Viscosity,
               int InletProfile,   int CollisionModel, int CurvedBoundaries,
               int OutletProfile,  int Iterations,     int AutosaveAfter,
               int AutosaveEvery,  int postproc_prog,  int CalculateDragLift,
               float ConvergenceCritVeloc, float ConvergenceCritRho)
{


    clock_t tStart;     // Time measurement: declaration, begin
    if(MYTHREAD==0)
        tStart = clock(); // BEGIN OF OVERALL TIME MEASUREMENT

    ////////////////////////////////////////////////////
    ///////////////////// Declare //////////////////////
    ////////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    /////////// LOCAL VARIABLES FOR EACH CPU ///////////
    ////////////////////////////////////////////////////

    int i, j, k, iter = 0;                  // variables for loops
    FILE* resid_file;                       // file for residuals
    FILE* log_file;                         // file for log
    FILE* TimeMeasurementFile;              // file for time measurement results
    char IterOutputFile[50];                // write results to this file after the iterations
    char AutosaveOutputFile[50];            // autosave filename
    char OutputFile[50];                    // initial data will be written to this file
    char FinalOutputFile[50];               // final data will be written to this file
    char logFile[] = "Results/logFile.log"; // path of the .log file
    int  AutosaveI = 1;                     // autosave i variable, will be incremented after every autosave
    int* ppp;                               // pointer of the postproc_prog variable

    CellProps *Cells;                // Struct for Cells

    // Time measurement variables
    float tInitialization  = 0.0; // Time measurement of Initialization
    float tIteration       = 0.0; // Time measurement of Iteration
    float tCollision       = 0.0; // Time measurement of Collision
    float tUpdateF         = 0.0; // Time measurement of UpdateF
    float tStreaming       = 0.0; // Time measurement of Streaming
    float tBoundaries      = 0.0; // Time measurement of Boundaries
    float tUpdateMacro     = 0.0; // Time measurement of Update Macroscopic vars
    float tResiduals       = 0.0; // Time measurement of calculating residuals
    float tWriting         = 0.0; // Time measurement of writing data
    float tBCells          = 0.0; // Time measurement of handling boundaries
    clock_t tInstant1, tInstant2; // Time measurement points: universal
    clock_t tIterStart, tIterEnd; // Time measurement points: main loop

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

    ////////////////////////////////////////////////////
    //////////////////// ALLOCATE //////////////////////
    ////////////////////////////////////////////////////

    // allocate residuals
    sumVel0   = Create1DArrayDouble(1);
    sumVel1   = Create1DArrayDouble(1);
    sumRho0   = Create1DArrayDouble(1);
    sumRho1   = Create1DArrayDouble(1);
    Residuals = Create1DArrayDouble(4);

    ////////////////////////////
    // THESE ARE SHARED STUFF //
    ////////////////////////////

    // allocate mesh properties  :: DEFINED IN ShellFunctions.h
    Delta          = (shared double*)upc_alloc(1*sizeof(double));
    m              = (shared int*)upc_alloc(1*sizeof(int));
    n              = (shared int*)upc_alloc(1*sizeof(int));
    o              = (shared int*)upc_alloc(1*sizeof(int));
    NumNodes       = (shared int*)upc_alloc(1*sizeof(int));
    NumConn        = (shared int*)upc_alloc(1*sizeof(int));
    MaxInletCoordY = (shared double*)upc_alloc(1*sizeof(double));
    MinInletCoordY = (shared double*)upc_alloc(1*sizeof(double));
    NumInletNodes  = (shared int*)upc_alloc(1*sizeof(int));

    ///////////////////////////////////////////////////////////////////////
    //// WE HAVE BETTER TO STORE THIS ON EACH CPU, MAKE THINGS EASIER  ////
    ///////////////////////////////////////////////////////////////////////

    // D2Q9 Variables of the lattice
    w  = Create1DArrayDouble(19); // weight values for the directions
    c  = Create1DArrayInt(19);    // 
    cx  = Create1DArrayInt(19);    // x coordinate of the discrete lattice directions
    cy  = Create1DArrayInt(19);    // y coordinate of the discrete lattice directions
    cz  = Create1DArrayInt(19);    // z coordinate of the discrete lattice directions
    opp = Create1DArrayInt(19);    // opposite vector

    ////////////////////////////////////////////////////
    ///////////////////// Read data ////////////////////
    ////////////////////////////////////////////////////

    Nodes   = ReadNodes(NodeDataFile);          // Read Node data
    BCconn  = ReadBCconn(BCconnectorDataFile);  // Read BCconn data
    CompDataNode(Nodes);
    //CompDataConn(BCconn);

    ////////////////////////////////////////////////////
    /////////////// Print info to log //////////////////
    ////////////////////////////////////////////////////

    if(MYTHREAD==0) // Print information to log file
    {
        // Check whether we got back what we wanted :), write to log file!
        log_file = fopen(logFile, "w");  // open log file
        ppp      = &postproc_prog;       // for convenience ppp points to postproc_prog
        fprintf(log_file,"This is the 3D lattice Boltzmann *.log file\n\n");
        fprintf(log_file,"\n:::: Imported variables from the *.ini file :::: \n");
        fprintf(log_file,">>> Uavg              : %3.6f\n", Uavg);
        fprintf(log_file,">>> Vavg              : %3.6f\n", Vavg);
        fprintf(log_file,">>> Wavg              : %3.6f\n", Wavg);
        fprintf(log_file,">>> Initial density   : %2.1f\n", rho_ini);
        fprintf(log_file,">>> Viscosity         : %3.8f\n", Viscosity);
        fprintf(log_file,">>> # of iterations   : %1.1d\n", Iterations);
        fprintf(log_file,">>> Autosave after    : %1.1d\n", AutosaveAfter);
        fprintf(log_file,">>> Autosave every    : %1.1d\n", AutosaveEvery);
        fprintf(log_file,">>> Convergence Veloc : %3.8f\n", ConvergenceCritVeloc);
        fprintf(log_file,">>> Convergence Rho   : %3.8f\n", ConvergenceCritRho);
        switch(CollisionModel)         // 1: BGKW, 2: TRT, 3: MRT
        {
            case 1: fprintf(log_file,">>> CollisionModel    : BGKW\n"); break;
            case 2: fprintf(log_file,">>> CollisionModel    : TRT\n" ); break;
            case 3: fprintf(log_file,">>> CollisionModel    : MRT\n" ); break;
        }
        switch(InletProfile)                      // 1:ON, 2:OFF
        {
            case 1: fprintf(log_file,">>> InletProfile      : ON\n" ); break;
            case 2: fprintf(log_file,">>> InletProfile      : OFF\n"); break;
        }
        switch(OutletProfile)                     // 1:ON, 2:OFF
        {
            case 1: fprintf(log_file,">>> OutletProfile     : ON\n" ); break;
            case 2: fprintf(log_file,">>> OutletProfile     : OFF\n"); break;
        }
        switch(CurvedBoundaries)                  // 1:ON, 2:OFF
        {
            case 1: fprintf(log_file,">>> CurvedBoundaries  : ON\n" ); break;
            case 2: fprintf(log_file,">>> CurvedBoundaries  : OFF\n"); break;
        }
        switch(postproc_prog)   // 1->Paraview (*.csv)     2->Tecplot
        {
            case 1: fprintf(log_file,">>> Results format    : Paraview (*.csv)\n" ); break;
            case 2: fprintf(log_file,">>> Results format    : Tecplot (*.dat)\n"); break;
        }
        if (CalculateDragLift != 0)
            fprintf(log_file,">>> Drag, lift @ BC   : %d\n", CalculateDragLift);
        else
            fprintf(log_file,">>> Drag, lift was not calculated\n");

        fprintf(log_file,"\n:::: Calculated variables from mesh :::: \n");
        fprintf(log_file,">>> Grid spacing        = %f\n", *Delta);
        fprintf(log_file,">>> # of nodes in x (n) = %d\n", *n);
        fprintf(log_file,">>> # of nodes in y (m) = %d\n", *m);
        fprintf(log_file,">>> NumInletNodes       = %d\n", *NumInletNodes);
        fprintf(log_file,">>> MaxInletCoordY      = %f\n", *MaxInletCoordY);
        fprintf(log_file,">>> MinInletCoordY      = %f\n", *MinInletCoordY);

        fprintf(log_file,"\n:::: Parallel properties :::: \n");
        fprintf(log_file,">>> # of threads        = %d\n", THREADS);
        fprintf(log_file,">>> BlockSize           = %d\n", BLOCKSIZE);

        // In case of no autosave
        sprintf(AutosaveOutputFile, "NOWHERE!");

    } // END OF THREAD ZERO

    ////////////////////////////////////////////////////
    ///////////////// INITIALIZE ///////////////////////
    ////////////////////////////////////////////////////

    // Fill up D3Q19 variables
    D3Q19Vars(w, cx, cy, cz, opp, c);


    //////////////////////////////////////////////////////
    // Allocate structure for the cell properties (see ShellFunctions.h)
    Cells = calloc(BLOCKSIZE+2*LAYER,sizeof(CellProps));
    //////////////////////////////////////////////////////

    if(MYTHREAD==0)
    {
        fprintf(log_file,"\n:::: Initializing ::::\n");
        printf("\n:::: Initializing ::::\n");
    } // END OF THREAD ZERO

    ////////////////////////////////////////////////////
    ///////////////// INITIALIZE CELLS /////////////////
    ////////////////////////////////////////////////////

    upc_barrier;         // Synchronise
    //loop = 0;            // This will measure that how much is done of initialization
    tInstant1 = clock(); // Measure time of initialization

    CellIni( Cells,
             Nodes,            // Nodes
             BCconn,           // BCconn
             Uavg,             // INPUT PARAMETER
             Vavg,             // INPUT PARAMETER
             Wavg,             // INPUT PARAMETER
             InletProfile,     // INPUT PARAMETER
             CollisionModel,   // INPUT PARAMETER
             opp,              // Opposite direction
             rho_ini);         // Initial density



    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /* FALTA IMPRIMIR TODO LO QUE PASA EN ESTA FUNCIÓN A LO GRANDE, COMO SI NO HUBIESE UN MAÑANA */
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    tInstant2 = clock(); // Measure time of initialization
    tInitialization = (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;


    //I will start from the end

    tInstant1 = clock(); // Start measuring time
    ComputeResiduals(Cells, Residuals,
                     sumVel0, sumVel1,
                     sumRho0, sumRho1,
                     CalculateDragLift, &iter, &Iterations);
    //ComputeResiduals(Cells, Residuals, sumVel0, sumVel1, sumRho0, sumRho1, CalculateDragLift, &iter, &Iterations);
    tInstant2 = clock(); // End of time measuremet
    tResiduals = tResiduals + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;



    upc_barrier;         // Synchronise

    upc_free(Delta);
    upc_free(m);
    upc_free(n);
    upc_free(MaxInletCoordY);
    upc_free(MinInletCoordY);
    upc_free(NumInletNodes);
    upc_free(NumNodes);
    upc_free(NumConn);

    free(Cells);
    free(w);
    free(cx);
    free(cy);
    free(cz);
    free(c);
    free(opp);
}
