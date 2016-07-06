#include <stdio.h>                      // printf();
#include <stdlib.h>                     // for calloc();
#include <time.h>                       // time functions
#include <upc.h>                // Required for UPC
#include <CellFunctions.h>
#include <ComputeResiduals.h>
#include <ShellFunctions.h>

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



    //Things to print for each thread
    //Point coordinates and values
    //Cells contain, all the nodes, on the different threads.
    upc_barrier;
    //if(MYTHREAD == 0)
        //out_cells_file = shared fopen("Results/outCells.dat","a");

    upc_barrier;
    print_cells_info(Cells);
    upc_barrier;
    print_boundary_type(Cells);
    upc_barrier;

    //For every boundary type, print which points correspond to them.

    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /* FALTA IMPRIMIR TODO LO QUE PASA EN ESTA FUNCIÓN A LO GRANDE, COMO SI NO HUBIESE UN MAÑANA */
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    tInstant2 = clock(); // Measure time of initialization
    tInitialization = (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;


    //I will start from the end

    tInstant1 = clock(); // Start measuring time
    /*ComputeResiduals(Cells, Residuals,
                     sumVel0, sumVel1,
                     sumRho0, sumRho1,
                     CalculateDragLift, &iter, &Iterations);*/
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

void print_cells_info(CellProps* Cells) {
    upc_barrier;
    FILE* out_cells_file;
    char out_cells_filename [50];
    sprintf(out_cells_filename,"Results/outCells/T_%i.dat",MYTHREAD);
    out_cells_file = fopen(out_cells_filename,"w");
    for (int t_to_print_from = 0; t_to_print_from<THREADS;t_to_print_from++){
        if (MYTHREAD == t_to_print_from)
        {
            printf("Going to Thread %i\n",MYTHREAD);
            fprintf(out_cells_file,"Printing from thread %i\n",MYTHREAD);
            fprintf(out_cells_file,"   ID   |  i |  j |  k ||    x    |    y    |    z\n");
            for(int cell_to_print = LAYER; cell_to_print< BLOCKSIZE+LAYER; cell_to_print++) {
                print_cell_line(out_cells_file,Cells+cell_to_print);
            }
        }
        upc_barrier;
    }
    fclose(out_cells_file);
    upc_barrier;
}

void print_boundary_type(CellProps* Cells) {
    main_thread
        printf("printing boundaries T%i\n",MYTHREAD);


    int count_B[4];
    count_B[0]=0;
    count_B[1]=0;
    count_B[2]=0;
    count_B[3]=0;

    int* N_B[4];

    N_B[0] = (int*)malloc(BLOCKSIZE*sizeof(int));
    N_B[1] = (int*)malloc(BLOCKSIZE*sizeof(int));
    N_B[2] = (int*)malloc(BLOCKSIZE*sizeof(int));
    N_B[3] = (int*)malloc(BLOCKSIZE*sizeof(int));

    upc_barrier;
    main_thread
        printf("Get cell boundary info\n");
    for (int node_to_look = LAYER;node_to_look<BLOCKSIZE+LAYER;node_to_look++) {
        int BT;
        //printf("TEST, thread %i\n",MYTHREAD);
        if (node_to_look == LAYER)    {
            printf("Thread: %i,Node: %i,BT: %i\n",MYTHREAD,node_to_look,(Cells+node_to_look)->Boundary);
            int index_n, index_i, index_j, index_k;
            index_n = (Cells+LAYER)->ID;
            index_i = index_n/(NM*NL);
            index_j = (index_n - index_i * NM* NL)/NM;
            index_k = (index_n - index_i * NM* NL - index_j * NM);
            printf("%7i |%3i |%3i |%3i || %5.5f | %5.5f | %5.5f \n",
                    index_n,
                    index_i,
                    index_j,
                    index_k,
                    (Cells+LAYER)->CoordX,
                    (Cells+LAYER)->CoordY,
                    (Cells+LAYER)->CoordZ
            );
        }



        if ((BT = (Cells+node_to_look)->Boundary-1) !=-1 ) {
            //printf("BT: %i\n",BT);
            N_B[BT][count_B[BT]] = (Cells+node_to_look)->ID;
            count_B[BT]++;
        }
    }
    upc_barrier;
    main_thread
        printf("Start saving to files\n");

    char* b_filename[4];
    char* b1_filename = "Results/boundary/solidplane_boundary.dat";
    char* b2_filename = "Results/boundary/fluidplane_boundary.dat";
    char* b3_filename = "Results/boundary/edge_boundary.dat";
    char* b4_filename = "Results/boundary/corner_boundary.dat";
    b_filename[0] = b1_filename;
    b_filename[1] = b2_filename;
    b_filename[2] = b3_filename;
    b_filename[3] = b4_filename;
    if (MYTHREAD == 0) {
    printf(b_filename[0]);
    printf("\n");
    printf(b1_filename);
    printf("\n");
    printf(b_filename[1]);
    printf("\n");
    printf(b2_filename);
    printf("\n");
    printf(b_filename[2]);
    printf("\n");
    printf(b3_filename);
    printf("\n");
    printf(b_filename[3]);
    printf("\n");
    printf(b4_filename);
    printf("\n");}

    if (MYTHREAD == 0){
        for(int BT = 0; BT<4; BT++) {
            FILE* b_file = fopen(b_filename[BT],"w");
            fclose(b_file);
        }
    }

    upc_barrier;
    for (int t_to_print_from = 0; t_to_print_from<THREADS;t_to_print_from++){
        if (MYTHREAD == t_to_print_from) {

            for(int BT = 0; BT<4; BT++) {
                FILE* b_file = fopen(b_filename[BT],"a");
                fprintf(b_file,"Printing from thread %i\n",MYTHREAD);
                fprintf(b_file,"   ID   |  i |  j |  k ||    x    |    y    |    z\n");
                for(int c=0;c<count_B[BT];c++){
                    print_cell_line(b_file,cell_from_id(Cells,N_B[BT][c]));
                }
                fclose(b_file);
            }
        }
        upc_barrier;
    }



    free(N_B[0]);
    free(N_B[1]);
    free(N_B[2]);
    free(N_B[3]);

}

void print_cell_line(FILE* file, const CellProps* Cell) {
    int index_n, index_i, index_j, index_k;
    index_n = Cell->ID;
    index_i = index_n/(NM*NL);
    index_j = (index_n - index_i * NM* NL)/NM;
    index_k = (index_n - index_i * NM* NL - index_j * NM);
    fprintf(file,"%7i |%3i |%3i |%3i || %5.5f | %5.5f | %5.5f \n",
            index_n,
            index_i,
            index_j,
            index_k,
            Cell->CoordX,
            Cell->CoordY,
            Cell->CoordZ
    );
}

CellProps* cell_from_id(CellProps* Cells, int ID){
    int n = ID + LAYER - MYTHREAD*BLOCKSIZE;
    return (Cells+n);
}