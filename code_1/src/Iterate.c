#include <stdio.h>                      // printf();
#include <stdlib.h>                     // for calloc();
#include <time.h>                       // time functions
#include <upc.h>                // Required for UPC
#include <CellFunctions.h>
#include <ComputeResiduals.h>
#include <ShellFunctions.h>
#include <FilesWriting.h>

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

void time_meas_vars_init();

void allocate_residuals();

void allocate_lattice_vars();

void allocate_shared();

void allocate_vars();

void read_data(const char *NodeDataFile, const char *BCconnectorDataFile);

void print_init_info_to_log(float Uavg, float Vavg, float Wavg, float rho_ini, float Viscosity, int InletProfile,
                            int CollisionModel, int CurvedBoundaries, int OutletProfile, int Iterations,
                            int AutosaveAfter, int AutosaveEvery, int postproc_prog, int CalculateDragLift,
                            float ConvergenceCritVeloc, float ConvergenceCritRho);

void free_vars();

void init_vars(int *postproc_prog);

void calc_collision_freq(float Viscosity);

void alloc_cells();

void export_data(int postproc_prog, const char *fnMemCopyRes);

void main_while_loop(int CollisionModel, int CurvedBoundaries, int OutletProfile, int *Iterations, int AutosaveAfter,
                     int AutosaveEvery, int postproc_prog, int CalculateDragLift, float ConvergenceCritVeloc,
                     float ConvergenceCritRho);

void save_init_data(int postproc_prog);

void write_cells_to_results(int postproc_prog);

void free_mesh_data_matrices();

void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg, float Vavg, float Wavg,
               float rho_ini, float Viscosity,
               int InletProfile, int CollisionModel, int CurvedBoundaries,
               int OutletProfile, int Iterations, int AutosaveAfter,
               int AutosaveEvery, int postproc_prog, int CalculateDragLift,
               float ConvergenceCritVeloc, float ConvergenceCritRho)
{


    main_thread
        tStart = clock(); // BEGIN OF OVERALL TIME MEASUREMENT



    init_vars(&postproc_prog);

    read_data(NodeDataFile, BCconnectorDataFile);

    print_init_info_to_log(Uavg, Vavg, Wavg, rho_ini, Viscosity, InletProfile, CollisionModel, CurvedBoundaries,
                           OutletProfile, Iterations, AutosaveAfter, AutosaveEvery, postproc_prog, CalculateDragLift,
                           ConvergenceCritVeloc, ConvergenceCritRho);

    ////////////////////////////////////////////////////
    ///////////////// INITIALIZE ///////////////////////
    ////////////////////////////////////////////////////

    // Fill up D3Q19 variables
    D3Q19Vars(w, cx, cy, cz, opp, c);
    calc_collision_freq(Viscosity);


    // Initialize variables for MRT Collision model, if used
    tm    = Create2DArrayDouble(9, 9);
    stmiv = Create2DArrayDouble(9, 9);

    if (CollisionModel == 3)
        MRTInitializer(tm, stmiv, Omega);


    // PUT THREAD INFO TO BCELLS (boundary cells)

    upc_forall(i = 0; i < 2*LAYER*THREADS; i++; &BCells[i])
    { (BCells+i)->ThreadNumber = MYTHREAD; }

    if(MYTHREAD==0)
    {
        fprintf(log_file,"\n:::: Initializing ::::\n");
        printf("\n:::: Initializing ::::\n");
    } // END OF THREAD ZERO

    ////////////////////////////////////////////////////
    ///////////////// INITIALIZE CELLS /////////////////
    ////////////////////////////////////////////////////

    upc_barrier;         // Synchronise

    init_measure_time;
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
    end_measure_time(tCellsInitialization);

    upc_barrier;
    print_cells_info(Cells);
    upc_barrier;
    print_boundary_type(Cells);
    upc_barrier;

    end_measure_time(tInitialization);
    //tInstant2 = clock(); // Measure time of initialization
    //tInitialization = (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

    upc_barrier;         // Synchronise

    // PUT INITIALIZED DATA TO BOUNDARIES
    putCellsToShared(Cells);

    // COPY CELLS TO WCELLS (writing cells) TO WRITE DATA
    putCellsToWCells(Cells);
    char fnMemCopyRes[50];


    write_cells_to_results(postproc_prog);


    free_mesh_data_matrices();  //BCNode

    upc_barrier;         // Synchronise

    ////////////////////////////////////////////////////
    /////////////// SAVE INITIALIZATION ////////////////
    ////////////////////////////////////////////////////

    if(MYTHREAD==0)
    {
        fprintf(log_file,"\n:::: Initialization done! ::::\n");
        printf("\n:::: Initialization done! ::::\n");


        save_init_data(postproc_prog);

        // Open residuals file
        resid_file = fopen("Results/residuals.dat", "w");
        fprintf(resid_file,"Iter Time Vel_res Rho_res Drag Lift\n");

        fprintf(log_file,"\n:::: Start Iterations ::::\n");
        printf("\n:::: Start Iterations ::::\n");
    } // END OF THREAD ZERO

    ////////////////////////////////////////////////////
    /////////////////// ITERATE ////////////////////////
    ////////////////////////////////////////////////////

    Residuals[0] = 100;
    Residuals[1] = 100;

    upc_barrier;         // Synchronise

    tIterStart = clock(); // Start measuring time of main loop

    main_while_loop(CollisionModel, CurvedBoundaries, OutletProfile, &Iterations, AutosaveAfter, AutosaveEvery,
                    postproc_prog, CalculateDragLift, ConvergenceCritVeloc, ConvergenceCritRho);

    tIterEnd = clock(); // End measuring time of main loop
    tIteration = (float)(tIterEnd - tIterStart ) / CLOCKS_PER_SEC;

    upc_barrier;         // Synchronise

    export_data(postproc_prog, fnMemCopyRes);


    ////////////////////////////////////////////////////
    ///////////////// End of line //////////////////////
    ////////////////////////////////////////////////////

    upc_barrier;         // Synchronise

    free_vars();

}

void free_mesh_data_matrices() {

// We dont need these matrices anymore
    free(Nodes);
    free(BCconn);
}

void write_cells_to_results(int postproc_prog) {

// Write boundary cells to Results to see how mesh was distributed
    char  fnMemCopyRes[50];
    if(MYTHREAD==0)
    {
        switch(postproc_prog)
        { case 1: sprintf(fnMemCopyRes, "Results/MyBoundaryCells.csv");   break;
            case 2: sprintf(fnMemCopyRes, "Results/MyBoundaryCells.dat");   break; }

        WriteBCells(fnMemCopyRes, ppp);
    }
}

void save_init_data(int postproc_prog) {// Write Initialized data
    switch(postproc_prog)
        { case 1: sprintf(OutputFile, "Results/InitialData.csv");   break;
            case 2: sprintf(OutputFile, "Results/InitialData.dat");   break; }

    WriteResults(OutputFile, ppp);
    printf("\nInitialized data was written to %s\n", OutputFile);
}

void main_while_loop(int CollisionModel, int CurvedBoundaries, int OutletProfile, int *Iterations, int AutosaveAfter,
int AutosaveEvery, int postproc_prog, int CalculateDragLift, float ConvergenceCritVeloc,
float ConvergenceCritRho) {
    while (Residuals[0] > ConvergenceCritVeloc && Residuals[1] > ConvergenceCritRho && iter < (*Iterations))
{

//////////////// COLLISION ////////////////
tInstant1 = clock();
CollisionStep(Cells, w, cx, cy,  opp,  Omega,  OmegaA, tm, stmiv, CollisionModel); ////////////////////// !!!!!!!!!!!!!!!!! CX CY!
tInstant2 = clock();
tCollision = tCollision + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

////////////// UPDATE DISTR ///////////////
tInstant1 = clock(); // Start measuring time
for(i = LAYER;  i < LAYER + BLOCKSIZE;  i++)
{ UpdateF(Cells, i); }
tInstant2 = clock();
tUpdateF = tUpdateF + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

// PUT THREAD-BOUNDARY CELLS TO SHARED
tInstant1 = clock(); // Start measuring time
putCellsToShared(Cells);
tInstant2 = clock(); // End of time measuremet
tBCells = tBCells + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

upc_barrier;         // Synchronise

//////////////// COPY SHARED BCELLS TO CELLS ////////////////
tInstant1 = clock();
getSharedToCells(Cells);
tInstant2 = clock();
tBCells = tBCells + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

////////////// STREAMING ///////////////
tInstant1 = clock(); // Start measuring time
StreamingStep(Cells, c);
tInstant2 = clock();
tStreaming = tStreaming + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

////////////// BOUNDARIES ///////////////
tInstant1 = clock(); // Start measuring time
HandleBoundariesStep(Cells, cx, cy, cz, c, opp, OutletProfile, CurvedBoundaries);
tInstant2 = clock();
tBoundaries = tBoundaries + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

// UPDATE VELOCITY AND DENSITY
tInstant1 = clock(); // Start measuring time
UpdateMacroscopicStep(Cells, cx, cy, CalculateDragLift);
tInstant2 = clock();
tUpdateMacro = tUpdateMacro + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

////////////// Residuals ///////////////
tInstant1 = clock(); // Start measuring time
ComputeResiduals(Cells, Residuals, sumVel0, sumVel1, sumRho0, sumRho1, CalculateDragLift, &iter, Iterations);
//ComputeResiduals(Cells, Residuals, sumVel0, sumVel1, sumRho0, sumRho1, CalculateDragLift, &iter, &Iterations);
tInstant2 = clock(); // End of time measuremet
tResiduals = tResiduals + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

if(MYTHREAD==0)
fprintf(resid_file,"%d %5.4e %5.4e %5.4e %f %f\n", iter, (iter+1.0)*(*Delta), Residuals[0], Residuals[1], Residuals[2], Residuals[3]);

iter++; // update loop variable

if(iter%100==0 && MYTHREAD==0){
printf("Iterations: %05d/%05d || ", iter, (*Iterations));
printf("Residuals: l2norm  %e; L2_norm_weighted  %e\n", Residuals[0], Residuals[1]);
}

////////////// Autosave ///////////////
if(iter == (AutosaveEvery*AutosaveI))
{
AutosaveI++;
if(iter>AutosaveAfter)
{
tInstant1 = clock(); // Start measuring time
switch(postproc_prog) {
case 1: sprintf(AutosaveOutputFile, "Results/autosave_iter%05d.csv", iter); break;
case 2: sprintf(AutosaveOutputFile, "Results/autosave_iter%05d.dat", iter); break; }
putCellsToWCells(Cells); // Put information to WCells and write (Write Cells)
if (MYTHREAD==0) // AUTOSAVE
WriteResults(AutosaveOutputFile, ppp);
tInstant2 = clock(); // End of time measurement
tWriting = tWriting + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;
}
}
//////////////////////////////////////////////////////
}     ////////////// END OF MAIN WHILE CYCLE ///////////////
//////////////////////////////////////////////////////

}

void export_data(int postproc_prog, const char *fnMemCopyRes) {
    if(MYTHREAD == 0) // EXPORT DATA, TIME MEASUREMENT RESULTS
    {
        // Close residuals file
        fclose(resid_file);

        clock_t tEnd = clock();
        float tOverall = (float)(tEnd - tStart) / CLOCKS_PER_SEC; // Calculate elapsed time

        fprintf(log_file,"\nOverall calculations took %f seconds\n", tOverall);
        fprintf(log_file,"Main while loop took %f seconds\n",        tIteration);
        fprintf(log_file,"Initialization took %f seconds\n",         tInitialization);
        fprintf(log_file,"Collision took %f seconds\n",              tCollision);
        fprintf(log_file,"UpdateF took %f seconds\n",                tUpdateF);
        fprintf(log_file,"Streaming took %f seconds\n",              tStreaming);
        fprintf(log_file,"Calculating Boundaries took %f seconds\n", tBoundaries);
        fprintf(log_file,"Update Macroscopic took %f seconds\n",     tUpdateMacro);
        fprintf(log_file,"Calculating Residuals took %f seconds\n",  tResiduals);
        fprintf(log_file,"Writing results took %f seconds\n",        tWriting);
        fprintf(log_file,"Copying boundary cells took %f seconds\n", tBCells);


        // end time measurement, close log file
        fprintf(log_file,"\n:::: Iterations done! ::::\n");
        fclose(log_file);

        TimeMeasurementFile = fopen("Results/ParallelTimeMeasuerment.dat","w");
        fprintf(TimeMeasurementFile,"tOverall %f\n",        tOverall);
        fprintf(TimeMeasurementFile,"tIteration %f\n",      tIteration);
        fprintf(TimeMeasurementFile,"tInitialization %f\n", tInitialization);
        fprintf(TimeMeasurementFile,"tCollision %f\n",      tCollision);
        fprintf(TimeMeasurementFile,"tUpdateF %f\n",        tUpdateF);
        fprintf(TimeMeasurementFile,"tStreaming %f\n",      tStreaming);
        fprintf(TimeMeasurementFile,"tBoundaries %f\n",     tBoundaries);
        fprintf(TimeMeasurementFile,"tUpdateMacro %f\n",    tUpdateMacro);
        fprintf(TimeMeasurementFile,"tResiduals %f\n",      tResiduals);
        fprintf(TimeMeasurementFile,"tWriting %f\n",        tWriting);
        fprintf(TimeMeasurementFile,"tBCells %f\n",         tBCells);
        fprintf(TimeMeasurementFile,"THREADS %d\n",         THREADS);
        fclose(TimeMeasurementFile);

        // Write final data
        switch(postproc_prog){
            case 1: sprintf(FinalOutputFile, "Results/FinalData.csv"); break;
            case 2: sprintf(FinalOutputFile, "Results/FinalData.dat"); break;
        }
        WriteResults(FinalOutputFile,  ppp);

        // Write information for user
        printf("\n\nLog was written to %s\n", logFile);
        printf("Last autosave result can be found at %s\n", AutosaveOutputFile);
        printf("Residuals were written to Results/residuals.dat\n");
        printf("Profiling results were written to Results/ParallelTimeMeasuerment.dat\n");
        printf("Final results were written to %s\n", FinalOutputFile);

        WriteBCells(fnMemCopyRes, ppp);
        puts("BCells were written!");
    } // END OF THREAD ZERO


}

void alloc_cells() {//////////////////////////////////////////////////////
    // Allocate structure for the cell properties (see ShellFunctions.h)
    WCells = (shared_block(BLOCKSIZE)   CellProps*)upc_all_alloc(THREADS, BLOCKSIZE*sizeof(CellProps));
    BCells = (shared_block(2*LAYER)     CellProps*)upc_all_alloc(THREADS,     LAYER*sizeof(CellProps));
    Cells = calloc(BLOCKSIZE+2*LAYER,sizeof(CellProps));
    //////////////////////////////////////////////////////

}

void calc_collision_freq(float Viscosity) {// Calculate collision frequency
    Omega  = 1.0/(3.*Viscosity+0.5);
    OmegaA = 8*(2-Omega)/(8-Omega);
}

void init_vars(int *postproc_prog) {
    iter = 0;                               // variables for loops
    sprintf(logFile,"Results/logFile.log"); // path of the .log file
    AutosaveI = 1;                          // autosave i variable, will be incremented after every autosave
    ppp      = postproc_prog;       // for convenience ppp points to postproc_prog

    time_meas_vars_init();
    allocate_vars();
}

void free_vars() {
    upc_free(Delta);
    upc_free(m);
    upc_free(n);
    //upc_free(MaxInletCoordY);
    //upc_free(MinInletCoordY);
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

void print_init_info_to_log(float Uavg, float Vavg, float Wavg, float rho_ini, float Viscosity, int InletProfile,
                            int CollisionModel, int CurvedBoundaries, int OutletProfile, int Iterations,
                            int AutosaveAfter, int AutosaveEvery, int postproc_prog, int CalculateDragLift,
                            float ConvergenceCritVeloc, float ConvergenceCritRho) {
////////////////////////////////////////////////////
/////////////// Print info to log //////////////////
////////////////////////////////////////////////////

    if(MYTHREAD==0) // Print information to log file
    {
// Check whether we got back what we wanted :), write to log file!
        log_file = fopen(logFile, "w");  // open log file
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
//fprintf(log_file,">>> MaxInletCoordY      = %f\n", *MaxInletCoordY);
//fprintf(log_file,">>> MinInletCoordY      = %f\n", *MinInletCoordY);

        fprintf(log_file,"\n:::: Parallel properties :::: \n");
        fprintf(log_file,">>> # of threads        = %d\n", THREADS);
        fprintf(log_file,">>> BlockSize           = %d\n", BLOCKSIZE);

// In case of no autosave
        sprintf(AutosaveOutputFile, "NOWHERE!");

    } // END OF THREAD ZERO

}

void read_data(const char *NodeDataFile, const char *BCconnectorDataFile) {////////////////////////////////////////////////////
    ///////////////////// Read data ////////////////////
    ////////////////////////////////////////////////////

    Nodes   = ReadNodes(NodeDataFile);          // Read Node data
    BCconn  = ReadBCconn(BCconnectorDataFile);  // Read BCconn data
    CompDataNode(Nodes);
    //CompDataConn(BCconn);

}

void allocate_vars() {////////////////////////////////////////////////////
    //////////////////// ALLOCATE //////////////////////
    ////////////////////////////////////////////////////

    allocate_residuals();
    allocate_shared();
    allocate_lattice_vars();
    alloc_cells();
}

void allocate_shared() {////////////////////////////
    // THESE ARE SHARED STUFF //
    ////////////////////////////

    // allocate mesh properties  :: DEFINED IN ShellFunctions.h
    Delta          = (shared double*)upc_alloc(1*sizeof(double));
    m              = (shared int*)upc_alloc(1*sizeof(int));
    n              = (shared int*)upc_alloc(1*sizeof(int));
    o              = (shared int*)upc_alloc(1*sizeof(int));
    NumNodes       = (shared int*)upc_alloc(1*sizeof(int));
    NumConn        = (shared int*)upc_alloc(1*sizeof(int));
    //MaxInletCoordY = (shared [] double*)upc_alloc(1*sizeof(double));--------------------------> DO I NEED THIS?
    //MinInletCoordY = (shared [] double*)upc_alloc(1*sizeof(double));--------------------------> DO I NEED THIS?
    NumInletNodes  = (shared int*)upc_alloc(1*sizeof(int));
}

void allocate_lattice_vars() {// D2Q9 Variables of the lattice
    w  = Create1DArrayDouble(19); // weight values for the directions
    c  = Create1DArrayInt(19);    //
    cx  = Create1DArrayInt(19);    // x coordinate of the discrete lattice directions
    cy  = Create1DArrayInt(19);    // y coordinate of the discrete lattice directions
    cz  = Create1DArrayInt(19);    // z coordinate of the discrete lattice directions
    opp = Create1DArrayInt(19);    // opposite vector

}

void allocate_residuals() {// allocate residuals
    sumVel0   = Create1DArrayDouble(1);
    sumVel1   = Create1DArrayDouble(1);
    sumRho0   = Create1DArrayDouble(1);
    sumRho1   = Create1DArrayDouble(1);
    Residuals = Create1DArrayDouble(4);
}

void time_meas_vars_init() {// Time measurement variables
    tInitialization  = 0.0; // Time measurement of Initialization
    tIteration       = 0.0; // Time measurement of Iteration
    tCollision       = 0.0; // Time measurement of Collision
    tUpdateF         = 0.0; // Time measurement of UpdateF
    tStreaming       = 0.0; // Time measurement of Streaming
    tBoundaries      = 0.0; // Time measurement of Boundaries
    tUpdateMacro     = 0.0; // Time measurement of Update Macroscopic vars
    tResiduals       = 0.0; // Time measurement of calculating residuals
    tWriting         = 0.0; // Time measurement of writing data
    tBCells          = 0.0; // Time measurement of handling boundaries

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