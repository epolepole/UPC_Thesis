#include <stdio.h>                      // printf();
#include <stdlib.h>                     // for calloc();
#include <time.h>                       // time functions
#include <CellFunctions.h>
#include <ComputeResiduals.h>
#include <FilesWriting.h>
#include <ShellFunctions.h>

////////////////////////////////////////////////////
////////////////// OWN HEADERS /////////////////////
////////////////////////////////////////////////////

#include "Iterate.h"            // Iteration takes place
#include "FilesReading.h"       // For reading files
#include "tests.h"


////////////////////////////////////////////////////
/////////////////// ITERATION //////////////////////
////////////////////////////////////////////////////




void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg, float Vavg, float Wavg,
               float rho_ini, float Viscosity,
               int InletProfile, int CollisionModel, int CurvedBoundaries,
               int OutletProfile, int Iterations, int AutosaveAfter,
               int AutosaveEvery, int postproc_prog, int CalculateDragLift,
               float ConvergenceCritVeloc, float ConvergenceCritRho)
{


    if(MYTHREAD == 0) {
        printf("Compiler defined parameters:\n");
        printf("NN = %i\n", NN);
        printf("NM = %i\n", NM);
        printf("NL = %i\n", NL);
        printf("BLOCKSIZE = %i\n", BLOCKSIZE);
    }

    main_thread
        tStart = clock(); // BEGIN OF OVERALL TIME MEASUREMENT


    init_tests();
    init_vars(&postproc_prog);

    read_data(NodeDataFile, BCconnectorDataFile);

    print_init_info_to_log(Uavg, Vavg, Wavg, rho_ini, Viscosity, InletProfile, CollisionModel, CurvedBoundaries,
                           OutletProfile, Iterations, AutosaveAfter, AutosaveEvery, postproc_prog, CalculateDragLift,
                           ConvergenceCritVeloc, ConvergenceCritRho);

    ////////////////////////////////////////////////////
    ///////////////// INITIALIZE ///////////////////////
    ////////////////////////////////////////////////////

    // Fill up D3Q19 variables
    D3Q19Vars();
    calc_collision_freq(Viscosity);


    // Initialize variables for MRT Collision model, if used
    tm    = Create2DArrayDouble(9, 9);
    stmiv = Create2DArrayDouble(9, 9);

    if (CollisionModel == 3)
        MRTInitializer(tm, stmiv, Omega);


    // PUT THREAD INFO TO BCELLS (boundary cells)

    /*upc_forall(i = 0; i < 2*LAYER*THREADS; i++; &BCells[i])
    { (BCells+i)->ThreadNumber = MYTHREAD; }*/

    int ii=0;
    upc_forall(ii = 0; ii < B_CELLS_SIZE*THREADS; ii++; &BCells[ii]) {
        (BCells+ii)->ThreadNumber = MYTHREAD;
    }

    if(MYTHREAD==0)
    {
        fprintf(log_file,"\n:::: Initializing ::::\n");
        printf("\n:::: Initializing ::::\n");
    } // END OF THREAD ZERO

    ////////////////////////////////////////////////////
    ///////////////// INITIALIZE CELLS /////////////////
    ////////////////////////////////////////////////////

    upc_barrier;         // Synchronise

    upc_barrier;

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

    main_thread
        printf("Cell intialization completed\n");
    upc_barrier;

    /*
    print_cells_info(Cells);
    upc_barrier;
    print_boundary_type(Cells);
    upc_barrier;
    */

    //end_measure_time(tInitialization);

    upc_barrier;         // Synchronise
    if (MYTHREAD == 0)
    {
        printf("Size of Cells = %i\n", (int)sizeof(CellProps));
    }


    setCubeType();
    setThingsToGet();

    upc_barrier;
    // PUT INITIALIZED DATA TO BOUNDARIES
    //printTest("Before put",0);
    putCellsToShared();
    upc_barrier;
    //printTest("After put/Before get",0);
    getSharedToCells();
    upc_barrier;
    // COPY CELLS TO WCELLS (writing cells) TO WRITE DATA
    putCellsToWCells();
    upc_barrier;
    //write_cells_to_results(postproc_prog);


    tInitialization = tInitialization + (float)(clock()-tStart) / CLOCKS_PER_SEC;
    //free_mesh_data_matrices();  //BCNode

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
    upc_barrier;         // Synchronise

    main_while_loop(CollisionModel, CurvedBoundaries, OutletProfile, &Iterations, AutosaveAfter, AutosaveEvery,
                    postproc_prog, CalculateDragLift, ConvergenceCritVeloc, ConvergenceCritRho);

    tIterEnd = clock(); // End measuring time of main loop
    tIteration = (float)(tIterEnd - tIterStart ) / CLOCKS_PER_SEC;
    main_thread
        printf("Finished main loop\n");

    upc_barrier;         // Synchronise

    end_tests();
    export_data(postproc_prog, AutosaveEvery);

    upc_barrier;         // Synchronise
    free_vars();
    upc_barrier;         // Synchronise

}

void main_while_loop(int CollisionModel, int CurvedBoundaries, int OutletProfile, int *Iterations, int AutosaveAfter,
                     int AutosaveEvery, int postproc_prog, int CalculateDragLift, float ConvergenceCritVeloc,
                     float ConvergenceCritRho) {

    while ((Residuals[0] > ConvergenceCritVeloc || Residuals[1] > ConvergenceCritRho) && iter < (*Iterations))
    //while ((shared_total_Residuals[0] > ConvergenceCritVeloc || shared_total_Residuals[1] > ConvergenceCritRho) && iter < (*Iterations))
    {






//////////////// COLLISION ////////////////
        init_measure_time;
        CollisionStep(CollisionModel); ////////////////////// !!!!!!!!!!!!!!!!! CX CY!
        end_measure_time(tCollision);


        ////////////// UPDATE DISTR ///////////////
        init_measure_time;
        UpdateStep();
        end_measure_time(tUpdateF);
        upc_barrier;
        SAVE_ITERATION;


//////////////// PUT THREAD-BOUNDARY CELLS TO SHARED WITH CUBES ////////////////
        init_measure_time;
        putCellsToShared();
        end_measure_time(tBCells);
        upc_barrier;


//////////////// COPY SHARED BCELLS TO CELLS WITH CUBES////////////////
        init_measure_time;
        getSharedToCells();
        end_measure_time(tBCells);

////////////// STREAMING ///////////////
        init_measure_time;
        StreamingStep();
        end_measure_time(tStreaming);
        SAVE_ITERATION;

////////////// BOUNDARIES ///////////////
        init_measure_time;
        HandleBoundariesStep(OutletProfile, CurvedBoundaries);
        end_measure_time(tBoundaries);
        SAVE_ITERATION;

// UPDATE VELOCITY AND DENSITY
        init_measure_time;
        UpdateMacroscopicStep(CalculateDragLift);
        end_measure_time(tUpdateMacro);

////////////// Residuals ///////////////
        init_measure_time;
        ComputeResiduals(Cells, Residuals, sumVel0, sumVel1, sumRho0, sumRho1, CalculateDragLift, &iter, Iterations);
        end_measure_time(tResiduals);

        if(MYTHREAD==0)
            fprintf(resid_file,"%d %5.4e %5.4e %5.4e %f %f\n", iter, (iter+1.0)*(*Delta), Residuals[0], Residuals[1], Residuals[2], Residuals[3]);

        iter++; // update loop variable

////////////// Autosave ///////////////
        auto_save(AutosaveAfter, AutosaveEvery, postproc_prog);
        SAVE_ITERATION;
    }
    upc_barrier;
    //printf("Reached with T%i, iteration %d\n",MYTHREAD,iter);
    init_measure_time;
    putCellsToWCells();
    end_measure_time(tWriting);

    upc_barrier;

    //printf("Reached with T%i, after putCellsToWCells\n",MYTHREAD);
    upc_barrier;

//////////////////////////////////////////////////////
////////////// END OF MAIN WHILE CYCLE ///////////////
//////////////////////////////////////////////////////

}

//Step functions

void putCellsToShared(){
    FillLocalBCells();
    //printTest("Between Fill Local and upc_memput to shared",0);
//                     DESTINATION                           SOURCE                            SIZE
    UPUT( &BCells[MYTHREAD*B_CELLS_SIZE],&L_B_Cells[0],    B_CELLS_SIZE*sizeof(CellProps));
    //printf (" C[%i].F[1] = %f\n",lID(2,1,1),Cells[lID(2,1,1)].F[1]);
    //if (MYTHREAD == 0 && iter == 0)
    //printf (" SBC[%i].F[1] = %f\n",4,BCells[4].F[1]);
}
void FillLocalBCells() {
    int i,j,k;

    int c_BC = 0;         //Local BCells counter
    int min=1;
    int max=LAT;

    int f, e;
    int Ax, Dir, Pos;
    int a,b,c;

    //printTest("0", 0);
    //Faces

    for (f = 0; f<6; f++){
        Ax=getAx_F(f);
        Dir=getDir_F(f);

        a = (min+(max-min)*Dir);

        for (c=1; c<LAT+1;c++) {
            for (b=1; b<LAT+1;b++){
                i = eq(Ax,0)*a + eq(Ax,1)*c + eq(Ax,2)*b;
                j = eq(Ax,1)*a + eq(Ax,2)*c + eq(Ax,0)*b;
                k = eq(Ax,2)*a + eq(Ax,0)*c + eq(Ax,1)*b;

                L_B_Cells[c_BC] = Cells[getLocalID_LocalIndex(i, j, k)];

                c_BC++;
            }
        }
    }


    //Edges
    for (e = 0; e<12; e++){
        Ax=getAx_E(e);
        Pos=getPos_E(e);

        b = (min + (max-min)*(Pos%2));
        c = (min + (max-min)*(Pos/2));


        for (a=1; a<LAT+1;a++){

            i= eq(Ax,1)*c + eq(Ax,2)*b + eq(Ax,0)*a;
            j= eq(Ax,2)*c + eq(Ax,0)*b + eq(Ax,1)*a;
            k= eq(Ax,0)*c + eq(Ax,1)*b + eq(Ax,2)*a;
            L_B_Cells[c_BC] = Cells[getLocalID_LocalIndex(i, j, k)];
            c_BC++;
        }
    }
}


void setCubeType(){
    int n=MYTHREAD;


    //                                                              Corners
    //                                                                XYZ
    if (n == c000){                                                 //000
        cur_corner=0;
    }
    else if (n == c100){                                            //100
        cur_corner=1;
    }
    else if (n == c010){                                            //010
        cur_corner=2;
    }
    else if (n == c110){                                            //110
        cur_corner=3;
    }
    else if (n == c001){                                            //001
        cur_corner=4;
    }
    else if (n == c101){                                            //101
        cur_corner=5;
    }
    else if (n == c011){                                            //011
        cur_corner=6;
    }
    else if (n == c111){                                            //111
        cur_corner=7;
    }
        //                                                          Edges
        //                                                            XYZ
    else if (n < NTDX) {                                            //X00
        cur_edge=0;
    }
    else if (n-c010 >= 0 && (n-c010) < NTDX) {                      //X10
        cur_edge=1;
    }
    else if (n-c001 >= 0 && (n-c001) < NTDX) {                      //X01
        cur_edge=2;
    }
    else if (n-c011 >= 0 && (n-c011) < NTDX) {                      //X11
        cur_edge=3;
    }
    else if (n < NTDX*NTDY && n%NTDX==0) {                              //0Y0
        cur_edge=4;
    }
    else if (n-c001>=0 && (n-c001) < NTDX*NTDY && (n-c001)%NTDX==0) {   //0Y1
        cur_edge=5;
    }
    else if (n-c100>=0 && (n-c100) < NTDX*NTDY && (n-c100)%NTDX==0) {   //1Y0
        cur_edge=6;
    }
    else if (n-c101>=0 && (n-c101) < NTDX*NTDY && (n-c101)%NTDX==0) {   //1Y1
        cur_edge=7;
    }
        //
    else if (n%(NTDX*NTDY) == 0) {                                  //00Z
        cur_edge=8;
    }
    else if (n-c100 >= 0 && (n-c100)%(NTDX*NTDY) == 0) {            //10Z
        cur_edge=9;
    }
    else if (n-c010 >= 0 && (n-c010)%(NTDX*NTDY) == 0) {            //01Z
        cur_edge=10;
    }
    else if (n-c110 >= 0 && (n-c110)%(NTDX*NTDY) == 0) {            //11Z
        cur_edge=11;
    }
        //                                                          FACES
        //                                                            XYZ
    else if (n%NTDX == 0){                                               //0YZ
        cur_face=0;
    }
    else if (n-c100>=0 && (n-c100)%NTDX == 0){                      //1YZ
        cur_face=1;
    }
    else if ((n%(NTDX*NTDY)) < NTDX){                               //X0Z
        cur_face=2;
    }
    else if (n-c010>=0 && ((n-c010)%(NTDX*NTDY)) < NTDX){           //X1Z
        cur_face=3;
    }
    else if (n < (NTDX*NTDY)){                                      //XY0
        cur_face=4;
    }
    else if (n-c001 >= 0 && (n-c001) < (NTDX*NTDY)){                //XY1
        cur_face=5;
    }

        //                                                          INTERIOR
    else {

    }

    upc_barrier;

    /*printf("T=%2i ->  ",MYTHREAD);
    if(cur_corner != -1)
        printf("Corner %2i",cur_corner);
    else if(cur_edge != -1)
        printf("Edge   %2i",cur_edge);
    else if(cur_face != -1)
        printf("Face   %2i",cur_face);
    else
        printf("Interior");

    printf("\n");*/

    //PRINTING
    //printf("T%i: Finished cube type\n",MYTHREAD);
}  //Checked
void setThingsToGet(){

    if(cur_corner!=-1) {
        int X = getX_C(cur_corner);
        int Y = getY_C(cur_corner);
        int Z = getZ_C(cur_corner);

        corners_to_get[7-cur_corner] = 1;

        edges_to_get[getE(0,3-Y-2*Z)] = 1;
        edges_to_get[getE(1,3-Z-2*X)] = 1;
        edges_to_get[getE(2,3-X-2*Y)] = 1;

        faces_to_get[getF(0,1-X)] = 1;
        faces_to_get[getF(1,1-Y)] = 1;
        faces_to_get[getF(2,1-Z)] = 1;
    }//Checked

    else if(cur_edge!=-1) {

        int Ax = getAx_E(cur_edge);
        int Pos = getPos_E(cur_edge);

        for (int t=0;t<2;t++) {
            int X[3];
            for (int r=0;r<3;r++){
                X[r] = eq(Ax,r)*t +
                       eq(Ax,(r+1)%3)*(1-Pos/2) +
                       eq(Ax,(r+2)%3)*(1-Pos%2);
            }

            corners_to_get[getC(X[0],X[1],X[2])] = 1;
        }//Checked

        edges_to_get[getE(Ax,3-Pos)] = 1;
        edges_to_get[getE((Ax+1)%3 , 3-Pos/2)] = 1;
        edges_to_get[getE((Ax+1)%3 , 1-Pos/2)] = 1;
        edges_to_get[getE((Ax+2)%3 , 3-2*(Pos%2))] = 1;
        edges_to_get[getE((Ax+2)%3 , 2-2*(Pos%2))] = 1;
        //Checked

        for (int f = 0;f<6;f++) {
            if (f!=getF((Ax+1)%3 , Pos%2) && f!=getF((Ax+2)%3 , Pos/2)){
                faces_to_get[f] = 1;
            }
        }//Checked
    }

    else if(cur_face!=-1) {
        int Ax = getAx_F(cur_face);
        int Dir = getDir_F(cur_face);

        for (int b=0;b<2;b++){
            for(int a=0;a<2;a++) {
                int X[3];
                for(int r=0; r<3;r++){
                    X[r] = eq(Ax,r)*(1-Dir) +
                           eq(Ax,(r+2)%3)*a +
                           eq(Ax,(r+1)%3)*b;
                }

                corners_to_get[getC(X[0],X[1],X[2])] = 1;
            }
        }//Checked

        for (int e = 0; e<12; e++) {
            if(e!=getE((Ax+1)%3 , 2*Dir + 0) &&
               e!=getE((Ax+1)%3 , 2*Dir + 1) &&
               e!=getE((Ax+2)%3 , Dir + 0) &&
               e!=getE((Ax+2)%3 , Dir + 2)) {

                edges_to_get[e] = 1;
            }
        }//Checked

        for (int f = 0; f<6; f++) {
            if(f!=getF(Ax, Dir)){
                faces_to_get[f] = 1;
            }
        }
    }

    else {
        for (int c=0; c<8;c++) {
            corners_to_get[c] = 1;
        }
        for (int e=0; e<12;e++) {
            edges_to_get[e] = 1;
        }
        for (int f=0; f<6;f++) {
            faces_to_get[f] = 1;
        }
    }

    upc_barrier;


    //PRINTING
    //printf("T%i: Finished what to get\n",MYTHREAD);
    upc_barrier;
} //Checked
void getSharedToCells() {

    int X[3];
    int dX[3];
    getCubeCoords_CubeID(MYTHREAD,&X[0]);

    for (int face=0; face<6; face++) {
        if (faces_to_get[face] == 1){
            int Ax = getAx_F(face);
            int Dir = getDir_F(face);
            int op = getF(Ax,1-Dir);


            for (int r = 0; r<3; r++) {
                dX[r] = eq(Ax,r)*(2*Dir-1);
            }


            int pos_local = face*(LAT*LAT);
            int pos_shared = getCubeID(X[0] + dX[0],
                                       X[1] + dX[1],
                                       X[2] + dX[2])*B_CELLS_SIZE + op*(LAT*LAT);


            UGET(&L_B_Cells[pos_local],
                 &BCells[pos_shared],
                 (LAT*LAT)*sizeof(CellProps));


        }
    }
    PRINTING
        printf("T%i: Got faces\n",MYTHREAD);

    for (int edge = 0; edge <12; edge ++) {
        if(edges_to_get[edge] == 1) {
            int Ax = getAx_E(edge);
            int Pos = getPos_E(edge);
            int op = getE(Ax,3-Pos);

            for (int r = 0; r<3; r++) {
                dX[r] = eq(Ax , (r+2)%3)*(2*(Pos%2) - 1) +  eq(Ax , (r+1)%3)*(2*(Pos/2) - 1);
            }

            int pos_local = 6*LAT*LAT + edge*LAT;
            int pos_shared = getCubeID(X[0] + dX[0], X[1] + dX[1], X[2] + dX[2])
                             *B_CELLS_SIZE + 6*LAT*LAT + op*LAT;


            UGET(&L_B_Cells[pos_local],
                 &BCells[pos_shared],
                 LAT*sizeof(CellProps));
        }
    }
    PRINTING
        printf("T%i: Got edges\n",MYTHREAD);


    FillCellsWithLBCells();
}
void FillCellsWithLBCells() {
    int i,j,k;
    int c_BC = 0;         //Local BCells counter
    int min=0;
    int max=LAT+1;

    int f, e;
    int Ax, Dir, Pos;
    int a, b, c;

    for (f = 0; f<6; f++){
        Ax=getAx_F(f);
        Dir=getDir_F(f);

        a = (min+(max-min)*Dir);

        for (c=1; c<LAT+1;c++) {
            for (b=1; b<LAT+1;b++){

                //Faces
                i = eq(Ax,0)*a + eq(Ax,1)*c + eq(Ax,2)*b;
                j = eq(Ax,1)*a + eq(Ax,2)*c + eq(Ax,0)*b;
                k = eq(Ax,2)*a + eq(Ax,0)*c + eq(Ax,1)*b;
                Cells[getLocalID_LocalIndex(i, j, k)] = L_B_Cells[c_BC];
                c_BC++;
            }
        }
    }

    PRINTING
        printf("T%i: Copied Faces\n",MYTHREAD);

    //Edges
    for (e = 0; e<12; e++){
        Ax=getAx_E(e);
        Pos=getPos_E(e);

        b = (min + (max-min)*(Pos%2));
        c = (min + (max-min)*(Pos/2));


        for (a=1; a<LAT+1;a++){
            i= i=eq(Ax,1)*c + eq(Ax,2)*b + eq(Ax,0)*a;
            j= j=eq(Ax,2)*c + eq(Ax,0)*b + eq(Ax,1)*a;
            k= k=eq(Ax,0)*c + eq(Ax,1)*b + eq(Ax,2)*a;


            Cells[getLocalID_LocalIndex(i, j, k)] = L_B_Cells[c_BC];
            c_BC++;
        }
    }
    PRINTING
        printf("T%i: Copied Edges\n",MYTHREAD);
}

void printTest(const char * text,int it) {
    if (MYTHREAD == 0 && iter == it)
        printf("******%s\n", text);

    //FACE
    upc_barrier;
    if (MYTHREAD == 1 && iter == it) {
        printf("  Face (0,0,0)\n");
        int Cpos = getLocalID_LocalIndex(1, 1, 1);
        int LBpos = 0;
        int Shpos = MYTHREAD*B_CELLS_SIZE + LBpos;

        printf("    C_%i[%i]->ID = %i\n", MYTHREAD, Cpos, Cells[Cpos].ID);
        printf("    L_B_%i[%i]->ID = %i\n", MYTHREAD, LBpos, L_B_Cells[LBpos].ID);
        printf("    SH[%i]->ID = %i\n",Shpos, BCells[Shpos].ID);
    }
    upc_barrier;

    if (MYTHREAD == 0 && iter == it) {
        int Cpos = getLocalID_LocalIndex(3, 1, 1);
        int LBpos = 16;
        printf("    L_B_%i[%i]->ID = %i\n", MYTHREAD, LBpos, L_B_Cells[LBpos].ID);
        printf("    C_%i[%i]->ID = %i\n", MYTHREAD, Cpos, Cells[Cpos].ID);
    }

    upc_barrier;
}
void putCellsToWCells(){
    UPUT( &WCells[MYTHREAD * CELL_TOT_SIZE], &Cells[0], CELL_TOT_SIZE * sizeof(CellProps));
}
void FillLocalWCells(){
    int c_WC = 0;
    for (int k = 1; k < LAT+1; k++) {
        for (int j = 1; j < LAT + 1; j++) {
            for (int i = 1; i < LAT + 1; i++) {
                //L_W_Cells[c_WC] = Cells[lID(i,j,k)];
                c_WC++;
            }
        }
    }
}

void CollisionStep(int CollisionModel){


    switch(CollisionModel)
    {
        // BGKW
        case 1:
            for (int l_rID = 0; l_rID< BLOCKSIZE; l_rID++){

                BGKW(LocalID[l_rID], Omega);
                //}
                //}
            }

            // TRT
        case 2:
            for (int k = 1; k < LAT+1; k++) {
                for (int j = 1; j < LAT + 1; j++) {
                    for (int i = 1; i < LAT + 1; i++) {

                        TRT (Cells, i, w, cx, cy, opp, Omega, OmegaA);
                    }
                }
            }
            break;

            // MRT
        case 3:
            for (int k = 1; k < LAT+1; k++) {
                for (int j = 1; j < LAT + 1; j++) {
                    for (int i = 1; i < LAT + 1; i++) {

                        MRT(Cells, i, tm, stmiv);
                    }
                }
            }
            break;
        default:
            break;
    }
}

void UpdateStep(){
    for (int l_rID = 0; l_rID< BLOCKSIZE; l_rID++){
        UpdateF(Cells, LocalID[l_rID]);
        //}
        //}
    }
}

void StreamingStep(){
    for (int l_rID = 0; l_rID< BLOCKSIZE; l_rID++){
        int lID = LocalID[l_rID];
        for (int l = 0; l < 19; l++) {

            if (((Cells + lID)->StreamLattice[l]) == 1) {
                (Cells + lID)->F[l] = (Cells + lID + c[l])->METAF[l];
            }

        }
    }
}

void HandleBoundariesStep(int OutletProfile, int CurvedBoundaries){

    // INLET
    for (int l_rID = 0; l_rID< BLOCKSIZE; l_rID++){
        int lID = LocalID[l_rID];
        InletBC(Cells, lID);
        WallBC(Cells, lID, opp);
        EdgeBC(Cells, lID);
        CornerBC(Cells, lID);
    }

    // OUTLET
    int i = 0;
    int j = 0;
    switch(OutletProfile)
    {
        // set profile in outlet
        case 1:
            OutletBoundaries(Cells, j, i);
            break;

            // OPEN BOUNDARY
        case 2 :
            break;
        default: break;
    }
}
void UpdateMacroscopicStep(int CalculateDragLift){
    for (int l_rID = 0; l_rID< BLOCKSIZE; l_rID++){
        UpdateMacroscopic(Cells, LocalID[l_rID], CalculateDragLift);
    }
}





//Init/Alloc/free functions
void init_vars(int *postproc_prog) {
    iter = 0;                               // variables for loops
    iter_counter = 0; //
    sprintf(logFile,"Results/logFile.log"); // path of the .log file
    sprintf(testingFileName,"Results/testingThings.txt");
    AutosaveI = 1;                          // autosave i variable, will be incremented after every autosave
    ppp      = postproc_prog;       // for convenience ppp points to postproc_prog


    /*shared_total_Residuals[0] = 100;
    shared_total_Residuals[1] = 100;
    shared_total_Residuals[2] = 100;
    shared_total_Residuals[3] = 100;
    shared_total_Residuals[4] = 100;*/

    allocate_vars();

    time_meas_vars_init();
    init_cube_vars();
    init_localID();
    upc_barrier;
}
void init_cube_vars(){
    cur_corner=-1;
    cur_edge=-1;
    cur_face=-1;

    for (int i=0;i<6;i++) {
        faces_to_get[i] = 0;
    }
    for (int i=0;i<12;i++) {
        edges_to_get[i] = 0;
    }
    for (int i=0;i<8;i++) {
        corners_to_get[i] = 0;
    }
}
void time_meas_vars_init() {// Time measurement variables
    tInitialization  = 0.0; // Time measurement of Initialization
    tCellsInitialization  = 0.0; // Time measurement of Initialization
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
void init_localID(){
    for (int l_rID = 0; l_rID< BLOCKSIZE; l_rID++) {
        LocalID[l_rID] = getLocalID_LocRealID(l_rID);
    }
}

void alloc_cells() {//////////////////////////////////////////////////////
    // Allocate structure for the cell properties (see ShellFunctions.h)

    // New approach
    WCells = (shared_block(CELL_TOT_SIZE)    CellProps*) upc_all_alloc(THREADS, CELL_TOT_SIZE*sizeof(CellProps));
    BCells = (shared_block(B_CELLS_SIZE)     CellProps*) upc_all_alloc(THREADS, B_CELLS_SIZE*sizeof(CellProps));
    Cells = (CellProps*)calloc(CELL_TOT_SIZE,sizeof(CellProps));
    L_B_Cells = (CellProps*)calloc(B_CELLS_SIZE,sizeof(CellProps));
    if (!Cells) {
        printf("Cells mem failure, exiting \n");
        exit(EXIT_FAILURE);
    }
    if (!L_B_Cells) {
        printf("L_B_Cells mem failure, exiting \n");
        exit(EXIT_FAILURE);
    }
    //////////////////////////////////////////////////////

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
    NumNodes       = (shared int*)upc_alloc(1*sizeof(int));
    NumConn        = (shared int*)upc_alloc(1*sizeof(int));
    MaxInletCoordY = (shared double*)upc_alloc(1*sizeof(double));//--------------------------> DO I NEED THIS?
    MinInletCoordY = (shared double*)upc_alloc(1*sizeof(double));//--------------------------> DO I NEED THIS?
    NumInletNodes  = (shared int*)upc_alloc(1*sizeof(int));
}
void allocate_lattice_vars() {// D2Q9 Variables of the lattice
    w  = Create1DArrayDouble(19); // weight values for the directions

    c  = Create1DArrayInt(19);    //
    cx  = Create1DArrayInt(19);    // x coordinate of the discrete lattice directions
    cy  = Create1DArrayInt(19);    // y coordinate of the discrete lattice directions
    cz  = Create1DArrayInt(19);    // z coordinate of the discrete lattice directions
    opp = Create1DArrayInt(19);    // opposite vector
    norm = Create2DArrayInt(3,6);
    j_wall_unknown = Create2DArrayInt(5,6);

    LocalID = Create1DArrayInt(BLOCKSIZE);

}
void allocate_residuals() {// allocate residuals
    sumVel0   = Create1DArrayDouble(1);
    sumVel1   = Create1DArrayDouble(1);
    sumRho0   = Create1DArrayDouble(1);
    sumRho1   = Create1DArrayDouble(1);
    Residuals = Create1DArrayDouble(5);
}


void free_vars() {
    upc_barrier;
    UFREE(Delta);
    UFREE(MaxInletCoordY);
    UFREE(MinInletCoordY);
    UFREE(NumInletNodes);
    UFREE(NumNodes);
    UFREE(NumConn);

    upc_barrier;

    UAFREE(WCells);
    UAFREE(BCells);

    free(Cells);
    free(L_B_Cells);


    free(w);
    free(cx);
    free(cy);
    free(cz);
    free(c);
    free(opp);
}
void free_mesh_data_matrices() {

// We dont need these matrices anymore
    free(Nodes);
    free(BCconn);
}



//IO functions


void read_data(const char *NodeDataFile, const char *BCconnectorDataFile) {////////////////////////////////////////////////////
    ///////////////////// Read data ////////////////////
    ////////////////////////////////////////////////////

    Nodes   = ReadNodes(NodeDataFile);          // Read Node data
    BCconn  = ReadBCconn(BCconnectorDataFile);  // Read BCconn data
    CompDataNode(Nodes);
    //CompDataConn(BCconn);

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
            default: break;
        }
        switch(InletProfile)                      // 1:ON, 2:OFF
        {
            case 1: fprintf(log_file,">>> InletProfile      : ON\n" ); break;
            case 2: fprintf(log_file,">>> InletProfile      : OFF\n"); break;
            default: break;
        }
        switch(OutletProfile)                     // 1:ON, 2:OFF
        {
            case 1: fprintf(log_file,">>> OutletProfile     : ON\n" ); break;
            case 2: fprintf(log_file,">>> OutletProfile     : OFF\n"); break;
            default: break;
        }
        switch(CurvedBoundaries)                  // 1:ON, 2:OFF
        {
            case 1: fprintf(log_file,">>> CurvedBoundaries  : ON\n" ); break;
            case 2: fprintf(log_file,">>> CurvedBoundaries  : OFF\n"); break;
            default: break;
        }
        switch(postproc_prog)   // 1->Paraview (*.csv)     2->Tecplot
        {
            case 1: fprintf(log_file,">>> Results format    : Paraview (*.csv)\n" ); break;
            case 2: fprintf(log_file,">>> Results format    : Tecplot (*.dat)\n"); break;
            default: break;
        }
        if (CalculateDragLift != 0)
            fprintf(log_file,">>> Drag, lift @ BC   : %d\n", CalculateDragLift);
        else
            fprintf(log_file,">>> Drag, lift was not calculated\n");

        fprintf(log_file,"\n:::: Calculated variables from mesh :::: \n");
        fprintf(log_file,">>> Grid spacing        = %f\n", *Delta);
        fprintf(log_file,">>> # of nodes in x (n) = %d\n", NN);
        fprintf(log_file,">>> # of nodes in y (m) = %d\n", NM);
        fprintf(log_file,">>> # of nodes in z (l) = %d\n", NL);
        //fprintf(log_file,">>> NumInletNodes       = %d\n", *NumInletNodes);
//fprintf(log_file,">>> MaxInletCoordY      = %f\n", *MaxInletCoordY);
//fprintf(log_file,">>> MinInletCoordY      = %f\n", *MinInletCoordY);

        fprintf(log_file,"\n:::: Parallel properties :::: \n");
        fprintf(log_file,">>> # of threads        = %d\n", THREADS);
        fprintf(log_file,">>> BlockSize           = %d\n", BLOCKSIZE);

// In case of no autosave
        sprintf(AutosaveOutputFile, "NOWHERE!");

    } // END OF THREAD ZERO

}


void auto_save(int AutosaveAfter, int AutosaveEvery, int postproc_prog) {
    if(iter == (AutosaveEvery * AutosaveI))
    {
        AutosaveI++;
        if(iter>AutosaveAfter) {
            init_measure_time;
            switch (postproc_prog) {
                case 1:
                    sprintf(AutosaveOutputFile, "Results/autosave_iter/autosave_iter%05d.csv", iter);
                    break;
                case 2:
                    sprintf(AutosaveOutputFile, "Results/autosave_iter/autosave_iter%05d.dat", iter);
                    break;
                default:
                    break;
            }
            putCellsToWCells(); // Put information to WCells and write (Write Cells)
            upc_barrier;
            if (MYTHREAD == 0){ // AUTOSAVE
                WriteResults(AutosaveOutputFile, ppp);
                char AutosaveTime[200];
                sprintf(AutosaveTime, "Results/autosave_times/ParallelTimeMeasuerment_iter%05d.dat", iter);
                print_times(AutosaveTime, AutosaveEvery);
                printf("Autosaved n%i Iter n%i\n", AutosaveI - 1, iter);
            }

            end_measure_time(tWriting);
        }

    }
}

void save_iteration(int postproc_prog, int AutosaveEvery) {

    //if ((iter * 20)%AutosaveEvery == 0) {
    //if (iter < 10 || (iter%5 == 0 && iter < 100)) {
    //printTest("Before saving",0);
    if (iter==0) {

        //if(MYTHREAD == 0)
        //printf("Saving iteration %i, Rho residual = %lf\n",iter,Residuals[1]);
        //printTest("Before macro",0);

        UpdateMacroscopicStep(0);
        init_measure_time;
        switch (postproc_prog) {
            case 1:
                sprintf(IterationOutputFile, "Results/iterations/iter.csv.%i", iter_counter);
                break;
            case 2:
                sprintf(IterationOutputFile, "Results/iterations/iter.dat.%i", iter_counter);
                break;
            default:
                break;
        }
        //printTest("Before WCells",0);
        putCellsToWCells();
        //printTest("Before Writing",0);
        if (MYTHREAD == 0) // AUTOSAVE
            WriteResults(IterationOutputFile, ppp);

        upc_barrier;
        end_measure_time(tWriting);


        iter_counter++;
    }
    //printTest("After saving",0);
}

void write_boundary_cells_to_results(int postproc_prog) {

// Write boundary cells to Results to see how mesh was distributed
    if(MYTHREAD==0)
    {
        switch(postproc_prog)  {
            case 1: sprintf(fnMemCopyRes, "Results/MyBoundaryCells.csv");   break;
            case 2: sprintf(fnMemCopyRes, "Results/MyBoundaryCells.dat");   break;
            default: break;
        }
        printf("%s\n",fnMemCopyRes);
        WriteBCells(fnMemCopyRes, ppp);
    }
}

void save_init_data(int postproc_prog) {// Write Initialized data
    switch(postproc_prog){
        case 1: sprintf(OutputFile, "Results/InitialData.csv");
            break;
        case 2: sprintf(OutputFile, "Results/InitialData.dat");
            break;
        default:
            break;
    }
    putCellsToWCells();
    WriteResults(OutputFile, ppp);
    printf("\nInitialized data was written to %s\n", OutputFile);
}
void export_data(int postproc_prog, int AutosaveEvery) {
    if(MYTHREAD == 0) // EXPORT DATA, TIME MEASUREMENT RESULTS
    {
        // Close residuals file
        fclose(resid_file);

        clock_t tEnd = clock();
        float tOverall = (float)(tEnd - tStart) / CLOCKS_PER_SEC; // Calculate elapsed time

        fprintf(log_file,"\nOverall calculations took %f seconds\n", tOverall);
        fprintf(log_file,"Main while loop took %f seconds\n",        tIteration);
        fprintf(log_file,"Cells init took %f seconds\n",             tCellsInitialization);
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

        print_times("Results/ParallelTimeMeasuerment.dat",AutosaveEvery);

        // Write final data
        switch(postproc_prog){
            case 1: sprintf(FinalOutputFile, "Results/FinalData.csv"); break;
            case 2: sprintf(FinalOutputFile, "Results/FinalData.dat"); break;
            default: break;
        }

        //putCellsToWCells();
        WriteResults(FinalOutputFile,  ppp);

        // Write information for user
        printf("\n\nLog was written to %s\n", logFile);
        printf("Last autosave result can be found at %s\n", AutosaveOutputFile);
        printf("Residuals were written to Results/residuals.dat\n");
        printf("Profiling results were written to Results/ParallelTimeMeasuerment.dat\n");
        printf("Final results were written to %s\n", FinalOutputFile);
        //write_boundary_cells_to_results(*ppp);

        //WriteBCells(fnMemCopyRes, ppp);
        //printf("BCells were written!\n");
        //puts("BCells were written!");
    } // END OF THREAD ZERO


}
void print_times(const char* fname, int AutosaveEvery) {



    clock_t tEnd = clock();
    float tOverall = (float)(tEnd - tStart) / CLOCKS_PER_SEC; // Calculate elapsed time

    TimeMeasurementFile = fopen(fname,"w");

    fprintf(TimeMeasurementFile,"tOverall %f\n",        tOverall);
    fprintf(TimeMeasurementFile,"tIteration %f\n",      tIteration);
    fprintf(TimeMeasurementFile,"tCellsInitialization %f\n", tCellsInitialization);
    fprintf(TimeMeasurementFile,"tInitialization %f\n", tInitialization);
    fprintf(TimeMeasurementFile,"tCollision %f\n",      tCollision);
    fprintf(TimeMeasurementFile,"tUpdateF %f\n",        tUpdateF);
    fprintf(TimeMeasurementFile,"tStreaming %f\n",      tStreaming);
    fprintf(TimeMeasurementFile,"tBoundaries %f\n",     tBoundaries);
    fprintf(TimeMeasurementFile,"tUpdateMacro %f\n",    tUpdateMacro);
    fprintf(TimeMeasurementFile,"tResiduals %f\n",      tResiduals);
    fprintf(TimeMeasurementFile,"tBCells %f\n",         tBCells);
    fprintf(TimeMeasurementFile,"tWriting %f\n\n",    tWriting);


    fprintf(TimeMeasurementFile,"n_iterations %i\n",    iter);
    fprintf(TimeMeasurementFile,"tIteration %f\n",      tIteration/iter);
    fprintf(TimeMeasurementFile,"tCollision %f\n",      tCollision/iter);
    fprintf(TimeMeasurementFile,"tUpdateF %f\n",        tUpdateF/iter);
    fprintf(TimeMeasurementFile,"tStreaming %f\n",      tStreaming/iter);
    fprintf(TimeMeasurementFile,"tBoundaries %f\n",     tBoundaries/iter);
    fprintf(TimeMeasurementFile,"tUpdateMacro %f\n",    tUpdateMacro/iter);
    fprintf(TimeMeasurementFile,"tResiduals %f\n",      tResiduals/iter);
    fprintf(TimeMeasurementFile,"tBCells %f\n\n",       tBCells/iter);

    fprintf(TimeMeasurementFile,"tWriting %f\n",        tWriting/AutosaveEvery);

    fclose(TimeMeasurementFile);

}




void print_cells_info(CellProps* Cells) {
    /*upc_barrier;
    FILE* out_cells_file;
    char out_cells_filename [50];
    sprintf(out_cells_filename,"Results/outCells/T_%i.dat",MYTHREAD);
    out_cells_file = fopen(out_cells_filename,"w");
    for (int t_to_print_from = 0; t_to_print_from<THREADS;t_to_print_from++){
        if (MYTHREAD == t_to_print_from)
        {
            //printf("Going to Thread %i\n",MYTHREAD);
            fprintf(out_cells_file,"Printing from thread %i\n",MYTHREAD);
            fprintf(out_cells_file,"   ID   |  i |  j |  k ||    x    |    y    |    z   | s14 |\n");
            for(int cell_to_print = LAYER; cell_to_print< BLOCKSIZE+LAYER; cell_to_print++) {
                print_cell_line(out_cells_file,Cells+cell_to_print);
            }
        }
        upc_barrier;
    }
    fclose(out_cells_file);
    upc_barrier;*/
    //printf("Cells info printed\n");
}
void print_boundary_type(CellProps* Cells) {
    /*main_thread
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
    for (int node_to_look = LAYER; node_to_look<BLOCKSIZE+LAYER; node_to_look++) {
        int BT;
        //printf("TEST, thread %i\n",MYTHREAD);
        if (node_to_look == LAYER)    {
            //printf("Thread: %i,Node: %i,BT: i\n",MYTHREAD,node_to_look,(Cells+node_to_look)->Boundary);
            int index_n, index_i, index_j, index_k;
            index_n = (Cells+LAYER)->ID;
            index_k = index_n/(NM*NN);
            index_j = (index_n - index_k * NM* NN)/NM;
            index_i = (index_n - index_k * NM* NN - index_j * NM);
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

    const char* b_filename[4];
    const char* b1_filename = "Results/boundary/solidplane_boundary.dat";
    const char* b2_filename = "Results/boundary/fluidplane_boundary.dat";
    const char* b3_filename = "Results/boundary/edge_boundary.dat";
    const char* b4_filename = "Results/boundary/corner_boundary.dat";
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
        printf("\n");}*/

    /*if (MYTHREAD == 0){
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
                fprintf(b_file,"   ID   |  i |  j |  k ||    x    |    y    |    z   | s14 |\n");
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
    free(N_B[3]);*/

}
void print_cell_line(FILE* file, const CellProps* Cell) {
    int index_n, index_i, index_j, index_k;
    index_n = Cell->ID;
    index_k = index_n/(NM*NN);
    index_j = (index_n - index_k * NM* NN)/NM;
    index_i = (index_n - index_k * NM* NN - index_j * NM);
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


//Other functions

void calc_collision_freq(float Viscosity) {// Calculate collision frequency
    Omega  = 1.0/(3.*Viscosity+0.5);
    OmegaA = 8*(2-Omega)/(8-Omega);
}
