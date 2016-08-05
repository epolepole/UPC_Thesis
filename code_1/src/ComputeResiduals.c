/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ComputeResiduals.h"
#include "tests.h"


void ComputeResiduals(CellProps *Cells, double* Residuals,
                      double* sumVel0, double* sumVel1,
                      double* sumRho0, double* sumRho1,
                      int ComputeDragLift, int* iter, int* Iterations)
{

    // Update variables
    *sumVel0=*sumVel1;
    *sumRho0=*sumRho1;

    // Create variables for residuals
    //double ResVel   = 0.0;
    //double ResRho   = 0.0;
    double ResDragX  = 0.0;
    double ResDragY  = 0.0;
    double ResLift  = 0.0;
    //double LsumVel1 = 0.0;
    //double LsumRho1 = 0.0;
    double PUTtmp[5];
    double GETtmp[5*THREADS];

    double L2n   = 0.0;  // L2 norm
    double L2n_w = 0.0;  // L2 norm weighted

    // Loop variables
    int i, j, k;

    *sumVel1 = 0;
    //*sumRho1 = 0;

    // sum up velocity and density
    for(i = LAYER;  i < LAYER+BLOCKSIZE;  i++)
    {

        //LsumVel1 = LsumVel1 + sqrt( ((Cells+i)->U)*((Cells+i)->U)  +  ((Cells+i)->V)*((Cells+i)->V)  );
        //LsumRho1 = LsumRho1 + (Cells+i)->Rho;

        for (k=0; k<19;k++)
        { L2n = L2n + pow(((Cells+i)->F[k]-(Cells+i)->METAF[k]),2); }


        //For now checking on every external surface, if ComputeDragLift == 1
        if ((Cells+i)->BoundaryID == ComputeDragLift)  // BE AWARE!! check in the STAR-CD files the ID of the surface where you want to check the drag/lift.
        {
            ResDragX += (Cells+i)->DragXF;
            ResDragY += (Cells+i)->DragYF;
            ResLift += (Cells+i)->LiftF;
        }
    }

    //The same?
    PUTtmp[0] = L2n;
    PUTtmp[1] = L2n;

    PUTtmp[2] = ResDragX;
    PUTtmp[3] = ResDragY;
    PUTtmp[4] = ResLift;

    *sumVel1 = 0;
    //*sumRho1 = 0;
    ResDragX  = 0;
    ResDragY  = 0;
    ResLift  = 0;

    UPUT( &sResiduals[0+5*MYTHREAD] , &PUTtmp[0],  5*sizeof(double) );
//  tInstant1 = clock(); // Start measuring time
    upc_barrier;
//  tInstant2 = clock(); // End of time measuremet
    UGET( &GETtmp[0] , &sResiduals[0] , (THREADS*5)*sizeof(double) );
//  *tSendRcv =  (tInstant2-tInstant1) ;

    //printf("iter %d TH%d Res Synced\n", iter, MYTHREAD);

    if(MYTHREAD==0)
    {
        //printf("Th%d residualcalc\n",MYTHREAD);
        for (i = 0; i < THREADS; i++)
        {

            *sumVel1 = *sumVel1 + GETtmp[0+5*i];
            //*sumRho1 = *sumRho1 + GETtmp[1+4*i];
            ResDragX  = ResDragX  + GETtmp[2+5*i];
            ResDragY  = ResDragY  + GETtmp[3+5*i];
            ResLift  = ResLift  + GETtmp[4+5*i];

        }

        L2n = *sumVel1;
        // Calculate residuals
        L2n_w = sqrt(L2n/(NN*NM*NL));
        L2n = sqrt(L2n);

        Residuals[0] = L2n;
        Residuals[1] = L2n_w;
        Residuals[2] = ResDragX;
        Residuals[3] = ResDragY;
        Residuals[4] = ResLift;

    }

    upc_barrier;
    if(L2n!=L2n) // if density residuals are NaN
    {
        printf("\nDIVERGENCE!\n");
        test_all(Cells,iter);
        upc_barrier;
        exit(1); // ERROR!
        *iter  = *Iterations+1;
        Residuals[0] = 1;
        Residuals[1] = 1;
        Residuals[2] = 0;
        Residuals[3] = 0;

        //upc_global_exit(1);
        //exit(1); // ERROR!
        end_tests();
    }
}*/



#include <stdio.h>
#include <math.h>
#include <CellFunctions.h>
//#include <Iterate.h>

#include "ShellFunctions.h"

void ComputeResiduals(CellProps *Cells, double* Residuals, double* sumVel0, double* sumVel1, double* sumRho0, double* sumRho1, int ComputeDragLift, int* iter, int* Iterations)
{
    // Update variables
    *sumVel0=*sumVel1;
    *sumRho0=*sumRho1;

    // Create variables for residuals
    double ResVel  = 0.0;
    double ResRho  = 0.0;
    double ResDrag = 0.0;
    double ResLift = 0.0;

    // Loop variables

    *sumVel1 = 0;
    *sumRho1 = 0;
    // sum up velocity and density

    for (int l_rID = 0; l_rID< BLOCKSIZE_NEW; l_rID++){
                int lID = LocalID[l_rID];
                *sumVel1 = *sumVel1 + sqrt(pow((Cells + lID)->U, 2) + pow((Cells + lID)->V, 2) + pow((Cells + lID)->W, 2));
                *sumRho1 = *sumRho1 + (Cells + lID)->Rho;
                /*if ((Cells+i)->Rho < 0) {
                    printf("Iteration %i\n",*iter);
                    printf("T=%i, rho[%i] = %f\n",MYTHREAD,i+LAYER+MYTHREAD*BLOCKSIZE,(Cells+i)->Rho);
                    *iter  = *Iterations+1;
                    //upc_global_exit(1);
                    //exit(1); // ERROR!
                }*/

                /*if ((Cells+i)->BoundaryID == ComputeDragLift)  // BE AWARE!! check in the STAR-CD files the ID of the surface where you want to check the drag/lift.
                {
                  ResDrag += (Cells+i)->DragF;
                  ResLift += (Cells+i)->LiftF;
                }*/

    }

    sResiduals[0+5*MYTHREAD] = *sumVel1;
    sResiduals[1+5*MYTHREAD] = *sumRho1;
    sResiduals[2+5*MYTHREAD] = ResDrag;
    sResiduals[3+5*MYTHREAD] = ResDrag;
    sResiduals[4+5*MYTHREAD] = ResLift;

    *sumVel1 = 0;
    *sumRho1 = 0;
    ResDrag  = 0;
    ResLift  = 0;


    //printf("iter %d TH%d Res Synced\n", *iter, MYTHREAD);

    upc_barrier;

    if(MYTHREAD==0) {
        //printf("Th%d residualcalc\n",MYTHREAD);
        for (int i = 0; i < THREADS; i++)
        {
            *sumVel1 = *sumVel1 + sResiduals[0+5*i];
            *sumRho1 = *sumRho1 + sResiduals[1+5*i];
            ResDrag  = ResDrag  + sResiduals[2+5*i];
            ResLift  = ResLift  + sResiduals[3+5*i];
        }

        // Calculate residuals
        ResVel = sqrt( pow( ((*sumVel0-*sumVel1)/max(*sumVel0,*sumVel1)) ,2) );
        ResRho = sqrt( pow( ((*sumRho0-*sumRho1)/max(*sumRho0,*sumRho1)) ,2) );

        Residuals[0] = ResVel;
        Residuals[1] = ResRho;
        Residuals[2] = ResDrag;
        Residuals[3] = ResLift;

    }

    upc_barrier;

    if(MYTHREAD == 0) {
        shared_total_Residuals[0] = ResVel;
        shared_total_Residuals[1] = ResRho;
        shared_total_Residuals[2] = ResDrag;
        shared_total_Residuals[2] = ResDrag;
        shared_total_Residuals[3] = ResLift;
    }
    /*if ((MYTHREAD == 0 && *iter%50 == 0) || ResRho != ResRho){
        printf("SumRho0 = %f\n",*sumRho0);
        printf("SumRho1 = %f\n",*sumRho1);
        printf("ResRho = %f\n",ResRho);
    }*/

    /*if (*iter==367){
        for (int t=0;t<THREADS;t++) {
            if (t==MYTHREAD) {
                printf("T%i: Vr=%8.8f , Rr=%8.8f\n", t, Residuals[0], Residuals[1]);
                printf("T%i: sVr=%8.8f , sRr=%8.8f\n", t, shared_total_Residuals[0], shared_total_Residuals[1]);
            }
            for (int i=0;i<10000;i++)
            upc_barrier;
        }
    }*/
    upc_barrier;

    if(ResRho!=ResRho) // if density residuals are NaN
    {
        if(MYTHREAD==0){
            printf("ResRho = %f\n", ResRho);
            printf("\nDensity divergence!\n");
        }
        *iter  = *Iterations+1;
        Residuals[0] = 1;
        Residuals[1] = 1;
        Residuals[2] = 0;
        Residuals[3] = 0;

        //upc_global_exit(1);
        //exit(1); // ERROR!
    }
}

