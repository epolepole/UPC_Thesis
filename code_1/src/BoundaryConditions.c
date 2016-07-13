#include <stdio.h>                  // printf();
#include <stdlib.h>                 // for malloc();
#include "CellFunctions.h"
#include <math.h>

#include "BoundaryConditions.h" // convenience


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
////////////////////// Inlet boundary treatment /////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void InletBC(CellProps *Cells, int i)
{
    MyReal Rho = 0.0;
    MyReal U, V, W, N1, N2, XMAX, YMAX, ZMAX;
    XMAX = 1;
    YMAX = 1;
    ZMAX = 1;

    U = (Cells + i)->Uo;
    V = (Cells + i)->Vo;
    W = (Cells + i)->Wo;

    if ((Cells + i)->Boundary == 2) // INLET BOUNDARY
    {

        if ((Cells + i)->CoordZ == 0) // INLET PLANE xy, z=0
        {
            //Unknowns: Rho, f5, f11, f12, f15 and f16

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[1] + (Cells + i)->F[2] + (Cells + i)->F[3] + (Cells + i)->F[4] +
                   (Cells + i)->F[7] + (Cells + i)->F[8] + (Cells + i)->F[9] + (Cells + i)->F[10] +
                   2 * ((Cells + i)->F[6] + (Cells + i)->F[13] + (Cells + i)->F[14] + (Cells + i)->F[17] +
                        (Cells + i)->F[18])) / (1 - W);

            N1 = 0.5 * ((Cells + i)->F[1] + (Cells + i)->F[7] + (Cells + i)->F[9] - (Cells + i)->F[2] - (Cells + i)->F[8] -
                        (Cells + i)->F[10]) - Rho * U / 3;
            N2 = 0.5 * ((Cells + i)->F[3] + (Cells + i)->F[7] + (Cells + i)->F[8] - (Cells + i)->F[4] - (Cells + i)->F[9] -
                        (Cells + i)->F[10]) - Rho * V / 3;

            (Cells + i)->F[5]  = (Cells + i)->F[6] + (Rho * W) / 3;
            (Cells + i)->F[11] = (Cells + i)->F[14] + Rho / 6 * (W + U) - N1;
            (Cells + i)->F[12] = (Cells + i)->F[13] + Rho / 6 * (W - U) + N1;
            (Cells + i)->F[15] = (Cells + i)->F[18] + Rho / 6 * (W + V) - N2;
            (Cells + i)->F[16] = (Cells + i)->F[17] + Rho / 6 * (W - V) + N2;
        }
        else if ((Cells + i)->CoordZ == ZMAX) // INLET PLANE xy, z=zmax
        {
            //Unknowns: Rho, f6, f13, f14, f17 and f18

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[1] + (Cells + i)->F[2] + (Cells + i)->F[3] + (Cells + i)->F[4] +
                   (Cells + i)->F[7] + (Cells + i)->F[8] + (Cells + i)->F[9] + (Cells + i)->F[10] +
                   2 * ((Cells + i)->F[5] + (Cells + i)->F[11] + (Cells + i)->F[12] + (Cells + i)->F[15] +
                        (Cells + i)->F[16])) / (1 + W);

            N1 = 0.5 * ((Cells + i)->F[1] + (Cells + i)->F[7] + (Cells + i)->F[9] - (Cells + i)->F[2] - (Cells + i)->F[8] -
                        (Cells + i)->F[10]) - Rho * U / 3;
            N2 = 0.5 * ((Cells + i)->F[3] + (Cells + i)->F[7] + (Cells + i)->F[8] - (Cells + i)->F[4] - (Cells + i)->F[9] -
                        (Cells + i)->F[10]) - Rho * V / 3;

            (Cells + i)->F[6]  = (Cells + i)->F[5] - (Rho * W) / 3;
            (Cells + i)->F[13] = (Cells + i)->F[12] + Rho / 6 * (-W + U) - N1;
            (Cells + i)->F[14] = (Cells + i)->F[11] + Rho / 6 * (-W - U) + N1;
            (Cells + i)->F[17] = (Cells + i)->F[16] + Rho / 6 * (-W + V) - N2;
            (Cells + i)->F[18] = (Cells + i)->F[15] + Rho / 6 * (-W - V) + N2;
            //printf("Always ZMAX\n");
        }

        else if ((Cells + i)->CoordX == 0) // INLET PLANE yz, x=0
        {
            //Unknowns: Rho, f1, f7, f9, f11, and f13

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[3] + (Cells + i)->F[4] + (Cells + i)->F[5] + (Cells + i)->F[6] +
                   (Cells + i)->F[15] + (Cells + i)->F[16] + (Cells + i)->F[17] + (Cells + i)->F[18] +
                   2 *
                   ((Cells + i)->F[2] + (Cells + i)->F[8] + (Cells + i)->F[10] + (Cells + i)->F[12] + (Cells + i)->F[14])) /
                  (1 + U);

            N1 = 0.5 *
                 ((Cells + i)->F[3] + (Cells + i)->F[15] + (Cells + i)->F[17] - (Cells + i)->F[4] - (Cells + i)->F[16] -
                  (Cells + i)->F[18]) - Rho * V / 3;
            N2 = 0.5 *
                 ((Cells + i)->F[5] + (Cells + i)->F[15] + (Cells + i)->F[16] - (Cells + i)->F[6] - (Cells + i)->F[17] -
                  (Cells + i)->F[18]) - Rho * W / 3;

            (Cells + i)->F[1]  = (Cells + i)->F[2]  + (Rho * U) / 3;
            (Cells + i)->F[7]  = (Cells + i)->F[10] + Rho / 6 * (U + V) - N1;
            (Cells + i)->F[9]  = (Cells + i)->F[8]  + Rho / 6 * (U - V) + N1;
            (Cells + i)->F[11] = (Cells + i)->F[14] + Rho / 6 * (U + W) - N2;
            (Cells + i)->F[13] = (Cells + i)->F[12] + Rho / 6 * (U - W) + N2;
        }

        else if ((Cells + i)->CoordX == XMAX) // INLET PLANE yz, x=xmax
        {
            //Unknowns: Rho, f2, f8, f10, f12 and f14

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[3] + (Cells + i)->F[4] + (Cells + i)->F[5] + (Cells + i)->F[6] +
                   (Cells + i)->F[15] + (Cells + i)->F[16] + (Cells + i)->F[17] + (Cells + i)->F[18] +
                   2 *
                   ((Cells + i)->F[1] + (Cells + i)->F[7] + (Cells + i)->F[9] + (Cells + i)->F[11] + (Cells + i)->F[13])) /
                  (1 - U);

            N1 = 0.5 *
                 ((Cells + i)->F[3] + (Cells + i)->F[15] + (Cells + i)->F[17] - (Cells + i)->F[4] - (Cells + i)->F[16] -
                  (Cells + i)->F[18]) - Rho * V / 3;
                  N2 = 0.5 *
                       ((Cells + i)->F[5] + (Cells + i)->F[15] + (Cells + i)->F[16] - (Cells + i)->F[6] - (Cells + i)->F[17] -
                        (Cells + i)->F[18]) - Rho * W / 3;

            (Cells + i)->F[2]  = (Cells + i)->F[1]  - (Rho * U) / 3;
            (Cells + i)->F[8]  = (Cells + i)->F[9]  + Rho / 6 * (-U + V) - N1;
            (Cells + i)->F[10] = (Cells + i)->F[7]  + Rho / 6 * (-U - V) + N1;
            (Cells + i)->F[12] = (Cells + i)->F[13] + Rho / 6 * (-U + W) - N2;
            (Cells + i)->F[14] = (Cells + i)->F[11] + Rho / 6 * (-U - W) + N2;
        }

        else if ((Cells + i)->CoordY == 0) // INLET PLANE xz, y=0
        {
            //Unknowns: Rho, f3, f7, f8, f15 and f17

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[1] + (Cells + i)->F[2] + (Cells + i)->F[5] + (Cells + i)->F[6] +
                   (Cells + i)->F[11] + (Cells + i)->F[12] + (Cells + i)->F[13] + (Cells + i)->F[14] +
                   2 *
                   ((Cells + i)->F[4] + (Cells + i)->F[9] + (Cells + i)->F[10] + (Cells + i)->F[16] + (Cells + i)->F[18])) /
                  (1 - V);

            N1 = 0.5 *
                 ((Cells + i)->F[1] + (Cells + i)->F[11] + (Cells + i)->F[13] - (Cells + i)->F[2] - (Cells + i)->F[12] -
                  (Cells + i)->F[14]) - Rho * U / 3;
            N2 = 0.5 *
                 ((Cells + i)->F[5] + (Cells + i)->F[11] + (Cells + i)->F[12] - (Cells + i)->F[6] - (Cells + i)->F[13] -
                  (Cells + i)->F[14]) - Rho * W / 3;

            (Cells + i)->F[3]  = (Cells + i)->F[4]  + (Rho * V) / 3;
            (Cells + i)->F[7]  = (Cells + i)->F[10] + Rho / 6 * (V + U) - N1;
            (Cells + i)->F[8]  = (Cells + i)->F[9]  + Rho / 6 * (V - U) + N1;
            (Cells + i)->F[15] = (Cells + i)->F[18] + Rho / 6 * (V + W) - N2;
            (Cells + i)->F[17] = (Cells + i)->F[16] + Rho / 6 * (V - W) + N2;
        }

        else if ((Cells + i)->CoordY == YMAX) // INLET PLANE xz, y=ymax
        {
            //Unknowns: Rho, f4, f9, f10, f16 and f18

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[1] + (Cells + i)->F[2] + (Cells + i)->F[5] + (Cells + i)->F[6] +
                   (Cells + i)->F[11] + (Cells + i)->F[12] + (Cells + i)->F[13] + (Cells + i)->F[14] +
                   2 *
                   ((Cells + i)->F[3] + (Cells + i)->F[7] + (Cells + i)->F[8] + (Cells + i)->F[15] + (Cells + i)->F[17])) /
                  (1 + V);

            N1 = 0.5 *
                 ((Cells + i)->F[1] + (Cells + i)->F[11] + (Cells + i)->F[13] - (Cells + i)->F[2] - (Cells + i)->F[12] -
                  (Cells + i)->F[14]) - Rho * U / 3;
            N2 = 0.5 *
                 ((Cells + i)->F[5] + (Cells + i)->F[11] + (Cells + i)->F[12] - (Cells + i)->F[6] - (Cells + i)->F[13] -
                  (Cells + i)->F[14]) - Rho * W / 3;

            (Cells + i)->F[4]  = (Cells + i)->F[3]  - (Rho * V) / 3;
            (Cells + i)->F[9]  = (Cells + i)->F[8]  + Rho / 6 * (-V + U) - N1;
            (Cells + i)->F[10] = (Cells + i)->F[7]  + Rho / 6 * (-V - U) + N1;
            (Cells + i)->F[16] = (Cells + i)->F[17] + Rho / 6 * (-V + W) - N2;
            (Cells + i)->F[18] = (Cells + i)->F[15] + Rho / 6 * (-V - W) + N2;
        }
    }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
////////////////////// Outlet boundary treatment ////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void OutletBoundaries(CellProps *Cells, int j, int i)
{
    double RhoE;
    /*if ((Cells+j*(*n)+i)->BC_ID[1]==3) // outlet boundary on the right side of the domain
    {

        RhoE = ((Cells+j*(*n)+i)->F[0]+(Cells+j*(*n)+i)->F[2]+(Cells+j*(*n)+i)->F[4]
                +2.0*((Cells+j*(*n)+i)->F[1]+(Cells+j*(*n)+i)->F[5]+(Cells+j*(*n)+i)->F[8]))/(1-(Cells+j*(*n)+i)->Uo);

        (Cells+j*(*n)+i)->F[3] = (Cells+j*(*n)+i)->F[1]-2*RhoE*((Cells+j*(*n)+i)->Uo)/3.0;

        (Cells+j*(*n)+i)->F[7] = (Cells+j*(*n)+i)->F[5]-RhoE*((Cells+j*(*n)+i)->Uo)/6.0;

        (Cells+j*(*n)+i)->F[6] = (Cells+j*(*n)+i)->F[8]-RhoE*((Cells+j*(*n)+i)->Uo)/6.0;
    }
    if ((Cells+j*(*n)+i)->BC_ID[2]==3)
    {
        // FILL!!
    }
    if ((Cells+j*(*n)+i)->BC_ID[3]==3)
    {
        // FILL!!
    }
    if ((Cells+j*(*n)+i)->BC_ID[4]==3)
    {
        // FILL!!
    }*/

}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////// Wall boundary treatment /////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void WallBC(CellProps *Cells, int i, int* opp)
{

    MyReal Rho = 0.0;
    MyReal U, V, W, N1, N2, XMAX, YMAX, ZMAX;
    XMAX = 1;
    YMAX = 1;
    ZMAX = 1;

    //      Velocity of the wall depends                Velocity of the wall is read      Velocity of the
    //      on the velocity of the particles            from SetUpData.ini                wall is 0

    //            U = (Cells + i)->U;          |          U = (Cells + i)->Uo;          |  U = 0;
    //            V = (Cells + i)->V;          |          V = (Cells + i)->Vo;          |  V = 0;
    //            W = (Cells + i)->W;          |          W = (Cells + i)->Wo;          |  W = 0;

    U = 0;
    V = 0;
    W = 0;


    if ((Cells + i)->Boundary == 1) // WALL BOUNDARY
    {
        if ((Cells + i)->CoordZ == 0) // INLET PLANE xy, z=0
        {
            //Unknowns: Rho, f5, f11, f12, f15 and f16

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[1] + (Cells + i)->F[2] + (Cells + i)->F[3] + (Cells + i)->F[4] +
                   (Cells + i)->F[7] + (Cells + i)->F[8] + (Cells + i)->F[9] + (Cells + i)->F[10] +
                   2 * ((Cells + i)->F[6] + (Cells + i)->F[13] + (Cells + i)->F[14] + (Cells + i)->F[17] +
                        (Cells + i)->F[18])) / (1 - W);

            N1 = 0.5 * ((Cells + i)->F[1] + (Cells + i)->F[7] + (Cells + i)->F[9] - (Cells + i)->F[2] - (Cells + i)->F[8] -
                        (Cells + i)->F[10]) - Rho * U / 3;
            N2 = 0.5 * ((Cells + i)->F[3] + (Cells + i)->F[7] + (Cells + i)->F[8] - (Cells + i)->F[4] - (Cells + i)->F[9] -
                        (Cells + i)->F[10]) - Rho * V / 3;

            (Cells + i)->F[5] = (Cells + i)->F[6] + (Rho * W) / 3;
            (Cells + i)->F[11] = (Cells + i)->F[14] + Rho / 6 * (W + U) - N1;
            (Cells + i)->F[12] = (Cells + i)->F[13] + Rho / 6 * (W - U) + N1;
            (Cells + i)->F[15] = (Cells + i)->F[18] + Rho / 6 * (W + V) - N2;
            (Cells + i)->F[16] = (Cells + i)->F[17] + Rho / 6 * (W - V) + N2;
        }
        else if ((Cells + i)->CoordZ == ZMAX) // INLET PLANE xy, z=zmax
        {
            //Unknowns: Rho, f6, f13, f14, f17 and f18

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[1] + (Cells + i)->F[2] + (Cells + i)->F[3] + (Cells + i)->F[4] +
                   (Cells + i)->F[7] + (Cells + i)->F[8] + (Cells + i)->F[9] + (Cells + i)->F[10] +
                   2 * ((Cells + i)->F[5] + (Cells + i)->F[11] + (Cells + i)->F[12] + (Cells + i)->F[15] +
                        (Cells + i)->F[16])) / (1 + W);

            N1 = 0.5 * ((Cells + i)->F[1] + (Cells + i)->F[7] + (Cells + i)->F[9] - (Cells + i)->F[2] - (Cells + i)->F[8] -
                        (Cells + i)->F[10]) - Rho * U / 3;
            N2 = 0.5 * ((Cells + i)->F[3] + (Cells + i)->F[7] + (Cells + i)->F[8] - (Cells + i)->F[4] - (Cells + i)->F[9] -
                        (Cells + i)->F[10]) - Rho * V / 3;

            (Cells + i)->F[6]  = (Cells + i)->F[5]  - (Rho * W) / 3;
            (Cells + i)->F[13] = (Cells + i)->F[12] + Rho / 6 * (-W + U) - N1;
            (Cells + i)->F[14] = (Cells + i)->F[11] + Rho / 6 * (-W - U) + N1;
            (Cells + i)->F[17] = (Cells + i)->F[16] + Rho / 6 * (-W + V) - N2;
            (Cells + i)->F[18] = (Cells + i)->F[15] + Rho / 6 * (-W - V) + N2;
        }

        else if ((Cells + i)->CoordX == 0) // INLET PLANE yz, x=0
        {
            //Unknowns: Rho, f1, f7, f9, f11 and f13

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[3] + (Cells + i)->F[4] + (Cells + i)->F[5] + (Cells + i)->F[6] +
                   (Cells + i)->F[15] + (Cells + i)->F[16] + (Cells + i)->F[17] + (Cells + i)->F[18] +
                   2 *
                   ((Cells + i)->F[2] + (Cells + i)->F[8] + (Cells + i)->F[10] + (Cells + i)->F[12] + (Cells + i)->F[14])) /
                  (1 + U);

            N1 = 0.5 *
                 ((Cells + i)->F[3] + (Cells + i)->F[15] + (Cells + i)->F[17] - (Cells + i)->F[4] - (Cells + i)->F[16] -
                  (Cells + i)->F[18]) - Rho * V / 3;
            N2 = 0.5 *
                 ((Cells + i)->F[5] + (Cells + i)->F[15] + (Cells + i)->F[16] - (Cells + i)->F[6] - (Cells + i)->F[17] -
                  (Cells + i)->F[18]) - Rho * W / 3;

            (Cells + i)->F[1]  = (Cells + i)->F[2]  + (Rho * U) / 3;
            (Cells + i)->F[7]  = (Cells + i)->F[10] + Rho / 6 * (U + V) - N1;
            (Cells + i)->F[9]  = (Cells + i)->F[8]  + Rho / 6 * (U - V) + N1;
            (Cells + i)->F[11] = (Cells + i)->F[14] + Rho / 6 * (U + W) - N2;
            (Cells + i)->F[13] = (Cells + i)->F[12] + Rho / 6 * (U - W) + N2;
        }

        else if ((Cells + i)->CoordX == XMAX) // INLET PLANE yz, x=xmax
        {
            //Unknowns: Rho, f2, f8, f10, f12 and f14

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[3] + (Cells + i)->F[4] + (Cells + i)->F[5] + (Cells + i)->F[6] +
                   (Cells + i)->F[15] + (Cells + i)->F[16] + (Cells + i)->F[17] + (Cells + i)->F[18] +
                   2 *
                   ((Cells + i)->F[1] + (Cells + i)->F[7] + (Cells + i)->F[9] + (Cells + i)->F[11] + (Cells + i)->F[13])) /
                  (1 - U);

            N1 = 0.5 *
                 ((Cells + i)->F[3] + (Cells + i)->F[15] + (Cells + i)->F[17] - (Cells + i)->F[4] - (Cells + i)->F[16] -
                  (Cells + i)->F[18]) - Rho * V / 3;
            N2 = 0.5 *
                 ((Cells + i)->F[5] + (Cells + i)->F[15] + (Cells + i)->F[16] - (Cells + i)->F[6] - (Cells + i)->F[17] -
                  (Cells + i)->F[18]) - Rho * W / 3;

            (Cells + i)->F[2]  = (Cells + i)->F[1]  - (Rho * U) / 3;
            (Cells + i)->F[8]  = (Cells + i)->F[9]  + Rho / 6 * (-U + V) - N1;
            (Cells + i)->F[10] = (Cells + i)->F[7]  + Rho / 6 * (-U - V) + N1;
            (Cells + i)->F[12] = (Cells + i)->F[13] + Rho / 6 * (-U + W) - N2;
            (Cells + i)->F[14] = (Cells + i)->F[11] + Rho / 6 * (-U - W) + N2;
        }

        else if ((Cells + i)->CoordY == 0) // INLET PLANE xz, y=0
        {
            //Unknowns: Rho, f3, f7, f8, f15 and f17

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[1] + (Cells + i)->F[2] + (Cells + i)->F[5] + (Cells + i)->F[6] +
                   (Cells + i)->F[11] + (Cells + i)->F[12] + (Cells + i)->F[13] + (Cells + i)->F[14] +
                   2 *
                   ((Cells + i)->F[4] + (Cells + i)->F[9] + (Cells + i)->F[10] + (Cells + i)->F[16] + (Cells + i)->F[18])) /
                  (1 - V);

            N1 = 0.5 *
                 ((Cells + i)->F[1] + (Cells + i)->F[11] + (Cells + i)->F[13] - (Cells + i)->F[2] - (Cells + i)->F[12] -
                  (Cells + i)->F[14]) - Rho * U / 3;
            N2 = 0.5 *
                 ((Cells + i)->F[5] + (Cells + i)->F[11] + (Cells + i)->F[12] - (Cells + i)->F[6] - (Cells + i)->F[13] -
                  (Cells + i)->F[14]) - Rho * W / 3;

            (Cells + i)->F[3]  = (Cells + i)->F[4]  + (Rho * V) / 3;
            (Cells + i)->F[7]  = (Cells + i)->F[10] + Rho / 6 * (V + U) - N1;
            (Cells + i)->F[8]  = (Cells + i)->F[9]  + Rho / 6 * (V - U) + N1;
            (Cells + i)->F[15] = (Cells + i)->F[18] + Rho / 6 * (V + W) - N2;
            (Cells + i)->F[17] = (Cells + i)->F[16] + Rho / 6 * (V - W) + N2;
        }

        else if ((Cells + i)->CoordY == YMAX) // INLET PLANE xz, y=ymax
        {
            //Unknowns: Rho, f4, f9, f10, f16 and f18

            Rho = ((Cells + i)->F[0] + (Cells + i)->F[1] + (Cells + i)->F[2] + (Cells + i)->F[5] + (Cells + i)->F[6] +
                   (Cells + i)->F[11] + (Cells + i)->F[12] + (Cells + i)->F[13] + (Cells + i)->F[14] +
                   2 *
                   ((Cells + i)->F[3] + (Cells + i)->F[7] + (Cells + i)->F[8] + (Cells + i)->F[15] + (Cells + i)->F[17])) /
                  (1 + V);

            N1 = 0.5 *
                 ((Cells + i)->F[1] + (Cells + i)->F[11] + (Cells + i)->F[13] - (Cells + i)->F[2] - (Cells + i)->F[12] -
                  (Cells + i)->F[14]) - Rho * U / 3;
            N2 = 0.5 *
                 ((Cells + i)->F[5] + (Cells + i)->F[11] + (Cells + i)->F[12] - (Cells + i)->F[6] - (Cells + i)->F[13] -
                  (Cells + i)->F[14]) - Rho * W / 3;

            (Cells + i)->F[4]  = (Cells + i)->F[3]  - (Rho * V) / 3;
            (Cells + i)->F[9]  = (Cells + i)->F[8]  + Rho / 6 * (-V + U) - N1;
            (Cells + i)->F[10] = (Cells + i)->F[7]  + Rho / 6 * (-V - U) + N1;
            (Cells + i)->F[16] = (Cells + i)->F[17] + Rho / 6 * (-V + W) - N2;
            (Cells + i)->F[18] = (Cells + i)->F[15] + Rho / 6 * (-V - W) + N2;
        }
    }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////// Curved wall boundary treatment //////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void CurvedWallBoundaries(CellProps *Cells, int j, int i, int* opp)
{
    int k=0;
    /*for(k=0;k<9;k++)
    {
        if ((Cells+j*(*n)+i)->BC_ID[k]==1) // if wall
        {
            if ((Cells+j*(*n)+i)->Q[k]<0.5) // if the distance from the boundary is less than 0.5?
            {
                (Cells+j*(*n)+i)->F[opp[k]] = 2*(Cells+j*(*n)+i)->Q[k]*(Cells+j*(*n)+i)->METAF[k]
                                              +(1-2*(Cells+j*(*n)+i)->Q[k])*(Cells+j*(*n)+i)->Fneighbours[k];
            }
            else
            {
                (Cells+j*(*n)+i)->F[opp[k]] = (Cells+j*(*n)+i)->METAF[k]/2/(Cells+j*(*n)+i)->Q[k]
                                              +(2*(Cells+j*(*n)+i)->Q[k]-1)/(2*(Cells+j*(*n)+i)->Q[k])*(Cells+j*(*n)+i)->METAF[opp[k]];
            }
        }
    }*/
}

void EdgeBC(CellProps *Cells, int i)
{
    MyReal N, XMAX, YMAX, ZMAX, Fburied;
    XMAX = 1;
    YMAX = 1;
    ZMAX = 1;
    Fburied = 0;

    if ((Cells + i)->Boundary == 3)
    {
        if ((Cells + i)->CoordX == 0 && (Cells + i)->CoordZ == 0) // EDGE 1: BOTTOM LEFT
        {
            //Unknowns: f1, f7, f9    xy-plane
            //          f5, f15, f16  yz-plane
            //          f11           xz-plane
            //          f12, f13      Buried Links

            (Cells + i)->F[1]  = (Cells + i)->F[2];
            (Cells + i)->F[5]  = (Cells + i)->F[6];
            (Cells + i)->F[11] = (Cells + i)->F[14];

            N = 0.25 * ((Cells + i)->F[3] - (Cells + i)->F[4]);

            (Cells + i)->F[7]  = (Cells + i)->F[10] - N;
            (Cells + i)->F[9]  = (Cells + i)->F[8]  + N;
            (Cells + i)->F[15] = (Cells + i)->F[18] - N;
            (Cells + i)->F[16] = (Cells + i)->F[17] + N;

            for (int j = 1; j < 12; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 14; j < 19; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[12] = Fburied / 22;
            (Cells + i)->F[13] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[12];
        }

        else if ((Cells + i)->CoordX == XMAX && (Cells + i)->CoordZ == 0) // EDGE 2: BOTTOM RIGHT
        {
            //Unknowns: f2, f8, f10   xy-plane
            //          f5, f15, f16  yz-plane
            //          f12           xz-plane
            //          f11, f14      Buried Links

            (Cells + i)->F[2]  = (Cells + i)->F[1];
            (Cells + i)->F[5]  = (Cells + i)->F[6];
            (Cells + i)->F[12] = (Cells + i)->F[13];

            N = 0.25 * ((Cells + i)->F[3] - (Cells + i)->F[4]);

            (Cells + i)->F[8]  = (Cells + i)->F[9]  - N;
            (Cells + i)->F[10] = (Cells + i)->F[7]  + N;
            (Cells + i)->F[15] = (Cells + i)->F[18] - N;
            (Cells + i)->F[16] = (Cells + i)->F[17] + N;

            for (int j = 1; j < 11; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 12; j < 14; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 15; j < 19; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[11] = Fburied / 22;
            (Cells + i)->F[14] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[11];
        }

        else if ((Cells + i)->CoordY == 0 && (Cells + i)->CoordZ == 0) // EDGE 3: BOTTOM FRONT
        {
            //Unknowns: f3, f7, f8    xy-plane
            //          f5, f11, f12  xz-plane
            //          f15           yz-plane
            //          f16, f17      Buried Links

            (Cells + i)->F[3]  = (Cells + i)->F[4];
            (Cells + i)->F[5]  = (Cells + i)->F[6];
            (Cells + i)->F[15] = (Cells + i)->F[18];

            N = 0.25 * ((Cells + i)->F[1] - (Cells + i)->F[2]);

            (Cells + i)->F[7]  = (Cells + i)->F[10] - N;
            (Cells + i)->F[8]  = (Cells + i)->F[9]  + N;
            (Cells + i)->F[11] = (Cells + i)->F[14] - N;
            (Cells + i)->F[12] = (Cells + i)->F[13] + N;

            for (int j = 1; j < 16; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            Fburied = Fburied + (Cells + i)->F[18];

            (Cells + i)->F[16] = Fburied / 22;
            (Cells + i)->F[17] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[16];
        }

        else if ((Cells + i)->CoordY == YMAX && (Cells + i)->CoordZ == 0) // EDGE 4: BOTTOM BACK
        {
            //Unknowns: f4, f9, f10   xy-plane
            //          f5, f11, f12  xz-plane
            //          f16           yz-plane
            //          f15, f18      Buried Links

            (Cells + i)->F[4]  = (Cells + i)->F[3];
            (Cells + i)->F[5]  = (Cells + i)->F[6];
            (Cells + i)->F[16] = (Cells + i)->F[17];

            N = 0.25 * ((Cells + i)->F[1] - (Cells + i)->F[2]);

            (Cells + i)->F[9]  = (Cells + i)->F[8]  - N;
            (Cells + i)->F[10] = (Cells + i)->F[7]  + N;
            (Cells + i)->F[11] = (Cells + i)->F[14] - N;
            (Cells + i)->F[12] = (Cells + i)->F[13] + N;

            for (int j = 1; j < 15; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 16; j < 18; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[15] = Fburied / 22;
            (Cells + i)->F[18] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[15];

        }

        else if ((Cells + i)->CoordX == 0 && (Cells + i)->CoordZ == ZMAX) // EDGE 5: TOP LEFT
        {
            //Unknowns: f1, f7, f9    xy-plane
            //          f6, f17, f18  yz-plane
            //          f13           xz-plane
            //          f11, f14      Buried Links

            (Cells + i)->F[1]  = (Cells + i)->F[2];
            (Cells + i)->F[6]  = (Cells + i)->F[5];
            (Cells + i)->F[13] = (Cells + i)->F[12];

            N = 0.25 * ((Cells + i)->F[3] - (Cells + i)->F[4]);

            (Cells + i)->F[7]  = (Cells + i)->F[10] - N;
            (Cells + i)->F[9]  = (Cells + i)->F[8]  + N;
            (Cells + i)->F[17] = (Cells + i)->F[16] - N;
            (Cells + i)->F[18] = (Cells + i)->F[15] + N;

            for (int j = 1; j < 11; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 12; j < 14; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 15; j < 19; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[11] = Fburied / 22;
            (Cells + i)->F[14] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[11];
        }

        else if ((Cells + i)->CoordX == XMAX && (Cells + i)->CoordZ == ZMAX) // EDGE 6: TOP RIGHT
        {
            //Unknowns: f2, f8, f10   xy-plane
            //          f6, f17, f18  yz-plane
            //          f14           xz-plane
            //          f12, f13      Buried Links

            (Cells + i)->F[2]  = (Cells + i)->F[1];
            (Cells + i)->F[6]  = (Cells + i)->F[5];
            (Cells + i)->F[14] = (Cells + i)->F[11];

            N = 0.25 * ((Cells + i)->F[3] - (Cells + i)->F[4]);

            (Cells + i)->F[8]  = (Cells + i)->F[9]  - N;
            (Cells + i)->F[10] = (Cells + i)->F[7]  + N;
            (Cells + i)->F[17] = (Cells + i)->F[16] - N;
            (Cells + i)->F[18] = (Cells + i)->F[15] + N;

            for (int j = 1; j < 12; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 14; j < 19; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[12] = Fburied / 22;
            (Cells + i)->F[13] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[12];
        }

        else if ((Cells + i)->CoordY == 0 && (Cells + i)->CoordZ == ZMAX) // EDGE 7: TOP FRONT
        {
            //Unknowns: f3, f7, f8    xy-plane
            //          f6, f13, f14  xz-plane
            //          f17           yz-plane
            //          f15, f18      Buried Links

            (Cells + i)->F[3]  = (Cells + i)->F[4];
            (Cells + i)->F[6]  = (Cells + i)->F[5];
            (Cells + i)->F[17] = (Cells + i)->F[16];

            N = 0.25 * ((Cells + i)->F[1] - (Cells + i)->F[2]);

            (Cells + i)->F[7]  = (Cells + i)->F[10] - N;
            (Cells + i)->F[8]  = (Cells + i)->F[9]  + N;
            (Cells + i)->F[13] = (Cells + i)->F[12] - N;
            (Cells + i)->F[14] = (Cells + i)->F[11] + N;

            for (int j = 1; j < 15; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 16; j < 18; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[15] = Fburied / 22;
            (Cells + i)->F[18] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[15];
        }

        else if ((Cells + i)->CoordY == YMAX && (Cells + i)->CoordZ == ZMAX) // EDGE 8: TOP BACK
        {
            //Unknowns: f4, f9, f10   xy-plane
            //          f6, f13, f14  xz-plane (ME HAN HECHO LA 13-14)
            //          f18           yz-plane
            //          f16, f17      Buried Links

            (Cells + i)->F[4]  = (Cells + i)->F[3];
            (Cells + i)->F[6]  = (Cells + i)->F[5];
            (Cells + i)->F[18] = (Cells + i)->F[15];

            N = 0.25 * ((Cells + i)->F[1] - (Cells + i)->F[2]);

            (Cells + i)->F[9]  = (Cells + i)->F[8]  - N;
            (Cells + i)->F[10] = (Cells + i)->F[7]  + N;
            (Cells + i)->F[13] = (Cells + i)->F[12] - N;
            (Cells + i)->F[14] = (Cells + i)->F[11] + N;

            for (int j = 1; j < 16; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            Fburied = Fburied + (Cells + i)->F[18];

            (Cells + i)->F[16] = Fburied / 22;
            (Cells + i)->F[17] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[16];
        }

        else if ((Cells + i)->CoordY == 0 && (Cells + i)->CoordX == 0) // EDGE 9: LEFT FRONT
        {
            //Unknowns: f1, f11, f13  xz-plane
            //          f3, f15, f17  yz-plane
            //          f7            yz-plane
            //          f8, f9        Buried Links

            (Cells + i)->F[1] = (Cells + i)->F[2];
            (Cells + i)->F[3] = (Cells + i)->F[4];
            (Cells + i)->F[7] = (Cells + i)->F[10];

            N = 0.25 * ((Cells + i)->F[5] - (Cells + i)->F[6]);

            (Cells + i)->F[11] = (Cells + i)->F[14] - N;
            (Cells + i)->F[13] = (Cells + i)->F[12] + N;
            (Cells + i)->F[15] = (Cells + i)->F[18] - N;
            (Cells + i)->F[17] = (Cells + i)->F[16] + N;

            for (int j = 1; j < 8; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 10; j < 19; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[8] = Fburied / 22;
            (Cells + i)->F[9] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[8];
        }

        else if ((Cells + i)->CoordY == 0 && (Cells + i)->CoordX == XMAX) // EDGE 10: RIGHT FRONT
        {
            //Unknowns: f2, f12, f14  xz-plane
            //          f3, f15, f17  yz-plane
            //          f8            yz-plane
            //          f7, f10       Buried Links

            (Cells + i)->F[2] = (Cells + i)->F[1];
            (Cells + i)->F[3] = (Cells + i)->F[4];
            (Cells + i)->F[8] = (Cells + i)->F[9];

            N = 0.25 * ((Cells + i)->F[5] - (Cells + i)->F[6]);

            (Cells + i)->F[12] = (Cells + i)->F[13] - N;
            (Cells + i)->F[14] = (Cells + i)->F[11] + N;
            (Cells + i)->F[15] = (Cells + i)->F[18] - N;
            (Cells + i)->F[17] = (Cells + i)->F[16] + N;

            for (int j = 1; j < 7; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 8; j < 10; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 11; j < 19; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[7] = Fburied / 22;
            (Cells + i)->F[10] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[7];
        }

        else if ((Cells + i)->CoordY == YMAX && (Cells + i)->CoordX == 0) // EDGE 11: LEFT BACK
        {
            //Unknowns: f1, f11, f13  xz-plane
            //          f4, f16, f18  yz-plane
            //          f9            yz-plane
            //          f7, f10       Buried Links

            (Cells + i)->F[1] = (Cells + i)->F[2];
            (Cells + i)->F[4] = (Cells + i)->F[3];
            (Cells + i)->F[9] = (Cells + i)->F[8];

            N = 0.25 * ((Cells + i)->F[5] - (Cells + i)->F[6]);

            (Cells + i)->F[11] = (Cells + i)->F[14] - N;
            (Cells + i)->F[13] = (Cells + i)->F[12] + N;
            (Cells + i)->F[16] = (Cells + i)->F[17] - N;
            (Cells + i)->F[18] = (Cells + i)->F[15] + N;

            for (int j = 1; j < 7; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 8; j < 10; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 11; j < 19; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[7] = Fburied / 22;
            (Cells + i)->F[10] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[7];
        }

        else if ((Cells + i)->CoordY == YMAX && (Cells + i)->CoordX == XMAX) // EDGE 12: RIGHT BACK
        {
            //Unknowns: f2, f12, f14  xz-plane
            //          f4, f16, f18  yz-plane
            //          f10            yz-plane
            //          f8, f9        Buried Links

            (Cells + i)->F[2] = (Cells + i)->F[1];
            (Cells + i)->F[4] = (Cells + i)->F[3];
            (Cells + i)->F[10] = (Cells + i)->F[7];

            N = 0.25 * ((Cells + i)->F[5] - (Cells + i)->F[6]);

            (Cells + i)->F[12] = (Cells + i)->F[13] - N;
            (Cells + i)->F[14] = (Cells + i)->F[11] + N;
            (Cells + i)->F[16] = (Cells + i)->F[17] - N;
            (Cells + i)->F[18] = (Cells + i)->F[15] + N;

            for (int j = 1; j < 8; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 10; j < 19; ++j)
            {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[8] = Fburied / 22;
            (Cells + i)->F[9] = Fburied / 22;
            (Cells + i)->F[0]  = 12 * (Cells + i)->F[8];
        }
    }
}

void CornerBC(CellProps *Cells, int i) {
    if ((Cells + i)->Boundary == 4) {
        MyReal XMAX, YMAX, ZMAX, Fburied;
        XMAX = 1;
        YMAX = 1;
        ZMAX = 1;
        Fburied = 0;

        if ((Cells + i)->CoordX == 0 && (Cells + i)->CoordY == 0 && (Cells + i)->CoordZ == 0) {
            (Cells + i)->F[1] = (Cells + i)->F[2];
            (Cells + i)->F[3] = (Cells + i)->F[4];
            (Cells + i)->F[5] = (Cells + i)->F[6];
            (Cells + i)->F[7] = (Cells + i)->F[10];
            (Cells + i)->F[11] = (Cells + i)->F[14];
            (Cells + i)->F[15] = (Cells + i)->F[18];

            for (int j = 1; j < 8; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 10; j < 12; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 14; j < 16; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            Fburied = Fburied + (Cells + i)->F[18];

            (Cells + i)->F[8] = Fburied / 18;
            (Cells + i)->F[9] = Fburied / 18;
            (Cells + i)->F[12] = Fburied / 18;
            (Cells + i)->F[13] = Fburied / 18;
            (Cells + i)->F[16] = Fburied / 18;
            (Cells + i)->F[17] = Fburied / 18;
            (Cells + i)->F[0] = 12 * (Cells + i)->F[8];

        }

        else if ((Cells + i)->CoordX == XMAX && (Cells + i)->CoordY == 0 && (Cells + i)->CoordZ == 0) {
            (Cells + i)->F[2] = (Cells + i)->F[1];
            (Cells + i)->F[3] = (Cells + i)->F[4];
            (Cells + i)->F[5] = (Cells + i)->F[6];
            (Cells + i)->F[8] = (Cells + i)->F[9];
            (Cells + i)->F[12] = (Cells + i)->F[13];
            (Cells + i)->F[15] = (Cells + i)->F[18];

            for (int j = 1; j < 7; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 8; j < 10; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 12; j < 14; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            Fburied = Fburied + (Cells + i)->F[15];
            Fburied = Fburied + (Cells + i)->F[18];

            (Cells + i)->F[7] = Fburied / 18;
            (Cells + i)->F[10] = Fburied / 18;
            (Cells + i)->F[11] = Fburied / 18;
            (Cells + i)->F[14] = Fburied / 18;
            (Cells + i)->F[16] = Fburied / 18;
            (Cells + i)->F[17] = Fburied / 18;
            (Cells + i)->F[0] = 12 * (Cells + i)->F[7];

        }

        else if ((Cells + i)->CoordX == 0 && (Cells + i)->CoordY == YMAX && (Cells + i)->CoordZ == 0) {
            (Cells + i)->F[1] = (Cells + i)->F[2];
            (Cells + i)->F[4] = (Cells + i)->F[3];
            (Cells + i)->F[5] = (Cells + i)->F[6];
            (Cells + i)->F[9] = (Cells + i)->F[8];
            (Cells + i)->F[11] = (Cells + i)->F[14];
            (Cells + i)->F[16] = (Cells + i)->F[17];

            for (int j = 1; j < 7; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 8; j < 10; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            Fburied = Fburied + (Cells + i)->F[11];
            Fburied = Fburied + (Cells + i)->F[14];

            for (int j = 16; j < 18; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[7] = Fburied / 18;
            (Cells + i)->F[10] = Fburied / 18;
            (Cells + i)->F[12] = Fburied / 18;
            (Cells + i)->F[13] = Fburied / 18;
            (Cells + i)->F[15] = Fburied / 18;
            (Cells + i)->F[18] = Fburied / 18;
            (Cells + i)->F[0] = 12 * (Cells + i)->F[7];

        }

        else if ((Cells + i)->CoordX == XMAX && (Cells + i)->CoordY == YMAX && (Cells + i)->CoordZ == 0) {
            (Cells + i)->F[2] = (Cells + i)->F[1];
            (Cells + i)->F[4] = (Cells + i)->F[3];
            (Cells + i)->F[5] = (Cells + i)->F[6];
            (Cells + i)->F[10] = (Cells + i)->F[7];
            (Cells + i)->F[12] = (Cells + i)->F[13];
            (Cells + i)->F[16] = (Cells + i)->F[17];

            for (int j = 1; j < 8; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            Fburied = Fburied + (Cells + i)->F[10];

            for (int j = 12; j < 14; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 16; j < 18; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[8] = Fburied / 18;
            (Cells + i)->F[9] = Fburied / 18;
            (Cells + i)->F[11] = Fburied / 18;
            (Cells + i)->F[14] = Fburied / 18;
            (Cells + i)->F[15] = Fburied / 18;
            (Cells + i)->F[18] = Fburied / 18;
            (Cells + i)->F[0] = 12 * (Cells + i)->F[8];

        }

        else if ((Cells + i)->CoordX == 0 && (Cells + i)->CoordY == 0 && (Cells + i)->CoordZ == ZMAX) {
            (Cells + i)->F[1] = (Cells + i)->F[2];
            (Cells + i)->F[3] = (Cells + i)->F[4];
            (Cells + i)->F[6] = (Cells + i)->F[5];
            (Cells + i)->F[7] = (Cells + i)->F[10];
            (Cells + i)->F[13] = (Cells + i)->F[12];
            (Cells + i)->F[16] = (Cells + i)->F[17];

            for (int j = 1; j < 8; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            Fburied = Fburied + (Cells + i)->F[10];

            for (int j = 12; j < 14; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 16; j < 18; ++i)
            {
              Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[8] = Fburied / 18;
            (Cells + i)->F[9] = Fburied / 18;
            (Cells + i)->F[11] = Fburied / 18;
            (Cells + i)->F[14] = Fburied / 18;
            (Cells + i)->F[15] = Fburied / 18;
            (Cells + i)->F[18] = Fburied / 18;
            (Cells + i)->F[0] = 12 * (Cells + i)->F[8];

        }

        else if ((Cells + i)->CoordX == XMAX && (Cells + i)->CoordY == 0 && (Cells + i)->CoordZ == ZMAX) {
            (Cells + i)->F[2] = (Cells + i)->F[1];
            (Cells + i)->F[3] = (Cells + i)->F[4];
            (Cells + i)->F[6] = (Cells + i)->F[5];
            (Cells + i)->F[8] = (Cells + i)->F[9];
            (Cells + i)->F[14] = (Cells + i)->F[11];
            (Cells + i)->F[17] = (Cells + i)->F[16];

            for (int j = 1; j < 7; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 8; j < 10; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            Fburied = Fburied + (Cells + i)->F[11];
            Fburied = Fburied + (Cells + i)->F[14];

            for (int j = 16; j < 18; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            (Cells + i)->F[7] = Fburied / 18;
            (Cells + i)->F[10] = Fburied / 18;
            (Cells + i)->F[12] = Fburied / 18;
            (Cells + i)->F[13] = Fburied / 18;
            (Cells + i)->F[15] = Fburied / 18;
            (Cells + i)->F[18] = Fburied / 18;
            (Cells + i)->F[0] = 12 * (Cells + i)->F[7];

        }

        else if ((Cells + i)->CoordX == 0 && (Cells + i)->CoordY == YMAX && (Cells + i)->CoordZ == ZMAX) {
            (Cells + i)->F[1] = (Cells + i)->F[2];
            (Cells + i)->F[4] = (Cells + i)->F[3];
            (Cells + i)->F[6] = (Cells + i)->F[5];
            (Cells + i)->F[9] = (Cells + i)->F[8];
            (Cells + i)->F[13] = (Cells + i)->F[12];
            (Cells + i)->F[18] = (Cells + i)->F[15];

            for (int j = 1; j < 7; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 8; j < 10; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 12; j < 14; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            Fburied = Fburied + (Cells + i)->F[15];
            Fburied = Fburied + (Cells + i)->F[18];

            (Cells + i)->F[7] = Fburied / 18;
            (Cells + i)->F[10] = Fburied / 18;
            (Cells + i)->F[11] = Fburied / 18;
            (Cells + i)->F[14] = Fburied / 18;
            (Cells + i)->F[16] = Fburied / 18;
            (Cells + i)->F[17] = Fburied / 18;
            (Cells + i)->F[0] = 12 * (Cells + i)->F[7];

        }

        else if ((Cells + i)->CoordX == XMAX && (Cells + i)->CoordY == YMAX && (Cells + i)->CoordZ == ZMAX) {
            (Cells + i)->F[2] = (Cells + i)->F[1];
            (Cells + i)->F[4] = (Cells + i)->F[3];
            (Cells + i)->F[6] = (Cells + i)->F[5];
            (Cells + i)->F[10] = (Cells + i)->F[7];
            (Cells + i)->F[14] = (Cells + i)->F[11];
            (Cells + i)->F[18] = (Cells + i)->F[15];

            for (int j = 1; j < 8; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 10; j < 12; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            for (int j = 14; j < 16; ++j) {
                Fburied = Fburied + (Cells + i)->F[j];
            }

            Fburied = Fburied + (Cells + i)->F[18];

            (Cells + i)->F[8] = Fburied / 18;
            (Cells + i)->F[9] = Fburied / 18;
            (Cells + i)->F[12] = Fburied / 18;
            (Cells + i)->F[13] = Fburied / 18;
            (Cells + i)->F[16] = Fburied / 18;
            (Cells + i)->F[17] = Fburied / 18;
            (Cells + i)->F[0] = 12 * (Cells + i)->F[8];

        }

    }
}


void generalWall(CellProps* Cells, int i) {

    double v[3];
    int wall_id = 0;


    if ((Cells + i)->Boundary == 2 || (Cells + i)->Boundary == 1){
        if ((Cells + i)->CoordX == 0)
        {
            wall_id = 0;
        }
        else if ((Cells + i)->CoordX == 1)
        {
            wall_id = 1;
        }
        else if ((Cells + i)->CoordY == 0)
        {
            wall_id = 2;
        }
        else if ((Cells + i)->CoordY == 1)
        {
            wall_id = 3;
        }
        else if ((Cells + i)->CoordZ == 0)
        {
            wall_id = 4;
        }
        else if ((Cells + i)->CoordZ == 1)
        {
            wall_id = 5;
        }

        if ((Cells + i)->Boundary == 2) {// INLET BOUNDARY
            v[0] = (Cells + i)->Uo;
            v[1] = (Cells + i)->Vo;
            v[2] = (Cells + i)->Wo;
        }
        else if ((Cells + i)->Boundary == 1) {
            v[0] = 0;
            v[1] = 0;
            v[2] = 0;
        }
        getWall(Cells,i,wall_id,v);
    }
}

void getWall(CellProps *Cells, int i, int wall_ID, double* v) {



    int cj[3];
    int ck[3];
    int j;
    double temp = 0;

    j = opp[j_wall_unknown[wall_ID][0]];

    cj[0] = cx[j];
    cj[1] = cy[j];
    cj[2] = cz[j];

    Cells[i].F[opp[j]] = Cells[i].F[j] - 1./3.*Cells[i].Rho * multVec(cj,v);

    //f(opp(j)) = f(j) - rho/6 *c(j)*V - rho/3 * t(j)*V + 1/2 * SUM k=0:18 (fk * (t(j)*c(k)) * (1 - |c(k)*n|))
    for (int l = 1; l<5; l++) {
        //Ci
        j = opp[j_wall_unknown[wall_ID][l]];

        //c(j)
        cj[0] = cx[j];
        cj[1] = cy[j];
        cj[2] = cz[j];
        getTan(cj, norm[wall_ID], tang);

        for (int k = 0; k<19; k++) {
            ck[0] = cx[k];
            ck[1] = cy[k];
            ck[2] = cz[k];
            temp += Cells[i].F[k]*multVec(tang,ck)*
                    (1 - fabs(multVec(ck,norm[wall_ID])));
        }

        Cells[i].F[opp[j]] = Cells[i].F[j]
                             - Cells[i].Rho/6 * multVec(cj,v)
                             - Cells[i].Rho/3 * multVec(tang,v)
                             + 0.5 * temp;

        /*if (Cells[i].CoordZ == 1 && l == 1 && i == 6260)
        {
            printf("i = %i\n",i);
            printf("Wall ID = %i\n",wall_ID);
            printf("cj =(%i,%i,%i)\n",cj[0],cj[1],cj[2]);
            printf("norm =(%i,%i,%i)\n",norm[wall_ID][0],norm[wall_ID][1],norm[wall_ID][2]);
            printf("tang =(%i,%i,%i)\n",tang[0],tang[1],tang[2]);
            printf("v =(%f,%f,%f)\n",v[0],v[1],v[2]);
            printf("cj * v = %f\n",multVec(cj,v));
            printf("v * v = %f\n",multVec(v,v));
        }*/

    }



}


void getTan(const int* c, const int* n, int* t) {
    int cn = c[0]*n[0] + c[1]*n[1] + c[2]*n[2];

    for (int i = 0; i<3; i++) {
        t[i] = c[i] - cn*n[i];
    }

}
