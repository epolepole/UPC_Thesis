#include <stdio.h>   // for calloc();
#include "math.h"
#include <upc_strict.h>                 // Required for UPC
#include <CellFunctions.h>

#include "ShellFunctions.h"


void WriteResults(char* OutputFile, int* postproc_prog)
{
    FILE * fp1;                 // file pointer to output file
    fp1=fopen(OutputFile, "w"); // open file
    //upc_barrier;


    switch(*postproc_prog)
    {
        case 1: // ParaView
            //fprintf(fp1, "X Column,Y Column,Z Column,u,v,w,vel_mag,"\
                    "f00,f01,f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18,"\
                    "f_m00,f_m01,f_m02,f_m03,f_m04,f_m05,f_m06,f_m07,f_m08,f_m09,f_m10,f_m11,f_m12,f_m13,f_m14,f_m15,f_m16,f_m17,f_m18,"\
                    "rho,press,fluid,ThID\n");
            /*for(i = 0; i < NODES; i++)
            {
                fprintf(fp1, "%f, %f, %f, %f, %f, %f, %f,",
                        (WCells+i)->CoordX, // x
                        (WCells+i)->CoordY, // y
                        (WCells+i)->CoordZ, // z
                        (WCells+i)->U,      // u
                        (WCells+i)->V,      // v
                        (WCells+i)->W,      // w
                        sqrt(pow((WCells+i)->U,2)+pow((WCells+i)->V,2)+pow((WCells+i)->W,2)));
                for (int j = 0; j<19; j++) {
                    fprintf(fp1," %f,",(WCells+i)->F[j]);
                }
                for (int j = 0; j<19; j++) {
                    fprintf(fp1," %f,",(WCells+i)->METAF[j]);
                }
                fprintf(fp1," %f, %f, %d, %d\n",
                        (WCells+i)->Rho,    // density
                        ((WCells+i)->Rho)/3,  // pressure
                        (WCells+i)->Fluid, // fluid or solid
                        (WCells+i)->ThreadNumber);  */
            fprintf(fp1, "Points:0,Points:1,Points:2,u,v,w,vel_mag,"\
                        // "f00,f01,f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18,"
                    "rho\n");
            for (int k = 0; k < NN; k++) {
                for (int j = 0; j < NM; j++) {
                    for (int i = 0; i < NL; i++) {
                        int tX[3] = {i,j,k}; //Global index

                        int CX[3];  //Cube index
                        int lX[3];  //local index
                        int cubeID;
                        int locID;

                        getCubeCoords_TotIndex(&tX[0],&CX[0]);
                        getLocalIndex_TotIndex(&tX[0],&lX[0]);
                        cubeID = getCubeID(CX[0],CX[1],CX[2]);
                        locID = lID(lX[0],lX[1],lX[2]);



                        int pos = CELL_TOT_SIZE*cubeID + locID;



                        /*int pos = i%LAT + (i/LAT)*LAT*LAT*LAT +
                                (j%LAT)*LAT + (j/LAT)*LAT*LAT*LAT*NTDX +
                                (k%LAT)*LAT*LAT + (k/LAT)*LAT*LAT*LAT*NTDX*NTDY;*/
                        //printf("pos(%i,%i,%i)=%i\n",i,j,k,pos);
                        fprintf(fp1, "%f, %f, %f, %f, %f, %f, %f,",
                                (WCells + pos)->CoordX, // x
                                (WCells + pos)->CoordY, // y
                                (WCells + pos)->CoordZ, // z
                                (WCells + pos)->U,      // u
                                (WCells + pos)->V,      // v
                                (WCells + pos)->W,      // w
                                sqrt(pow((WCells + pos)->U, 2) + pow((WCells + pos)->V, 2) +
                                     pow((WCells + pos)->W, 2)));
                        /*for (int l = 0; l<19; l++) {
                        fprintf(fp1," %f,",(WCells+pos)->F[l]);
                        }*/

                        fprintf(fp1, " %f\n", (WCells + pos)->Rho);    // density
                    }
                }
            }

            fclose(fp1);
            break;

        case 2: // TECPLOT
            fprintf(fp1, "Title = \"LBM results\"\n");
            fprintf(fp1, "Variables = \"x\",\"y\",\"z\",\"u\",\"v\",\"w\",\"u mag\",\"rho\",\"press\",\"fluid\"\n");
            fprintf(fp1, "Zone i=%d, j=%d, k=%d, f=point\n",NN,NM,NL);

            for(int i = 0; i < NODES; i++)
            {
                fprintf(fp1, "%f %f %f %f %f %f %f %f %f\n",
                        (WCells+i)->CoordX, // x
                        (WCells+i)->CoordY, // y
                        (WCells+i)->CoordZ, // z
                        (WCells+i)->U,      // u
                        (WCells+i)->V,      // v
                        (WCells+i)->W,      // w
                        sqrt(pow((WCells+i)->U,2)+pow((WCells+i)->V,2)+pow((WCells+i)->W,2)), // velocity magnitude
                        (WCells+i)->Rho,    // density
                        ((WCells+i)->Rho)/3);  // pressure
                //(WCells+i)->Fluid); // fluid or solid
            }

            fclose(fp1);
            break;
        default: break;
    }


}


void WriteBCells(char* OutputFile, int* postproc_prog)
{
    int i, j;                   // Loop variables
    FILE * fp1;                 // file pointer to output file
    fp1=fopen(OutputFile, "w"); // open file
    switch(*postproc_prog)
    {
        case 1: // ParaView
            fprintf(fp1, "x,y,z,u,v,w,vel_mag,"
                    "f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,"
                    "rho,press,ThID,\n");
            for(i=0;i<(2*THREADS*LAYER);i++)
            {
                fprintf(fp1, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %d\n",
                        (BCells+i)->CoordX, // x
                        (BCells+i)->CoordY, // y
                        (BCells+i)->CoordZ, // z
                        (BCells+i)->U,      // u
                        (BCells+i)->V,      // v
                        (BCells+i)->W,      // w
                        sqrt(pow((BCells+i)->U,2)+pow((BCells+i)->V,2)+pow((BCells+i)->W,2)),
                        (BCells+i)->Rho,    // density
                        ((BCells+i)->Rho)/3,  // pressure
                        (BCells+i)->ThreadNumber);
            }

            fclose(fp1);
            break;

        case 2: // TECPLOT
            fprintf(fp1, "Title = \"LBM results\"\n");
            fprintf(fp1, "Variables = \"x\",\"y\",\"z\",\"u\",\"v\",\"w\",\"u mag\",\"rho\",\"press\",\"fluid\"\n");
            fprintf(fp1, "Zone i=%d, j=%d, k=%d, f=point\n",NN,NM,NL);

            for(i=0;i<(2*THREADS*LAYER);i++)
            {
                fprintf(fp1, "%f %f %f %f %f %f %f %f %f\n",
                        (BCells+i)->CoordX, // x
                        (BCells+i)->CoordY, // y
                        (BCells+i)->CoordZ, // z
                        (BCells+i)->U,      // u
                        (BCells+i)->V,      // v
                        (BCells+i)->W,      // w
                        sqrt(pow((BCells+i)->U,2)+pow((BCells+i)->V,2)+pow((BCells+i)->W,2)), // velocity magnitude
                        (BCells+i)->Rho,    // density
                        ((BCells+i)->Rho)/3);  // pressure
                //(BCells+i)->Fluid); // fluid or solid
            }

            fclose(fp1);
            break;
        default: break;
    }
}
