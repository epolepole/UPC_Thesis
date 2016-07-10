#include <stdio.h>   // for calloc();
#include "math.h"
#include <upc_relaxed.h>                 // Required for UPC

#include "ShellFunctions.h"


void WriteResults(char* OutputFile, int* postproc_prog)
{
    int i;                      // Loop variable
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

            fprintf(fp1, "X Column,Y Column,Z Column,u,v,w,vel_mag,"\
                         "f00,f01,f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18,"\
                         "rho\n");

            for(i = 0; i < NODES; i++)
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
                
                fprintf(fp1," %f\n", (WCells+i)->Rho);    // density
            }

            fclose(fp1);
            break;

        case 2: // TECPLOT
            fprintf(fp1, "Title = \"LBM results\"\n");
            fprintf(fp1, "Variables = \"x\",\"y\",\"z\",\"u\",\"v\",\"w\",\"u mag\",\"rho\",\"press\",\"fluid\"\n");
            fprintf(fp1, "Zone i=%d, j=%d, k=%d, f=point\n",NN,NM,NL);

            for(i = 0; i < NODES; i++)
            {
                fprintf(fp1, "%f %f %f %f %f %f %f %f %f %d\n",
                        (WCells+i)->CoordX, // x
                        (WCells+i)->CoordY, // y
                        (WCells+i)->CoordZ, // z
                        (WCells+i)->U,      // u
                        (WCells+i)->V,      // v
                        (WCells+i)->W,      // w
                        sqrt(pow((WCells+i)->U,2)+pow((WCells+i)->V,2)+pow((WCells+i)->W,2)), // velocity magnitude
                        (WCells+i)->Rho,    // density
                        ((WCells+i)->Rho)/3,  // pressure
                        (WCells+i)->Fluid); // fluid or solid
            }

            fclose(fp1);
            break;
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
            fprintf(fp1, "x,y,z,u,v,w,vel_mag,rho,press,fluid,ThID\n");
            for(i=0;i<(2*THREADS*LAYER);i++)
            {
                fprintf(fp1, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %d, %d\n",
                        (BCells+i)->CoordX, // x
                        (BCells+i)->CoordY, // y
                        (BCells+i)->CoordZ, // z
                        (BCells+i)->U,      // u
                        (BCells+i)->V,      // v
                        (BCells+i)->W,      // w
                        sqrt(pow((BCells+i)->U,2)+pow((BCells+i)->V,2)+pow((BCells+i)->W,2)), // velocity magnitude
                        (BCells+i)->Rho,    // density
                        ((BCells+i)->Rho)/3,  // pressure
                        (BCells+i)->Fluid, // fluid or solid
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
                fprintf(fp1, "%f %f %f %f %f %f %f %f %f %d\n",
                        (BCells+i)->CoordX, // x
                        (BCells+i)->CoordY, // y
                        (BCells+i)->CoordZ, // z
                        (BCells+i)->U,      // u
                        (BCells+i)->V,      // v
                        (BCells+i)->W,      // w
                        sqrt(pow((BCells+i)->U,2)+pow((BCells+i)->V,2)+pow((BCells+i)->W,2)), // velocity magnitude
                        (BCells+i)->Rho,    // density
                        ((BCells+i)->Rho)/3,  // pressure
                        (BCells+i)->Fluid); // fluid or solid
            }

            fclose(fp1);
            break;
    }
}

void WriteResults2(char* OutputFile)
{
    int i,index;                     // Loop variable
    index = 2*(2*LAYER + BLOCKSIZE);
    FILE * fp1;                 // file pointer to output file
    fp1=fopen(OutputFile, "w"); // open file
    //upc_barrier;

    fprintf(fp1, "X Column,Y Column,Z Column,u,v,w,vel_mag,"\
                 "f00,f01,f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18,"\
                 "rho\n");

    for(i = 0; i < index; i++)
    {
        fprintf(fp1, "%f, %f, %f, %f, %f, %f, %f,",
                (WCells2+i)->CoordX, // x
                (WCells2+i)->CoordY, // y
                (WCells2+i)->CoordZ, // z
                (WCells2+i)->U,      // u
                (WCells2+i)->V,      // v
                (WCells2+i)->W,      // w
                sqrt(pow((WCells2+i)->U,2)+pow((WCells2+i)->V,2)+pow((WCells2+i)->W,2)));
        for (int j = 0; j<19; j++) 
        {
            fprintf(fp1," %f,",(WCells2+i)->F[j]);
        }
        
        fprintf(fp1," %f\n", (WCells2+i)->Rho);    // density
    }
    fclose(fp1);
}

