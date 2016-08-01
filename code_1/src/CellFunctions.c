
#include <math.h>                   // for sin,cos,pow... compile with -lm
#include <stdio.h>
#include <ShellFunctions.h>
#include <Iterate.h>

#include "CellFunctions.h"

/*==================================================
=========Initialization for the MRT model===========
==================================================*/
// This function fills up tm and stimv with variables
void MRTInitializer(double** tm, double** stmiv, double Omega)
{
    // RETURN THESE VALUES:
    // float tm[9][9];
    // float stmiv[9][9];

    ///////////// Declarations ////////////////
    int i, j, l;  // loop variables

    // declarations for this collision model
    float sumcc;
    float sm[9];
    float ev[9][9];
    /*
    const float a1=1./36.;
    float tminv[9][9]=
    {
        {4.*a1, -4.*a1,  4.*a1,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0},
        {4.*a1,    -a1, -2.*a1,  6.*a1, -6.*a1,    0.0,    0.0,  9.*a1,    0.0},
        {4.*a1,    -a1, -2.*a1,    0.0,    0.0,  6.*a1, -6.*a1, -9.*a1,    0.0},
        {4.*a1,    -a1, -2.*a1, -6.*a1,  6.*a1,    0.0,    0.0,  9.*a1,    0.0},
        {4.*a1,    -a1, -2.*a1,    0.0,    0.0, -6.*a1,  6.*a1, -9.*a1,    0.0},
        {4.*a1,  2.*a1,     a1,  6.*a1,  3.*a1,  6.*a1,  3.*a1,    0.0,  9.*a1},
        {4.*a1,  2.*a1,     a1, -6.*a1, -3.*a1,  6.*a1,  3.*a1,    0.0, -9.*a1},
        {4.*a1,  2.*a1,     a1, -6.*a1, -3.*a1, -6.*a1, -3.*a1,    0.0,  9.*a1},
        {4.*a1,  2.*a1,     a1,  6.*a1,  3.*a1, -6.*a1, -3.*a1,    0.0, -9.*a1}
    };*/

    float tminv[9][9] =
            {
                    {4.*(1./36), -4.*(1./36),  4.*(1./36),         0.0,         0.0,         0.0,         0.0,         0.0,         0.0},
                    {4.*(1./36),    -(1./36), -2.*(1./36),  6.*(1./36), -6.*(1./36),         0.0,         0.0,  9.*(1./36),         0.0},
                    {4.*(1./36),    -(1./36), -2.*(1./36),         0.0,         0.0,  6.*(1./36), -6.*(1./36), -9.*(1./36),         0.0},
                    {4.*(1./36),    -(1./36), -2.*(1./36), -6.*(1./36),  6.*(1./36),         0.0,         0.0,  9.*(1./36),         0.0},
                    {4.*(1./36),    -(1./36), -2.*(1./36),         0.0,         0.0, -6.*(1./36),  6.*(1./36), -9.*(1./36),         0.0},
                    {4.*(1./36),  2.*(1./36),     (1./36),  6.*(1./36),  3.*(1./36),  6.*(1./36),  3.*(1./36),         0.0,  9.*(1./36)},
                    {4.*(1./36),  2.*(1./36),     (1./36), -6.*(1./36), -3.*(1./36),  6.*(1./36),  3.*(1./36),         0.0, -9.*(1./36)},
                    {4.*(1./36),  2.*(1./36),     (1./36), -6.*(1./36), -3.*(1./36), -6.*(1./36), -3.*(1./36),         0.0,  9.*(1./36)},
                    {4.*(1./36),  2.*(1./36),     (1./36),  6.*(1./36),  3.*(1./36), -6.*(1./36), -3.*(1./36),         0.0, -9.*(1./36)}
            };



    float temp[9][9] = {
            {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0},
            {-4.,-1.,-1.,-1.,-1.,2.0,2.0,2.0,2.0},
            {4.0,-2.,-2.,-2.,-2.,1.0,1.0,1.0,1.0},
            {0.0,1.0,0.0,-1.,0.0,1.0,-1.,-1.,1.0},
            {0.0,-2.,0.0,2.0,0.0,1.0,-1.,-1.,1.0},
            {0.0,0.0,1.0,0.0,-1.,1.0,1.0,-1.,-1.},
            {0.0,0.0,-2.,0.0,2.0,1.0,1.0,-1.,-1.},
            {0.0,1.0,-1.,1.0,-1.,0.0,0.0,0.0,0.0},
            {0.0,0.0,0.0,0.0,0.0,1.0,-1.,1.0,-1.}
    };


    ///////////// Fill up variables ////////////////

    // Filling up tm
    for (i = 0; i < 9; i++)
    {
        for (j = 0; j < 9; j++)
            tm[i][j] = temp[i][j];
    }

    // Filling up stimv
    sm[0] = 1.0;
    sm[1] = 1.4;
    sm[2] = 1.4;
    sm[3] = 1.0;
    sm[4] = 1.2;
    sm[5] = 1.0;
    sm[6] = 1.2;
    sm[7] = Omega;
    sm[8] = Omega;

    for(i=0;i<9;i++)
    {
        for(j=0;j<9;j++)
        {
            sumcc=0.0;
            for(l=0;l<9;l++)
            {
                sumcc=sumcc+tminv[i][l]*tm[l][j];
            }

            ev[i][j]=sumcc;
        }
    }

    for(i=0;i<9;i++)
    {
        for(j=0;j<9;j++)
        {
            stmiv[i][j]=tminv[i][j]*sm[j];
        }
    }

} // End of function


/*==================================================
============Initialization of the cells=============
==================================================*/

void CellIni(CellProps *Cells,
             float  **Nod,
             float  **Con,
             float  Uavg,
             float  Vavg,
             float  Wavg,
             int    InletProfile,
             int    CollisionModel,
             int*   opp,
             float  rho_ini)
{
    ///////////// DECLARATION /////////////////
    int j, l, k;  // loop variables
    float Qlat[19];  // Q lattice, length of the discrete directions
    Qlat[0]  =    0   ;
    Qlat[1]  =    1   ;
    Qlat[2]  =    1   ;
    Qlat[3]  =    1   ;
    Qlat[4]  =    1   ;
    Qlat[5]  =    1   ;
    Qlat[6]  =    1   ;
    Qlat[7]  = (float)sqrt(2);
    Qlat[8]  = (float)sqrt(2);
    Qlat[9]  = (float)sqrt(2);
    Qlat[10] = (float)sqrt(2);
    Qlat[11] = (float)sqrt(2);
    Qlat[12] = (float)sqrt(2);
    Qlat[13] = (float)sqrt(2);
    Qlat[14] = (float)sqrt(2);
    Qlat[15] = (float)sqrt(2);
    Qlat[16] = (float)sqrt(2);
    Qlat[17] = (float)sqrt(2);
    Qlat[18] = (float)sqrt(2);

    int index_Cell;

    for (int i = MYTHREAD * BLOCKSIZE; i < BLOCKSIZE * (MYTHREAD + 1); i++ )
    {
        index_Cell = LAYER - MYTHREAD * BLOCKSIZE + i;

        // FIND ID of the actual cell
        (Cells + index_Cell)->ID = i;



        // FIND X, Y and Z of the actual cell
        (Cells + index_Cell)->CoordX = Nod[i][3];
        (Cells + index_Cell)->CoordY = Nod[i][4];
        (Cells + index_Cell)->CoordZ = Nod[i][5];

        // CHECK FLUID OR NOT
        //(Cells + index_Cell)->Fluid  = (int) Nod[i][6];

        // Which thread does it belongs to
        (Cells + index_Cell)->ThreadNumber = MYTHREAD;

        // INITIALIZE VARIABLEs
        (Cells + index_Cell)->U   = Uavg;
        (Cells + index_Cell)->V   = Vavg;
        (Cells + index_Cell)->W   = Wavg;
        (Cells + index_Cell)->Rho = rho_ini;

        // REMEMBER BCCONNECTOR FILE
        //  _________________________________________________________________________________________________________________
        // |        0      |       1      |       2      |     3     |    4    |    5     |    6     |    7     |      8     |
        // |  node index i | node index j | node index k | latticeID | BC type | x coord  | y coord  | z coord  | Boundary ID|
        // |_______________|______________|______________|___________|_________|__________|__________|__________|____________|
        //
        // lattice ID is based on the following speed model and it depends on the BC
        //  ID lattice
        //
        //   |       xy - plane          |       xz - plane          |       yz - plane          |
        //   |                           |                           |                           |
        //   |     8       3       7     |     12      5      11     |     16      5       15    |
        //   |       \     |     /       |       \     |     /       |       \     |     /       |
        //   |         \   |   /         |         \   |   /         |         \   |   /         |
        //   |           \ | /           |           \ | /           |           \ | /           |
        //   |     2 - - - 0 - - - 1     |     2 - - - 0 - - - 1     |     4 - - - 0 - - - 3     |
        //   |           / | \           |           / | \           |           / | \           |
        //   |         /   |   \         |         /   |   \         |         /   |   \         |
        //   |       /     |     \       |       /     |     \       |       /     |     \       |
        //   |     10      4       9     |     14      6      13     |     18      6       17    |

        // BC types: *1->wall; *2->inlet; *3->outlet

        (Cells + index_Cell)->BoundaryID = 0;  // IT IS NOT BOUNDARY NODE

        for(j = 0; j < 19; j++)
        {
            (Cells + index_Cell)->BC_ID[j]= 0  ; // IT IS NOT BOUNDARY LATTICE
            (Cells + index_Cell)->Q[j]    = 0.5;
        }

        //SEARCH FOR BC TYPE, BOUNDARY ID AND DISTANCES IN THE LATTICE
        for(j = 0; j < *NumConn; j++)
        {

            //When current cell is on the con iteration
            if (
                    ( (int)Con[j][0] == (int)Nod[i][0] )\
                && ( (int)Con[j][1] == (int)Nod[i][1] )\
                && ( (int)Con[j][2] == (int)Nod[i][2] ) )
            {
                for(k = 1; k < 19; k++)
                {

                    //Checks the direction cell it connects to
                    if ( Con[j][3] == k )
                    {
                        (Cells + index_Cell)->BC_ID[k]   = (int) Con[j][4];         //Gets the connection type(1:static, 2:mobile)
                        (Cells + index_Cell)->BoundaryID = (int) Con[j][8];         //Sets the direction as boundary node

                        // find distance from the boundary
                        (Cells + index_Cell)->Q[k] = sqrt(
                                pow(Con[j][5]- ((Cells + index_Cell)->CoordX),2)\
                                + pow(Con[j][6]-((Cells + index_Cell)->CoordY),2)\
                                + pow(Con[j][7]-((Cells + index_Cell)->CoordZ),2)) / ((*Delta)*Qlat[k]);
                    }
                }
            }
        }


        //////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////
        //////////WITH THE ACTUAL MESH THERE ARE NO CORNERS///////////////
        //////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////

        // CORNER DETERMINATION
        (Cells + index_Cell)->Boundary = 0;  // IT IS NOT BOUNDARY NODE   (0->Fluid / 1-> Wall / 2->Inlet or outlet / 3-> Edge / 4-> Corner

        int number_solid_lattices;
        number_solid_lattices = 0;
        int boundary;
        boundary = 0;

        for (j = 0; j < 19; ++j)
        {

            if ((Cells + index_Cell)->BC_ID[j]!=0)
            {
                boundary = (Cells + index_Cell)->BC_ID[j];
                ++number_solid_lattices;
            }
        }

        switch(number_solid_lattices)
        {
            case 5:
                if (boundary == 1)                    // SOLID PLANE
                {
                    (Cells + index_Cell)->Boundary = 1;
                }
                else if (boundary ==2)                // FLOW PLANE
                {
                    (Cells + index_Cell)->Boundary = 2;
                }
                else
                {
                    printf("Attention! Node %i is not in a wall nor an inlet/outlet!! (See CellFunctions.c, CORNER DETERMINATION)\n", (Cells + index_Cell)->ID);
                }
                break;
            case 9: (Cells + index_Cell)->Boundary = 3;   break; // EDGE
            case 12: (Cells + index_Cell)->Boundary = 4;  break; // CORNER
        }

        /*<--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // (0->Fluid / 1-> Wall / 2->Inlet or outlet / 3-> Edge (Wall-Wall) / 4-> Edge (Wall-Inlet or outlet) / 5-> Corner (Wall-Wall) / 6-> Corner (Wall-Inlet or outlet)                             |

        bool different_types_of_boundaries;
        different_types_of_boundaries = false;
                                                                                                                                                                                                   |
        switch(number_solid_lattices)                                                                                                                                                                  |
        {                                                                                                                                                                                              |
          case 5:                                                                                                                                                                                      |
            if (boundary == 1)                  // SOLID PLANE                                                                                                                                         |
            {                                                                                                                                                                                          |
              (Cells + index_Cell)->Boundary = 1;                                                                                                                                                      |
            }                                                                                                                                                                                          |
            else if (boundary ==2)              // FLOW PLANE                                                                                                                                          |
            {                                                                                                                                                                                          |
              (Cells + index_Cell)->Boundary = 2;                                                                                                                                                      |
            }                                                                                                                                                                                          |
            else                                                                                                                                                                                       |
            {                                                                                                                                                                                          |
              printf("Attention! Node %i is not in a wall nor an inlet/outlet!! (See CellFunctions.c, CORNER DETERMINATION)\n", (Cells + index_Cell)->ID);                                             |
            }                                                                                                                                                                                          |
            break;                                                                                                                                                                                     |
          case 9:                                                                                                                                                                                      |
            if (different_types_of_boundaries)  // SOLID EDGE                                                                                                                                          |
            {                                                                                                                                                                                          |
              (Cells + index_Cell)->Boundary = 3;                                                                                                                                                      |
            }                                                                                                                                                                                          |
            else                                // SOLID - FLOW EDGE                                                                                                                                   |
            {                                                                                                                                                                                          |
              (Cells + index_Cell)->Boundary = 4;                                                                                                                                                      |
            }                                                                                                                                                                                          |
            break;                                                                                                                                                                                     |
          case 12:                                                                                                                                                                                     |
            if (different_types_of_boundaries)  // SOLID CORNERS                                                                                                                                       |
            {                                                                                                                                                                                          |
              (Cells + index_Cell)->Boundary = 5;                                                                                                                                                      |
            }                                                                                                                                                                                          |
            else                                // SOLID - FLOW CORNER                                                                                                                                 |
            {                                                                                                                                                                                          |
              (Cells + index_Cell)->Boundary = 6;                                                                                                                                                      |
            }                                                                                                                                                                                          |
            break;                                                                                                                                                                                     |
          default: printf("Attention! Node %i is not in a wall, inlet/outlet, edge or corner!! (See CellFunctions.c, CORNER DETERMINATION)\n", (Cells + index_Cell)->ID);                              |
        } ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        /*
        for(j = 1; j < 19; j++)
        {
          if ((Cells + index_Cell)->BC_ID[j] != 0)
          { // Lattice-direction j is not surrounded by fluid

            if ((Cells + index_Cell)->Boundary == 0)
            {// The BC of the node becomes the BC of the first lattice direction which doesn't point to fluid
              (Cells + index_Cell)->Boundary = (Cells + index_Cell)->BC_ID[j];
            }

            else
            {// In the same node two lattice directions have different BC (corners). This automatically means that the node is a corner.
              if (((Cells + index_Cell)->Boundary) != ((Cells + index_Cell)->BC_ID[j]) )
              {
                (Cells + index_Cell)->Boundary=1;
                (Cells + index_Cell)->Corner=1;
              }
              /*
              if (((Cells + index_Cell)->Boundary) < ((Cells + index_Cell)->BC_ID[j]) )
              {
                (Cells + index_Cell)->Boundary=1;
                (Cells + index_Cell)->Corner=1;
              }
              if ((Cells + index_Cell)->Boundary > (Cells + index_Cell)->BC_ID[j])
              {
                (Cells + index_Cell)->Boundary=1;
                (Cells + index_Cell)->Corner=1;
              }
            }
          }*/


        /*

        //////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////
        //////////WITH THE ACTUAL MESH THERE ARE NO CORNERS///////////////
        //////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////

        WITH THE NEW APPROACH, CORNERS ARE TREATED AS WALLS

        //Boundary=0 NO BC
        //Boundary=1 WALL
        //Boundary=2 INLET
        //Boundary=3 OUTLET


        // IF TWO PERPENDICULAR LATTICE DIRECTIONS ARE NOT SURROUNDED BY FLUID IT MEANS THAT IT IS
        // A CORNER. IT IS MANDATORY THAT THE DIAGONAL LATTICE DIRECTION BETWEEN THEM CORRESPONDS
        // TO A WALL. HENCE, ALL THE CASES ARE CHECKED AND DIAGONAL DIRECTIONS EQUALED TO 1
        if ( ((Cells + index_Cell)->Corner) == 1)
        {
          if ((Cells + index_Cell)->BC_ID[1] != 0 && (Cells + index_Cell)->BC_ID[3] != 0 )
          {
            (Cells + index_Cell)->BC_ID[7]=1;
          }
          if ((Cells + index_Cell)->BC_ID[1] != 0 && (Cells + index_Cell)->BC_ID[4] != 0 )
          {
            (Cells + index_Cell)->BC_ID[9]=1;
          }
          if ((Cells + index_Cell)->BC_ID[1] != 0 && (Cells + index_Cell)->BC_ID[5] != 0 )
          {
            (Cells + index_Cell)->BC_ID[11]=1;
          }
          if ((Cells + index_Cell)->BC_ID[1] != 0 && (Cells + index_Cell)->BC_ID[6] != 0 )
          {
            (Cells + index_Cell)->BC_ID[13]=1;
          }
          if ((Cells + index_Cell)->BC_ID[2] != 0 && (Cells + index_Cell)->BC_ID[4] != 0 )
          {
            (Cells + index_Cell)->BC_ID[10]=1;
          }
          if ((Cells + index_Cell)->BC_ID[2] != 0 && (Cells + index_Cell)->BC_ID[3] != 0 )
          {
            (Cells + index_Cell)->BC_ID[8]=1;
          }
          if ((Cells + index_Cell)->BC_ID[2] != 0 && (Cells + index_Cell)->BC_ID[5] != 0 )
          {
            (Cells + index_Cell)->BC_ID[12]=1;
          }
          if ((Cells + index_Cell)->BC_ID[2] != 0 && (Cells + index_Cell)->BC_ID[6] != 0 )
          {
            (Cells + index_Cell)->BC_ID[14]=1;
          }
          if ((Cells + index_Cell)->BC_ID[3] != 0 && (Cells + index_Cell)->BC_ID[5] != 0 )
          {
            (Cells + index_Cell)->BC_ID[15]=1;
          }
          if ((Cells + index_Cell)->BC_ID[3] != 0 && (Cells + index_Cell)->BC_ID[6] != 0 )
          {
            (Cells + index_Cell)->BC_ID[17]=1;
          }
          if ((Cells + index_Cell)->BC_ID[4] != 0 && (Cells + index_Cell)->BC_ID[5] != 0 )
          {
            (Cells + index_Cell)->BC_ID[16]=1;
          }
          if ((Cells + index_Cell)->BC_ID[4] != 0 && (Cells + index_Cell)->BC_ID[6] != 0 )
          {
            (Cells + index_Cell)->BC_ID[18]=1;
          }
        }
        */

        // INITIALIZE STREAMING (STREAM EVERYWHERE)
        for(j = 0; j < 19; j++)
        {
            (Cells + index_Cell)->StreamLattice[j] = 1;
        }


        // DON'T STREAM TO SOLIDS
        for(j = 0; j < 19; j++)
        {

            if ((Cells + index_Cell)->BC_ID[j]!=0)
            {
                (Cells + index_Cell)->StreamLattice[opp[j]]= 0 ;

                /*if (MYTHREAD == 0 && i == 30) {
                    printf("opp[%i] = %i\n", j, opp[j]);
                    printf("Str(%i): %i\n",opp[j],(Cells + index_Cell)->StreamLattice[opp[j]]);
                }*/
            }

            //x de 0 a 10
            //j = 2
            //z 0 o 1
            /*if(MYTHREAD == 0 && i == 30) {
                printf("i: %i,  ", i);
                printf("BC_ID(%i): %i, ",j,(Cells + index_Cell)->BC_ID[j]);
                printf("Str(%i): %i ",j,(Cells + index_Cell)->StreamLattice[j]);
                printf("Str(%i): %i\n",opp[j],(Cells + index_Cell)->StreamLattice[opp[j]]);
            }*/
        }


        // INLET VELOCITY // THIS IS CRAPPY, NOT USED!
        switch(InletProfile)
        {
            case 1:
                ////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////
                ////////////////////////TO DEVELOPE/////////////////////////////////
                ////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////
                //Uo=1.5*Uavg*(1-4*(pow((CoordY-MinInletCoordY-0.5*(MaxInletCoordY-MinInletCoordY)),2)));
                //Uo=4*1.5*Uavg*CoordY*(41-CoordY)/(41*41);
                //(Cells + index_Cell)->Uo = 4*1.5*Uavg*(((Cells + index_Cell)->CoordY)-(*MinInletCoordY))*(((*MaxInletCoordY)-(*MinInletCoordY))-(((Cells + index_Cell)->CoordY)-(*MinInletCoordY)))/(((*MaxInletCoordY)-(*MinInletCoordY))*((*MaxInletCoordY)-(*MinInletCoordY)));
                //(Cells + index_Cell)->Vo = Vavg;
                //(Cells + index_Cell)->Wo = ?????????
                break;
            case 2:
                (Cells + index_Cell)->Uo = Uavg;
                (Cells + index_Cell)->Vo = Vavg;
                (Cells + index_Cell)->Wo = Wavg;
                break;

        }


        //(Cells + index_Cell)->U = (Cells + index_Cell)->Uo;
        //(Cells + index_Cell)->V = (Cells + index_Cell)->Vo;
        //(Cells + index_Cell)->W = (Cells + index_Cell)->Wo;
        if((Cells + index_Cell)->CoordZ==1) {

            (Cells + index_Cell)->U = (Cells + index_Cell)->Uo;
            (Cells + index_Cell)->V = (Cells + index_Cell)->Vo;
            (Cells + index_Cell)->W = (Cells + index_Cell)->Wo;
        }
        else {
            (Cells + index_Cell)->U = 0;
            (Cells + index_Cell)->V = 0;
            (Cells + index_Cell)->W = 0;
        }

    } // END OF for LOOP
} // END OF FUNCTION


void CellIni_NEW(CellProps *Cells,
                 float  **Nod,
                 float  **Con,
                 float  Uavg,
                 float  Vavg,
                 float  Wavg,
                 int    InletProfile,
                 int    CollisionModel,
                 int*   opp,
                 float  rho_ini)
{
    ///////////// DECLARATION /////////////////
    int l, n;  // loop variables
    int i_r, j_r, k_r;  // interior counting variables
    float Qlat[19];  // Q lattice, length of the discrete directions
    Qlat[0]  =    0   ;
    Qlat[1]  =    1   ;
    Qlat[2]  =    1   ;
    Qlat[3]  =    1   ;
    Qlat[4]  =    1   ;
    Qlat[5]  =    1   ;
    Qlat[6]  =    1   ;
    Qlat[7]  = (float)sqrt(2);
    Qlat[8]  = (float)sqrt(2);
    Qlat[9]  = (float)sqrt(2);
    Qlat[10] = (float)sqrt(2);
    Qlat[11] = (float)sqrt(2);
    Qlat[12] = (float)sqrt(2);
    Qlat[13] = (float)sqrt(2);
    Qlat[14] = (float)sqrt(2);
    Qlat[15] = (float)sqrt(2);
    Qlat[16] = (float)sqrt(2);
    Qlat[17] = (float)sqrt(2);
    Qlat[18] = (float)sqrt(2);

    int index_Cell;         //Local global index inside a block
    int lID;          //Local real index inside a block
    for (int k = 1; k < LAT+1; k++) {
        for (int j = 1; j< LAT+1; j++) {
            for (int i = 1; i< LAT+1; i++) {
                i_r = i-1;
                j_r = j-1;
                k_r = k-1;


                index_Cell = i + j*(LAT+2) + k*(LAT+2)*(LAT+2);
                lID = i_r + j_r*(LAT) + k_r*(LAT)*(LAT);

                //local move
                int l_i = lID%LAT;
                int l_j = ((lID/LAT)%LAT)*NN;
                int l_k = (lID/(LAT*LAT))*NN*NM;

                //global move
                int g_i = (MYTHREAD%NTDX)*LAT;
                int g_j = ((MYTHREAD/NTDX)%NTDY)*NTDX*LAT*LAT;
                int g_k = (MYTHREAD/(NTDX*NTDY))*NTDX*NTDY*LAT*LAT*LAT;


                int ID  = l_i + l_j + l_k + g_i + g_j + g_k;

                /*if (MYTHREAD == 0) {
                    printf("T=%i: (%i,%i,%i), (%i,%i,%i), index_Cell = %i, lID = %i\n",
                           MYTHREAD,i, j, k, i_r, j_r, k_r, index_Cell, lID);
                    printf("T=%i: %i + %i + %i + %i + %i + %i = %i\n",
                           MYTHREAD, l_i, l_j, l_k, g_i, g_j, g_k, ID);
                }*/
                upc_barrier;
                //index_Cell = LAYER - MYTHREAD * BLOCKSIZE + i;

                // FIND ID of the actual cell
                (Cells + index_Cell)->ID = ID;



                // FIND X, Y and Z of the actual cell
                (Cells + index_Cell)->CoordX = Nod[ID][3];
                (Cells + index_Cell)->CoordY = Nod[ID][4];
                (Cells + index_Cell)->CoordZ = Nod[ID][5];
                if (MYTHREAD==0 && i_r==0 && j_r==0 && k_r ==0) {
                    printf("Nod[%i] = (%f,%f,%f)\n",ID,Nod[ID][3],Nod[ID][4],Nod[ID][5]);
                }

                // CHECK FLUID OR NOT
                //(Cells + index_Cell)->Fluid  = (int) Nod[i][6];

                // Which thread does it belongs to
                (Cells + index_Cell)->ThreadNumber = MYTHREAD;

                // INITIALIZE VARIABLEs
                (Cells + index_Cell)->U   = Uavg;
                (Cells + index_Cell)->V   = Vavg;
                (Cells + index_Cell)->W   = Wavg;
                (Cells + index_Cell)->Rho = rho_ini;

                // REMEMBER BCCONNECTOR FILE
                //  _________________________________________________________________________________________________________________
                // |        0      |       1      |       2      |     3     |    4    |    5     |    6     |    7     |      8     |
                // |  node index i | node index j | node index k | latticeID | BC type | x coord  | y coord  | z coord  | Boundary ID|
                // |_______________|______________|______________|___________|_________|__________|__________|__________|____________|
                //
                // lattice ID is based on the following speed model and it depends on the BC
                //  ID lattice
                //
                //   |       xy - plane          |       xz - plane          |       yz - plane          |
                //   |                           |                           |                           |
                //   |     8       3       7     |     12      5      11     |     16      5       15    |
                //   |       \     |     /       |       \     |     /       |       \     |     /       |
                //   |         \   |   /         |         \   |   /         |         \   |   /         |
                //   |           \ | /           |           \ | /           |           \ | /           |
                //   |     2 - - - 0 - - - 1     |     2 - - - 0 - - - 1     |     4 - - - 0 - - - 3     |
                //   |           / | \           |           / | \           |           / | \           |
                //   |         /   |   \         |         /   |   \         |         /   |   \         |
                //   |       /     |     \       |       /     |     \       |       /     |     \       |
                //   |     10      4       9     |     14      6      13     |     18      6       17    |

                // BC types: *1->wall; *2->inlet; *3->outlet

                (Cells + index_Cell)->BoundaryID = 0;  // IT IS NOT BOUNDARY NODE

                for(l = 0; l < 19; l++)
                {
                    (Cells + index_Cell)->BC_ID[l]= 0  ; // IT IS NOT BOUNDARY LATTICE
                    (Cells + index_Cell)->Q[l]    = 0.5;
                }

                //SEARCH FOR BC TYPE, BOUNDARY ID AND DISTANCES IN THE LATTICE
                for(l = 0; l < *NumConn; l++)
                {

                    //When current cell is on the con iteration
                    if (
                            ( (int)Con[l][0] == (int)Nod[ID][0] )\
                            && ( (int)Con[l][1] == (int)Nod[ID][1] )\
                            && ( (int)Con[l][2] == (int)Nod[ID][2] ) )
                    {
                        for(n = 1; n < 19; n++)
                        {

                            //Checks the direction cell it connects to
                            if ( Con[l][3] == n )
                            {
                                (Cells + index_Cell)->BC_ID[n]   = (int) Con[l][4];         //Gets the connection type(1:static, 2:mobile)
                                (Cells + index_Cell)->BoundaryID = (int) Con[l][8];         //Sets the direction as boundary node

                                // find distance from the boundary
                                (Cells + index_Cell)->Q[n] = sqrt(
                                        pow(Con[l][5]- ((Cells + index_Cell)->CoordX),2)\
                                + pow(Con[l][6]-((Cells + index_Cell)->CoordY),2)\
                                + pow(Con[l][7]-((Cells + index_Cell)->CoordZ),2)) / ((*Delta)*Qlat[n]);
                            }
                        }
                    }
                }


                //////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////
                //////////WITH THE ACTUAL MESH THERE ARE NO CORNERS///////////////
                //////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////

                // CORNER DETERMINATION
                (Cells + index_Cell)->Boundary = 0;  // IT IS NOT BOUNDARY NODE   (0->Fluid / 1-> Wall / 2->Inlet or outlet / 3-> Edge / 4-> Corner

                int number_solid_lattices;
                number_solid_lattices = 0;
                int boundary;
                boundary = 0;

                for (l = 0; l < 19; ++l)
                {

                    if ((Cells + index_Cell)->BC_ID[l]!=0)
                    {
                        boundary = (Cells + index_Cell)->BC_ID[l];
                        ++number_solid_lattices;
                    }
                }

                switch(number_solid_lattices)
                {
                    case 5:
                        if (boundary == 1)                    // SOLID PLANE
                        {
                            (Cells + index_Cell)->Boundary = 1;
                        }
                        else if (boundary ==2)                // FLOW PLANE
                        {
                            (Cells + index_Cell)->Boundary = 2;
                        }
                        else
                        {
                            printf("Attention! Node %i is not in a wall nor an inlet/outlet!! (See CellFunctions.c, CORNER DETERMINATION)\n", (Cells + index_Cell)->ID);
                        }
                        break;
                    case 9: (Cells + index_Cell)->Boundary = 3;   break; // EDGE
                    case 12: (Cells + index_Cell)->Boundary = 4;  break; // CORNER
                }

                // INITIALIZE STREAMING (STREAM EVERYWHERE)
                for(l = 0; l < 19; l++)
                {
                    (Cells + index_Cell)->StreamLattice[l] = 1;
                }


                // DON'T STREAM TO SOLIDS
                for(l = 0; l < 19; l++)
                {

                    if ((Cells + index_Cell)->BC_ID[l]!=0)
                    {
                        (Cells + index_Cell)->StreamLattice[opp[l]]= 0 ;

                    }
                }


                // INLET VELOCITY // THIS IS CRAPPY, NOT USED!
                switch(InletProfile)
                {
                    case 1:
                        ////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////
                        ////////////////////////TO DEVELOPE/////////////////////////////////
                        ////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////
                        //Uo=1.5*Uavg*(1-4*(pow((CoordY-MinInletCoordY-0.5*(MaxInletCoordY-MinInletCoordY)),2)));
                        //Uo=4*1.5*Uavg*CoordY*(41-CoordY)/(41*41);
                        //(Cells + index_Cell)->Uo = 4*1.5*Uavg*(((Cells + index_Cell)->CoordY)-(*MinInletCoordY))*(((*MaxInletCoordY)-(*MinInletCoordY))-(((Cells + index_Cell)->CoordY)-(*MinInletCoordY)))/(((*MaxInletCoordY)-(*MinInletCoordY))*((*MaxInletCoordY)-(*MinInletCoordY)));
                        //(Cells + index_Cell)->Vo = Vavg;
                        //(Cells + index_Cell)->Wo = ?????????
                        break;
                    case 2:
                        (Cells + index_Cell)->Uo = Uavg;
                        (Cells + index_Cell)->Vo = Vavg;
                        (Cells + index_Cell)->Wo = Wavg;
                        break;

                }
                //For lid driven initialize uper nodes
                if((Cells + index_Cell)->CoordZ==1) {

                    (Cells + index_Cell)->U = (Cells + index_Cell)->Uo;
                    (Cells + index_Cell)->V = (Cells + index_Cell)->Vo;
                    (Cells + index_Cell)->W = (Cells + index_Cell)->Wo;
                }
                else {
                    (Cells + index_Cell)->U = 0;
                    (Cells + index_Cell)->V = 0;
                    (Cells + index_Cell)->W = 0;
                }


            }
        }
    } // END OF for LOOP
} // END OF FUNCTION

/*==================================================
========Creating constant lattice parameters========
==================================================*/

//void D3Q19Vars(double* w, int* cx, int* cy, int* cz, int* opp, int* c)
void D3Q19Vars()
{
    // Fill up variables with constants
    //  ID lattice
    //
    //  |       xy - plane          |       xz - plane          |       yz - plane          |
    //  |                           |                           |                           |
    //  |     8       3       7     |     12      5      11     |     16      5       15    |
    //  |       \     |     /       |       \     |     /       |       \     |     /       |
    //  |         \   |   /         |         \   |   /         |         \   |   /         |
    //  |           \ | /           |           \ | /           |           \ | /           |
    //  |     2 - - - 0 - - - 1     |     2 - - - 0 - - - 1     |     4 - - - 0 - - - 3     |
    //  |           / | \           |           / | \           |           / | \           |
    //  |         /   |   \         |         /   |   \         |         /   |   \         |
    //  |       /     |     \       |       /     |     \       |       /     |     \       |
    //  |     10      4       9     |     14      6      13     |     18      6       17    |

    int i;

    // D3Q19 properties
    w[0]=12./36.;

    for (i=1; i<7; i++ )
        w[i]=2./36.;

    for (i=7; i<19; i++ )
        w[i]=1./36.;

    cx[0]  =  0;
    cx[1]  =  1;
    cx[2]  = -1;
    cx[3]  =  0;
    cx[4]  =  0;
    cx[5]  =  0;
    cx[6]  =  0;
    cx[7]  =  1;
    cx[8]  = -1;
    cx[9]  =  1;
    cx[10] = -1;
    cx[11] =  1;
    cx[12] = -1;
    cx[13] =  1;
    cx[14] = -1;
    cx[15] =  0;
    cx[16] =  0;
    cx[17] =  0;
    cx[18] =  0;

    cy[0]  =  0;
    cy[1]  =  0;
    cy[2]  =  0;
    cy[3]  =  1;
    cy[4]  = -1;
    cy[5]  =  0;
    cy[6]  =  0;
    cy[7]  =  1;
    cy[8]  =  1;
    cy[9]  = -1;
    cy[10] = -1;
    cy[11] =  0;
    cy[12] =  0;
    cy[13] =  0;
    cy[14] =  0;
    cy[15] =  1;
    cy[16] = -1;
    cy[17] =  1;
    cy[18] = -1;

    cz[0]  =  0;
    cz[1]  =  0;
    cz[2]  =  0;
    cz[3]  =  0;
    cz[4]  =  0;
    cz[5]  =  1;
    cz[6]  = -1;
    cz[7]  =  0;
    cz[8]  =  0;
    cz[9]  =  0;
    cz[10] =  0;
    cz[11] =  1;
    cz[12] =  1;
    cz[13] = -1;
    cz[14] = -1;
    cz[15] =  1;
    cz[16] =  1;
    cz[17] = -1;
    cz[18] = -1;


    opp[0]  =  0;
    opp[1]  =  2;
    opp[2]  =  1;
    opp[3]  =  4;
    opp[4]  =  3;
    opp[5]  =  6;
    opp[6]  =  5;
    opp[7]  = 10;
    opp[8]  =  9;
    opp[9]  =  8;
    opp[10] =  7;
    opp[11] = 14;
    opp[12] = 13;
    opp[13] = 12;
    opp[14] = 11;
    opp[15] = 18;
    opp[16] = 17;
    opp[17] = 16;
    opp[18] = 15;

    // n: number of cols (x-direction)
    // m: number of rows (y-direction)
    // LAYER = n*m

    // Streaming comes from which node?
    int moveY = (LAT+2);
    int moveZ = (LAT+2)*(LAT+2);

    c[0]  =          0              ;
    c[1]  = -1                      ; // (i-1)
    c[2]  = +1                      ; // (i+1)
    c[3]  =     -1*moveY            ; //         (j-1)
    c[4]  =     +1*moveY            ; //         (j+1)
    c[5]  =                 -1*moveZ; //                 (k-1)
    c[6]  =                 +1*moveZ; //                 (k+1)
    c[7]  = -1  -1*moveY            ; // (i-1)   (j-1)
    c[8]  = +1  -1*moveY            ; // (i+1)   (j-1)
    c[9]  = -1  +1*moveY            ; // (i-1)   (j+1)
    c[10] = +1  +1*moveY            ; // (i+1)   (j+1)
    c[11] = -1              -1*moveZ; // (i-1)           (k-1)
    c[12] = +1              -1*moveZ; // (i+1)           (k-1)
    c[13] = -1              +1*moveZ; // (i-1)           (k+1)
    c[14] = +1              +1*moveZ; // (i+1)           (k+1)
    c[15] =     -1*moveY    -1*moveZ; //         (j-1)   (k-1)
    c[16] =     +1*moveY    -1*moveZ; //         (j+1)   (k-1)
    c[17] =     -1*moveY    +1*moveZ; //         (j-1)   (k+1)
    c[18] =     +1*moveY    +1*moveZ; //         (j+1)   (k+1)



    norm[0][0] = -1;
    norm[0][1] = 0;
    norm[0][2] = 0;

    norm[1][0] = 1;
    norm[1][1] = 0;
    norm[1][2] = 0;


    norm[2][0] = 0;
    norm[2][1] = -1;
    norm[2][2] = 0;

    norm[3][0] = 0;
    norm[3][1] = 1;
    norm[3][2] = 0;


    norm[4][0] = 0;
    norm[4][1] = 0;
    norm[4][2] = -1;

    norm[5][0] = 0;
    norm[5][1] = 0;
    norm[5][2] = 1;


    j_wall_unknown[0][0] = 1;
    j_wall_unknown[0][1] = 7;
    j_wall_unknown[0][2] = 9;
    j_wall_unknown[0][3] = 11;
    j_wall_unknown[0][4] = 13;

    j_wall_unknown[1][0] = 2;
    j_wall_unknown[1][1] = 8;
    j_wall_unknown[1][2] = 10;
    j_wall_unknown[1][3] = 12;
    j_wall_unknown[1][4] = 14;

    j_wall_unknown[2][0] = 3;
    j_wall_unknown[2][1] = 7;
    j_wall_unknown[2][2] = 8;
    j_wall_unknown[2][3] = 15;
    j_wall_unknown[2][4] = 17;

    j_wall_unknown[3][0] = 4;
    j_wall_unknown[3][1] = 9;
    j_wall_unknown[3][2] = 10;
    j_wall_unknown[3][3] = 16;
    j_wall_unknown[3][4] = 18;

    j_wall_unknown[4][0] = 5;
    j_wall_unknown[4][1] = 11;
    j_wall_unknown[4][2] = 12;
    j_wall_unknown[4][3] = 15;
    j_wall_unknown[4][4] = 16;

    j_wall_unknown[5][0] = 6;
    j_wall_unknown[5][1] = 13;
    j_wall_unknown[5][2] = 14;
    j_wall_unknown[5][3] = 17;
    j_wall_unknown[5][4] = 18;


}


// End of function

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  BELOW HERE ////////////////////////////////////////
/////////////////////////////// IS STILL IN ////////////////////////////////////////
///////////////////////////////      2D     ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


/*==================================================
=============Function for the BGKW model============
==================================================*/

void BGKW(int i, double Omega)
{

    double T1 = ((Cells+i)->U) * ((Cells+i)->U) + ((Cells+i)->V) * ((Cells+i)->V) + ((Cells+i)->W) * ((Cells+i)->W);
    double T2 = 0.0;
    int k;
    for (k = 0; k < 19; k++)
    {
        T2                  = ((Cells+i)->U) * cx[k] + ((Cells+i)->V)*cy[k] + ((Cells+i)->W)*cz[k];
        (Cells+i)->Feq[k]   = ((Cells+i)->Rho) * w[k] * (1.0 + 3.0*T2 + 4.5*T2*T2 - 1.5*T1);
        (Cells+i)->METAF[k] = Omega*((Cells+i)->Feq[k])+(1.0-Omega)*((Cells+i)->F[k]);
    }
}

/*==================================================
=============Function for the TRT model=============
==================================================*/

void TRT(CellProps *Cells, int i, double* w, int* cx, int* cy, int* opp, double Omega, double OmegaA)
{
    double T1 = ((Cells +i)->U)*((Cells +i)->U)+((Cells +i)->V)*((Cells +i)->V);
    double T2 = 0.0;
    int k;
    for (k=0; k<9; k++)
    {
        T2                  = ((Cells +i)->U)   *cx[k] + ((Cells +i)->V)*cy[k];
        (Cells +i)->Feq[k]  = ((Cells +i)->Rho) *w[k]*(1.0+3.0*T2+4.5*T2*T2-1.5*T1);
    }

    float F_p[9],Feq_p[9],F_m[9],Feq_m[9];
    for (k=0; k<9; k++)
    {
        F_p[k]   = 0.5*( (Cells +i)->F[k]   + (Cells +i)->F[opp[k]]   );
        Feq_p[k] = 0.5*( (Cells +i)->Feq[k] + (Cells +i)->Feq[opp[k]] );
        F_m[k]   = 0.5*( (Cells +i)->F[k]   - (Cells +i)->F[opp[k]]   );
        Feq_m[k] = 0.5*( (Cells +i)->Feq[k] - (Cells +i)->Feq[opp[k]] );

        (Cells +i)->METAF[k] = (Cells +i)->F[k] -
                               (F_p[k]-Feq_p[k])*Omega -
                               (F_m[k]-Feq_m[k])*OmegaA;
    }
}


/*==================================================
=============Function for the MRT model=============
==================================================*/

void MRT(CellProps *Cells, int i, double** tm, double** stmiv)
{
    int k, l;
    double fmom[9],fmeq[9];
    double U, V;
    double suma,sumb;

    U = (Cells +i)->U;
    V = (Cells +i)->V;

    fmeq[0] =  (Cells +i)->Rho;
    fmeq[1] = ((Cells +i)->Rho)*(-2.0+3.0*((Cells +i)->Rho)*(U*U+V*V));
    fmeq[2] = ((Cells +i)->Rho)*(1.0-3.0*((Cells +i)->Rho)*(U*U+V*V));
    fmeq[3] = ((Cells +i)->Rho)*((Cells +i)->U);
    fmeq[4] =-((Cells +i)->Rho)*((Cells +i)->U);
    fmeq[5] = ((Cells +i)->Rho)*((Cells +i)->V);
    fmeq[6] =-((Cells +i)->Rho)*((Cells +i)->V);
    fmeq[7] = ((Cells +i)->Rho)*(U*U-V*V);
    fmeq[8] = ((Cells +i)->Rho)*U*V;

    for (k=0; k<9;k++)
    {
        suma=0.0;
        for (l=0; l<9;l++)
            suma = suma + tm[k][l]*((Cells +i)->F[l]);

        fmom[k]=suma;
    }

    for (k=0; k<9;k++)
    {
        sumb=0.0;
        for (l=0; l<9;l++)
            sumb = sumb + stmiv[k][l]*(fmom[l]-fmeq[l]);

        (Cells +i)->METAF[k] = ((Cells +i)->F[k]) - sumb;
    }
}


/*==================================================
======Function to update the distribution fct.======
==================================================*/

void UpdateF(CellProps *Cells, int i)
{
    int k;
    for(k = 0; k < 19; k++)
        (Cells+i)->F[k] = (Cells+i)->METAF[k];
}


/*==================================================
======Function to update the macroscopic var.=======
==================================================*/

void UpdateMacroscopic(CellProps *Cells, int i, int CalculateDragLift)
{
    double Ssum, Usum, Vsum, Wsum;
    int k;

    //if ((Cells+i)->Fluid==1)
    //{
    Ssum=0.0;
    for (k=0; k<19; k++)
        Ssum = Ssum+(Cells+i)->F[k];

    (Cells+i)->Rho = Ssum;

    //printf("Rho[%d][%d] = %f\n", j, i, Ssum);

    Usum = 0.0;
    Vsum = 0.0;
    Wsum = 0.0;
    for (k=0; k<19; k++)
    {
        Usum = Usum + ((Cells+i)->F[k])*cx[k];
        Vsum = Vsum + ((Cells+i)->F[k])*cy[k];
        Wsum = Wsum + ((Cells+i)->F[k])*cz[k];
    }
    (Cells+i)->U = Usum/((Cells+i)->Rho);
    (Cells+i)->V = Vsum/((Cells+i)->Rho);
    (Cells+i)->W = Wsum/((Cells+i)->Rho);
    //}

    /*if ((Cells+i)->BC_ID[1]==3) // for outlet on the right
    {
        (Cells+i)->V=0.0;
    }*/

    //   DRAG/LIFT FORCE
    if (CalculateDragLift != 0 && (Cells+i)->BoundaryID==CalculateDragLift)
    {
        (Cells+i)->DragXF = ((Cells+i)->Rho)/3*(20-(Cells+i)->CoordX)/5;
        (Cells+i)->DragYF = ((Cells+i)->Rho)/3*(20-(Cells+i)->CoordY)/5;
        (Cells+i)->LiftF  = ((Cells+i)->Rho)/3*(20-(Cells+i)->CoordZ)/5;
    }

}

int getAx_F(int face) {
    return face/2;
}
int getDir_F(int face) {
    return face%2;
}
int getF(int Ax, int dir) {
    return Ax*2 + dir;
}

int getAx_E(int edge) {
    return edge/4;
}
int getPos_E(int edge) {
    return edge%4;
}

int getE(int Ax, int pos) {
    return Ax*4 + pos;
}

int getZ_C(int corner) {
    return corner/4;
}
int getY_C(int corner) {
    return (corner - getZ_C(corner) * 4)/2;
}
int getX_C(int corner) {
    return corner - getY_C(corner)*2 - getZ_C(corner)*4;
}
int getC(int X, int Y, int Z) {
    return X + Y*2 + Z*4;
}

int getCubeID(int x, int y, int z) {
    return x + y*NTDX + z*NTDX*NTDY;
}

void getCubeCoords(int ID, int X[]) {
    X[2] = ID/(NTDX*NTDY);
    X[1] = (ID - X[2] * NTDX*NTDY) / NTDX;
    //printf("(%i - %i * %i * %i)/%i = %i\n",ID, X[2],NTDX,NTDY,NTDX,3/NTDX);
    //printf("3/%i = %i\n",NTDX,100/NTDX);
    X[0] = ID - X[1]*NTDX - X[2] * NTDX*NTDY;
}



int getThread(int index) {
    return index/BLOCKSIZE;
}