#ifndef BLOCKSIZEDEFINER_H
#define BLOCKSIZEDEFINER_H

#include <math.h>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))



#ifndef NN
#define NN 30
#endif  //NN

#ifndef NM
#define NM 30
#endif  //NM

#ifndef NL
#define NL 30
#endif  //NL



#define LAYER NN*NM
#define NODES NN*NM*NL
#define LAYERS_PER_THREAD ((int)(NL/THREADS))
#define BLOCKSIZE LAYER*LAYERS_PER_THREAD

//Number of threads per dimension
#define NTD (int) pow(THREADS,1/3.)

//Lateral of cube
#define LAT (int) pow(NODES/THREADS,1/3.)

#define NTDX (int)NN/LAT
#define NTDY (int)NM/LAT
#define NTDZ (int)NL/LAT


#define BLOCKSIZE_NEW (int) pow(LAT,3.)
        //           Faces          Edges    Corners
#define B_CELLS_SIZE (int)(6*pow(LAT,2) + 12*LAT + 8)


#endif
