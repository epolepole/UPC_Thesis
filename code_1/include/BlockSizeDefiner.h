#ifndef BLOCKSIZEDEFINER_H
#define BLOCKSIZEDEFINER_H


#define NN 30
#define NM 30
#define NL 30
#define LAYER NN*NM
#define NODES NN*NM*NL
#define BLOCKSIZE NN*((int)(LAYER/THREADS))
#define NUMBER_OF_LAYERS ((int)(NL/THREADS))

#endif
