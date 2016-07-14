#ifndef BLOCKSIZEDEFINER_H
#define BLOCKSIZEDEFINER_H


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

#endif
