//
// Created by epolepole on 08/07/16.
//

#include <math.h>
#include <CellFunctions.h>

#include "tests.h"

void init_tests() {
    sprintf(tests_file_name,"Results/tests.txt");
}

void test_all(CellProps* Cells,int* iter) {
    test_stream(LAYER +(15 + 15*30 + 10*LAYER),Cells,*iter);
    test_collision(Cells,*iter);
}
void end_tests() {
}

void test_stream(int ID,CellProps* Cells,int iter) {
    double E, L2;
    double temp = 0;
    E = 0;
    L2 = 0;

    if (MYTHREAD == 0) {
        printf(tests_file);
        printf("\n");
        tests_file = fopen(tests_file_name, "w");
        fprintf(tests_file, "****Streaming Test****\n\n");
        fprintf(tests_file, "Iteration %i\n", iter);
        fprintf(tests_file, "    ID = %i:\n", ID);
        fclose(tests_file);
        printf("****Streaming Test****\n\n");
        printf("Iteration %i\n", iter);
        printf("    ID = %i:\n", ID);
    }
    upc_barrier;

    for(int T = 0; T<THREADS;T++) {
        if(MYTHREAD == T) {
            tests_file = fopen(tests_file_name, "a");
            fprintf(tests_file,"        Thread %i\n",T);
            printf("        Thread %i\n",T);
            for (int j = 0; j<19; j++) {
                temp = pow(Cells[LAYER + ID +c[j]].METAF[j]-Cells[LAYER + ID].F[j],2);
                E = E+ temp;
                fprintf(tests_file,"        Error[%i] = %f\n", j, temp);
                printf("        Error[%i] = %f\n", j, temp);
            }
            fclose(tests_file);
        }
        upc_barrier;
    }
}

void test_collision(CellProps* Cells,int iter) {

    int ID;
    double temp;
    if (MYTHREAD == 0) {
        tests_file = fopen(tests_file_name, "a");
        fprintf(tests_file, "****Collision Test****\n\n");
        fprintf(tests_file, "Iteration %i\n", iter);
        fclose(tests_file);
    }

    ID = (15 + 0 * NN + 10 * LAYER);
    if(MYTHREAD == 0) {
        tests_file = fopen(tests_file_name, "a");
        fprintf(tests_file, "    ID = %i:\n", ID);
        fclose(tests_file);
    }

    for(int T = 0; T<THREADS;T++) {
        if(MYTHREAD == T) {
            tests_file = fopen(tests_file_name, "a");
            fprintf(tests_file,"        Thread %i\n",T);
            int j [5] = {4,9,10,16,18};
            for (int k = 0; k<5;k ++) {
                temp = pow(Cells[LAYER + ID].F[j[k]]-Cells[LAYER + ID].F[opp[k]],2);
                fprintf(tests_file,"        Error[%i] = %f\n", j[k], temp);
            }
            fclose(tests_file);
        }
        upc_barrier;
    }


}