//
// Created by epolepole on 08/07/16.
//

#ifndef THESIS_UPC_TESTS_H
#define THESIS_UPC_TESTS_H

#include "ShellFunctions.h"
#include <stdio.h>

char tests_file_name [50];
FILE* tests_file;


void init_tests();
void test_all(CellProps* Cells,int* iter);
void end_tests();
void test_stream(int ID, CellProps* Cells, int iter);
void test_collision(CellProps* Cells, int iter);







#endif //THESIS_UPC_TESTS_H
