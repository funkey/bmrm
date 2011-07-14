#include <stdlib.h>
#include <stdio.h>
#ifndef LAP_HPP
#define LAP_HPP

typedef double cost;

#define INF 99999999.9

void lap(int col_size, int row_size, cost **C, int *col, int *row);
#endif
