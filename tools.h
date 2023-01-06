#ifndef __TOOLS_H__
#define __TOOLS_H__

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <sys/time.h>
#include <limits.h>
#include "mpi.h"
#include "stdarg.h"

// random real in [0,1]
#define RAND01 ((double)random() / RAND_MAX)

// A position of a cell in the grid.
typedef struct
{
    int x, y;
} position;

// A grid.
typedef struct
{
    int X, Y;       // dimensions: X and Y
    int **value;    // cell values: value[i][j], 0<=i<X, 0<=j<Y
    int **mark;     // cell markings: mark[i][j], 0<=i<X, 0<=j<Y
    position start; // position of the source
    position end;   // position of the destination
} grid;

// Possible values for the cells of a grid for the .value and .mark fields.
// The order is important: it must be consistent with the color[] arrays (from tools.c)
// and weight[] (from a_star.c).
enum
{
    // for .value
    V_FREE = 0, // empty cell
    V_WALL,     // Wall
    V_SAND,     // Sand
    V_WATER,    // Water
    V_MUD,      // Mud
    V_GRASS,    // Grass
    V_TUNNEL,   // Tunnel

    // for .mark
    M_NULL,  // unmarked vertex
    M_USED,  // vertex marked in P
    M_FRONT, // vertex marked in Q
    M_PATH,  // vertex in the path
};

// Drawing and grid construction routines. The (0,0) point of the grid is the top left corner.
// For more details on the functions, see tools.c

void saveGridValueFile(grid G, char *filename); // saves the grid values to a file
void saveGridMarkFile(grid G, char *filename);  // saves the grid markings to a file
// void addMarkToGrid(grid G, char *filename); // adds markings to a grid from a file

grid initGridLaby(int, int, int w);             // labyrinth x,y, w = corridor width
grid initGridPoints(int, int, int t, double p); // pts of texture t with proba p
grid initGridFile(char *);                      // builds a grid from a file
position randomPosition(grid, int t);           // random position on texture type t
void freeGrid(grid);                            // frees the memory allocated by the initGridXXX() functions
void debug(int rank, char *format, ...);        // debug function

#endif
