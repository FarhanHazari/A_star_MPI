#ifndef HEAP_H
#define HEAP_H

#include <stdbool.h>
#include "tools.h"

// "Node" structure for the min Q heap.
typedef struct node
{
    position pos;        // position (.x,.y) of node u
    double cost;         // cost[u]
    double score;        // score[u] = cost[u] + h(u,end)
    struct node *parent; // parent[u] = pointer to father, NULL for start
} *node;

// "MPI_Node" structure to send accross processses
typedef struct
{
    position pos;
    double cost;
    double score;
    int parent_rank;
    int parent_win_i;
} mpi_node;

// Binary heap structure:
//
//  array = storage array for objects starting at index 1 (instead of 0)
//  n     = number of objects (which are void*) stored in the heap
//  nmax  = maximum number of objects that can be stored in the heap
//  f     = comparison function for two objects (min, max, ..., see man qsort)
//
// Warning! "heap" is defined as a pointer to optimize calls (pushing
// one word (= 1 pointer) instead of 4 otherwise).
typedef struct heap
{
    void **array;
    int n, nmax;

    int (*f)(const void *, const void *);
} *heap;

// Creates a heap that can hold at most k>0 objects with a predefined
// comparison function f(). NB: The size of an object pointed to by a
// pointer h is sizeof(*h).
heap heap_create(int k, int (*f)(const void *, const void *));

// Destroys heap h. We will assume h!=NULL. Warning! This is about
// freeing what has been allocated by heap_create(). NB: The objects
// stored in the heap do not have to be freed.
void heap_destroy(heap h);

// Returns true if heap h is empty, false otherwise. We will assume
// h!=NULL.
int heap_empty(heap h);

// Adds an object to heap h. We will assume h!=NULL. Returns true if
// there is not enough space, and false otherwise.
bool heap_add(heap h, void *object);

// Returns the object at the top of heap h, that is, the minimal
// element according to f(), without deleting it. We will assume
// h!=NULL. Returns NULL if the heap is empty.
void *heap_top(heap h);

// Like heap_top() except that the object is also deleted from the
// heap. Returns NULL if the heap is empty.
void *heap_pop(heap h);

#endif
