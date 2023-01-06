#include "heap.h"
#include <stdlib.h>

heap heap_create(int k, int (*f)(const void *, const void *))
{
  heap h = malloc(sizeof(struct heap));
  h->array = malloc((k + 1) * sizeof(void *)); // k+1 because we start at 1
  h->n = 0;
  h->nmax = k;
  h->f = f;
  return h;
}

void heap_destroy(heap h)
{
  free(h->array);
  free(h);
}

int heap_empty(heap h)
{
  if (h->n == 0 || !h->array)
    return 1;
  return 0;
}

bool heap_add(heap h, void *obj)
{
  if (h->n == h->nmax)
  {
    int k = h->nmax * 2;
    void **array = realloc(h->array, (k + 1) * sizeof(void *));
    if (array == NULL)
    {
      return true;
    }
    else
    {
      h->array = array;
      h->nmax = k;
    }
  }

  // store object at the end
  h->array[++h->n] = obj;

  int son = h->n;
  int father = son / 2;
  while (father > 0 && h->f(h->array[father], h->array[son]) > 0)
  {
    node swap = h->array[father];
    h->array[father] = h->array[son];
    h->array[son] = swap;

    son = father;
    father = son / 2;
  }
  return false;
}

void *heap_top(heap h)
{
  if (!h->array)
    return NULL;

  return h->array[1];
}

void *heap_pop(heap h)
{
  if (h->n == 0 || !h->array)
    return NULL;

  // swap min and last
  node swap = h->array[h->n];
  h->array[h->n] = h->array[1];
  h->array[1] = swap;

  // store deleted (min)
  node deleted = h->array[h->n--];

  int father = 1;
  int left_son = (father * 2 > h->n) ? father : father * 2;          // in case n < 2 (left_son is empty)
  int right_son = (father * 2 + 1 > h->n) ? father : father * 2 + 1; // in case n < 3 (right_son is empty)
  int son = (h->f(h->array[left_son], h->array[right_son]) > 0) ? right_son : left_son;
  while (h->f(h->array[father], h->array[son]) > 0)
  {
    swap = h->array[father];
    h->array[father] = h->array[son];
    h->array[son] = swap;

    father = son;
    left_son = (father * 2 > h->n) ? father : father * 2;
    right_son = (father * 2 + 1 > h->n) ? father : father * 2 + 1;

    son = (h->f(h->array[left_son], h->array[right_son]) > 0) ? right_son : left_son;
  }

  return deleted;
}