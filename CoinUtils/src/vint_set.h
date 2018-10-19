/** set of integers implemented as a vector
    fast for queries
    slow for multiple modifications  **/

#ifndef VINT_SET_H_INCLUDED
#define VINT_SET_H_INCLUDED

#include <vector>

typedef struct _IntSet IntSet;

IntSet* vint_set_create();

void vint_set_free(IntSet **_iset);

IntSet* vint_set_clone(const IntSet *iset);

const std::vector<int>& vint_set_get_elements( IntSet *iset );

const int vint_set_find(IntSet *iset, int key);

const int vint_set_intersection( int intersection[], const int size, const int elements[], IntSet *is2 );

int vint_set_size(IntSet *iset);

void vint_set_add( IntSet *iset, const int elements[], int size );

void vint_set_add_using_original_indexes( IntSet *iset, const int elements[], int size, const int orig[] );

/**
 * returns 1 if these int_set are equal, 0 otherwise
 **/
int vint_set_equals( IntSet *is1, IntSet *is2 );

int bsearch_int(const std::vector<int> &v, const int key);

void vint_insert_sort( const int key, std::vector<int> &v );

#endif

