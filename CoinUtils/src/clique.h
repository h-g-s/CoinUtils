#ifndef CLIQUE_H_INCLUDED
#define CLIQUE_H_INCLUDED

#include <vector>
#include "cgraph.h"

typedef struct _CliqueSet CliqueSet;

CliqueSet *clq_set_clone( const CliqueSet *clqSet );

/***
 * creates a cliqueset
 **/
CliqueSet *clq_set_create( );

/***
 * creates a cliqueset loading its contents from a file
 **/
CliqueSet *clq_set_load( const char *fileName );

/**
 * saves a cliqueset
 **/
void clq_set_save( const CGraph *cgraph, const CliqueSet *clqSet, const char *fileName );

/***
 * weight sum of all cliques in the set
 **/
int clq_set_weight_sum( const CliqueSet *clqSet );


/***
 * adds a new clique - returns 1 if success, 0 if the clique was alread inserted
 * nodes are sorted before insertion
 **/
int clq_set_add( CliqueSet *clqSet, const int *idxs, const int size, const int w );
int clq_set_add( CliqueSet *clqSet, const std::vector<int> &idxs, const int w );

int clq_set_add_cliques( CliqueSet *clqs_target, const CliqueSet *clqs_source );

/**
 * clears all clique set
 * contents
 **/
void clq_set_clear( CliqueSet *clqSet );

/***
 * gets clique size
 **/
int clq_set_clique_size( const CliqueSet *clqSet, const int clique );

/***
 * gets clique size
 **/
const int* clq_set_clique_elements( const CliqueSet *clqSet, const int clique );

/**
 * finds an element in clique
 **/
int clq_set_clique_has_element( const CliqueSet *clqSet, const int clique, const int element );

/**
 * return weight of clique
 **/
int clq_set_weight( const CliqueSet *clqSet, const int clique );

/**
 * for debugging purposes
 **/
void clq_set_print( const CliqueSet *clqSet );

int clq_set_number_of_cliques( const CliqueSet *clqSet );

/***
 * frees clique set memory
 **/
void clq_set_free( CliqueSet **clqSet );

void clq_set_add_using_original_indexes( CliqueSet *target, const CliqueSet *source, const int *orig );


/***
 * checks if the clique really is a clique
 * returns  1 if everything is ok, 0 otherwize
 * if it is NOT a clique, it informs
 * which nodes are not neoighbors in n1 and n2
 ***/
int clq_validate( const CGraph *cgraph, const int *idxs, const int size, int *n1, int *n2 );


#endif

