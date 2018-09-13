#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include "cgraph.h"
#include "memory.h"
#include "vectormgm.h"
#include "clique.h"
#include "macros.h"
#include "containers.h"

#define INI_CLQ_CAP 512

#define MAX_NAME_SIZE 64

#define MIN_CLIQUE_ROW 256 /* mininum size for a row to be considered a clique row */

#define HASH_SIZE 16

struct _CGraph
{
    /* per node stored conflicts
       at position i there is a set with
       all node conflicts */
    ISet **nodeConflicts;
    ISet **nodeCliques;  /* all cliques in which a node appears */
    int nodeSize;         /* number of nodes considered */

    CliqueSet *clqSet;

    /* degree of each node */
    int *degree;

    int minDegree;
    int maxDegree;
    int lowDegree;

    int *origIdx; /* if it is a preprocessed graph,
                    indicates for each node i its original node */

    /* node names and weights are optional information */
    char **nodeNames;

    int *w;
};

typedef struct
{
    int node;
    int cost;
} NodeCost;

struct _NeighIterator
{
    /* binary heap vector */
    NodeCost *bhv;
    int bhvSize;
    int bhvCap;
    int node;

    /* previously returned node */
    int previousNode;
};

typedef struct CGArc
{
    int tail;
    int head;
} CGArc;

void cgraph_node_resize( CGraph *cgraph, int requiredNodeSize );
void cgraph_check_clique_size( CGraph *cgraph, int requiredCliqueSize );
void cgraph_add_clique_as_normal_conflicts( CGraph *cgraph, const int nEl, const int elClique[] );

/** for every tail add all heads as neighbors. vector arcs should have capacit of nArcs+1 **/
void cgraph_add_vector_arcs( CGraph *cgraph, const CGArc *arcs, const int nArcs, int neighs[] );

CGraph* cgraph_clone(const CGraph *cg)
{
    CGraph *clone = xmalloc(sizeof(CGraph));

    clone->nodeSize = cg->nodeSize;
    clone->nodeConflicts = xmalloc(sizeof(ISet*) * clone->nodeSize);
    clone->nodeCliques = xmalloc(sizeof(ISet*) * clone->nodeSize);
    for(int i = 0; i < clone->nodeSize; i++) {
        if(cg->nodeConflicts[i]) {
            clone->nodeConflicts[i] = iset_create(HASH_SIZE);
            iset_cpy(clone->nodeConflicts[i], cg->nodeConflicts[i]);
        } else {
            clone->nodeConflicts[i] = NULL;
        }

        if(cg->nodeCliques[i]) {
            clone->nodeCliques[i] = iset_create(HASH_SIZE);
            iset_cpy(clone->nodeCliques[i], cg->nodeCliques[i]);
        } else {
            clone->nodeCliques[i] = NULL;
        }
    }

    clone->clqSet = clq_set_clone(cg->clqSet);

    clone->degree = xmalloc(sizeof(int)* clone->nodeSize);
    memcpy(clone->degree, cg->degree, sizeof(int)* clone->nodeSize);
    clone->minDegree = cg->minDegree;
    clone->maxDegree = cg->maxDegree;
    clone->lowDegree = cg->lowDegree;

    if(cg->origIdx)
    {
        clone->origIdx = xmalloc(sizeof(int) * clone->nodeSize);
        memcpy(clone->origIdx, cg->origIdx, sizeof(int) * clone->nodeSize);
    }
    else
        clone->origIdx = NULL;

    if(cg->nodeNames)
    {
        clone->nodeNames = xmalloc(sizeof(char*) * clone->nodeSize);
        clone->nodeNames[0] = xmalloc(sizeof(char) * MAX_NAME_SIZE * clone->nodeSize);
        for(int i = 1; i < clone->nodeSize; i++)
            clone->nodeNames[i] = clone->nodeNames[i-1] + MAX_NAME_SIZE;
        for(int i = 0; i < clone->nodeSize; i++)
            strncpy(clone->nodeNames[i], cg->nodeNames[i], MAX_NAME_SIZE);
    }
    else
        clone->nodeNames = NULL;

    if(cg->w)
    {
        clone->w = xmalloc(sizeof(int) * clone->nodeSize);
        memcpy(clone->w, cg->w, sizeof(int) * clone->nodeSize);
    }
    else
        clone->w = NULL;

    return clone;
}

CGraph *cgraph_preprocess( const CGraph *cgraph, int nindexes[] )
{
    const int gSize = cgraph_size( cgraph );
    /* pre processed graph size */
    int ppgSize = 0;
    char recomputeDegree = 0;

    int nConflicts;
    const int confSpace = gSize*100;
    int *conflicts = malloc( sizeof(int)*confSpace );

    /* checking how many columns will remain */
    {
        int j;
        for ( j=0 ; ( j<gSize ) ; ++j )
        {
            nindexes[j] = -1;

            if ( cgraph_degree(cgraph, j)<3 )
                continue;

            nindexes[ j ] = ppgSize++;
        }
    }

    CGraph *result = NULL;
    result = cgraph_create( ppgSize );

    /* adding new conflicts */
    {
        int j;
        for ( j=0 ; (j<gSize) ; ++j )
        {
            const int newNode = nindexes[j];
            if ( newNode == -1 )
                continue;

            if(!cgraph->nodeConflicts[j])
                continue;

            const int gConfSize = iset_n_elements(cgraph->nodeConflicts[j]);
            {
                int k;
                for ( k=0 ; (k<gConfSize) ; ++k )
                {
                    const int gConf = iset_element(cgraph->nodeConflicts[j], k);
                    const int newNodeConf = nindexes[ gConf ];
                    if ( newNodeConf == -1 )
                        continue;

                    iset_add(result->nodeConflicts[newNode], newNodeConf);
                }
            }
        }
    }

    /* adding new cliques */
    {
        int i;

        for ( i=0 ; (i<clq_set_number_of_cliques(cgraph->clqSet)) ; ++i )
        {
            nConflicts = 0;

            const int nEl = clq_set_clique_size( cgraph->clqSet,i );
            const int *elClique = clq_set_clique_elements( cgraph->clqSet,i );
            {
                int j;
                for ( j=0 ; (j<nEl) ; ++j )
                {
                    const int newNode = nindexes[ elClique[j] ];
                    if ( newNode == -1 )
                        continue;

                    conflicts[nConflicts++] = newNode;
                }
            }

            if ( nConflicts < 2 )
                continue;

            if ( nConflicts < MIN_CLIQUE_ROW )
                /* not adding as clique anymore */
                cgraph_add_clique_as_normal_conflicts( result, nEl, elClique );
            else
                /* large clique, adding as it is */
                cgraph_add_clique( result, conflicts, nConflicts );
            recomputeDegree = 1;
        }
    }

#ifdef DEBUG
    assert( result != NULL );
    assert( cgraph_size(result) == ppgSize );
    {
        int i;
        for ( i=0 ; (i<cgraph_size(result)) ; ++i )
        {
            nConflicts = cgraph_get_all_conflicting( result, i, conflicts, confSpace );
            assert( nConflicts >= 3 );
            int j;
            for ( j=0 ; (j<nConflicts) ; ++j )
            {
                assert( conflicts[j] < ppgSize );
                assert( conflicts[j] >= 0 );
            }
            assert( nConflicts >= 3);
        }
    }
#endif

    free(conflicts);

    if(recomputeDegree)
        cgraph_recompute_degree(result);
    else
        cgraph_update_min_max_degree( result );

    return result;
}

int cgraph_cmp_int( const void *e1, const void *e2 )
{
    int i1 = (*((const int*) e1));
    int i2 = (*((const int*) e2));

    return (i1-i2);
}

CGraph *cgraph_create( int initialColumns )
{
    CGraph *result = xmalloc( sizeof(CGraph) );

    result->nodeSize = initialColumns;
    result->nodeConflicts = xmalloc( sizeof(ISet*)*result->nodeSize );
    result->nodeCliques   = xmalloc( sizeof(ISet*)*result->nodeSize );
    {
        int i;
        for ( i=0 ; (i<result->nodeSize) ; ++i ) {
            result->nodeCliques[i] = NULL;
            result->nodeConflicts[i] = NULL;
        }
    }

    result->degree = xmalloc( sizeof(int)*initialColumns );
    memset( result->degree, 0, sizeof(int)*initialColumns );

    result->clqSet = clq_set_create();

    result->nodeNames = NULL;

    result->w = NULL;

    result->origIdx = NULL;

    result->lowDegree = 15;

    return result;
}

void cgraph_add_node_conflict( CGraph *cgraph, const int node1, const int node2 )
{
    int lastNode = MAX(node1, node2);
    ++lastNode;

    if(cgraph->nodeSize < lastNode)
        cgraph_node_resize( cgraph, lastNode );

    /*-----------node1 to node2-----------*/
    //check if iset was not created
    if(!cgraph->nodeConflicts[node1])
        cgraph->nodeConflicts[node1] = iset_create(HASH_SIZE);

    //computes the correct degree
    if (!cgraph_conflicting_nodes( cgraph, node1, node2 ))
        cgraph->degree[node1]++;

    //adding conflict
    iset_add(cgraph->nodeConflicts[node1], node2);
    /*------------------------------------*/

    /*-----------node2 to node1-----------*/
    //check if iset was not created
    if(!cgraph->nodeConflicts[node2])
        cgraph->nodeConflicts[node2] = iset_create(HASH_SIZE);

    //computes the correct degree
    if (!cgraph_conflicting_nodes( cgraph, node2, node1 ))
        cgraph->degree[node2]++;

    //adding conflict
    iset_add(cgraph->nodeConflicts[node2], node1);
    /*------------------------------------*/
}

void cgraph_add_node_conflicts( CGraph *cgraph, const int node, const int conflicts[], const int size )
{
	int i;

    /* checking node capacity */
    int lastNode = node;
    for ( i=0 ; (i<size) ; ++i )
        lastNode = MAX( lastNode, conflicts[i] );
    ++lastNode;

    if(cgraph->nodeSize < lastNode)
        cgraph_node_resize( cgraph, lastNode );

    //check if iset was not created
    if(!cgraph->nodeConflicts[node])
        cgraph->nodeConflicts[node] = iset_create(HASH_SIZE);

    /* computes the correct degree */
    for ( i=0 ; (i<size) ; ++i )
    {
        //check if iset was not created
        if(!cgraph->nodeConflicts[conflicts[i]])
            cgraph->nodeConflicts[conflicts[i]] = iset_create(HASH_SIZE);

        if (!cgraph_conflicting_nodes( cgraph, node, conflicts[i] ))
            cgraph->degree[node]++;
        if (!cgraph_conflicting_nodes( cgraph, conflicts[i], node ))
            cgraph->degree[conflicts[i]]++;

        //adding conflicts
        iset_add(cgraph->nodeConflicts[node], conflicts[i]);
        iset_add(cgraph->nodeConflicts[conflicts[i]], node);
    }
}

void cgraph_add_node_conflicts_no_sim( CGraph *cgraph, const int node, const int conflicts[], const int size )
{
    int lastNode = node;
    {
        int i;
        for ( i=0 ; (i<size) ; ++i )
            lastNode = MAX( lastNode, conflicts[i] );
    }
    ++lastNode;

    if(cgraph->nodeSize < lastNode)
        cgraph_node_resize( cgraph, lastNode );

    //check if iset was not created
    if(!cgraph->nodeConflicts[node])
        cgraph->nodeConflicts[node] = iset_create(HASH_SIZE);

    /* computes the correct degree */
    int i;
    for ( i=0 ; (i<size) ; ++i ) {
        //check if iset was not created
        if(!cgraph->nodeConflicts[conflicts[i]])
            cgraph->nodeConflicts[conflicts[i]] = iset_create(HASH_SIZE);

        if (!cgraph_conflicting_nodes(cgraph, node, conflicts[i]))
            cgraph->degree[node]++;

        iset_add(cgraph->nodeConflicts[node], conflicts[i]);
    }
}

void cgraph_add_clique( CGraph *cgraph, const int *clique, const int size )
{
    /* checking node capacity */
    int lastNode = -1;
    {
        int i;
        for ( i=0 ; (i<size) ; ++i )
            lastNode = MAX( lastNode, clique[i] );
    }
    ++lastNode;

    if(cgraph->nodeSize < lastNode)
        cgraph_node_resize( cgraph, lastNode );

    if (!(clq_set_add( cgraph->clqSet, size, clique, 0 )))
        return;

    /* computes the correct degree */
    /*const int sizeM1 = size - 1;
    int i, j;
    for ( i=0 ; (i<sizeM1) ; ++i )
       for ( j=i+1 ; (j<size) ; ++j )
       {
          if (!cgraph_conflicting_nodes( cgraph, clique[i], clique[j] ) )
             cgraph->degree[clique[i]]++;
          if (!cgraph_conflicting_nodes( cgraph, clique[j], clique[i] ) )
             cgraph->degree[clique[j]]++;
       }*/

    /* computes a estimated degree */
    {
        int i;
        for ( i=0 ; (i<size) ; ++i )
            cgraph->degree[ clique[i] ] += size-1;
    }

    /* this clique will be related to all nodes inside it */
    {
        int i;
        const int idxClique = clq_set_number_of_cliques(cgraph->clqSet) - 1;
        for ( i=0 ; (i<size) ; ++i )
        {
            if(cgraph->nodeSize < lastNode)
                cgraph_node_resize( cgraph, clique[i]+1 );

            //check if iset was not created
            if(!cgraph->nodeCliques[clique[i]])
                cgraph->nodeCliques[clique[i]] = iset_create(HASH_SIZE);

            iset_add(cgraph->nodeCliques[clique[i]], idxClique);
        }
    }
}

int cgraph_size( const CGraph *cgraph )
{
    return cgraph->nodeSize;
}

void cgraph_node_resize( CGraph *cgraph, int requiredNodeSize )
{
    int newNodeSize = (int)(((double)cgraph->nodeSize) * 1.5);
    if ( newNodeSize < requiredNodeSize )
        newNodeSize = requiredNodeSize;

    /* printf(">>>> node size changed to %d\n", cgraph->nodeSize ); */

    cgraph->nodeConflicts = xrealloc( cgraph->nodeConflicts, sizeof(ISet*)*newNodeSize );
    cgraph->nodeCliques = xrealloc( cgraph->nodeCliques, sizeof(ISet*)*newNodeSize );
    cgraph->degree = xrealloc( cgraph->degree, sizeof(int)*newNodeSize );
    memset( cgraph->degree+cgraph->nodeSize, 0, newNodeSize-cgraph->nodeSize );
    if (cgraph->w)
        cgraph->w = xrealloc( cgraph->w, sizeof(int)*newNodeSize );
    if (cgraph->nodeNames)
    {
        cgraph->nodeNames = xrealloc( cgraph->nodeNames, sizeof(char*)*newNodeSize );
        cgraph->nodeNames[0] = xrealloc( cgraph->nodeNames[0], sizeof(char)*newNodeSize*MAX_NAME_SIZE );
        int i;
        for ( i=1 ; (i<newNodeSize); ++i )
            cgraph->nodeNames[i] = cgraph->nodeNames[i-1] + MAX_NAME_SIZE;
    }

    cgraph->nodeSize = newNodeSize;
}

int cgraph_conflicting_nodes( const CGraph *cgraph, const int i, const int j )
{
    if(i == j)
        return 0;

    const ISet *cliques;
#ifndef NDEBUG
    assert( i>=0 );
    assert( j>=0 );
    assert( i<cgraph_size(cgraph) );
#endif
    if(!cgraph->nodeConflicts[i] || iset_n_elements(cgraph->nodeConflicts[i]) == 0)
        goto FIND_IN_CLIQUES;

    if(iset_has(cgraph->nodeConflicts[i], j ))
        return 1;

FIND_IN_CLIQUES:
    cliques = cgraph->nodeCliques[i];

    if(!cliques)
        return 0;

    const int nCliques = iset_n_elements(cliques);

    {
        int idx;
        for ( idx=0 ; (idx<nCliques ) ; ++idx ) {
            const int idxClique = iset_element(cliques, idx);
            if (clq_set_clique_has_element(cgraph->clqSet, idxClique, j)) {
#ifndef NDEBUG
                if (!clq_set_clique_has_element( cgraph->clqSet, idxClique, i ))
                {
                    fprintf( stderr, "Error:  node %d does not appears in clique %d\n", i, idxClique );
                }
                assert( clq_set_clique_has_element( cgraph->clqSet, idxClique, i ) );
                assert( clq_set_clique_has_element( cgraph->clqSet, idxClique, j ) );
#endif
                return 1;
            }
        }
    }

    return 0;
}

int cgraph_get_all_conflicting( const CGraph *cgraph, int node, int neighs[], int maxSize )
{
    int i, j;
    const int size = cgraph->nodeSize;
    char *iv = (char*) calloc(size, sizeof(char));
    const ISet *isnc = cgraph->nodeConflicts[node];
    int nConfs = (isnc) ? iset_n_elements(isnc) : 0;
	int currClique = 0, clqSize = 0;

    if(nConfs > maxSize)
    	goto NO_SPACE;

    for(i = 0; i < nConfs; i++)
    {
        const int conf = iset_element(isnc, i);
    	iv[conf] = 1;
    	neighs[i] = conf;
    }

    /* now filling information from cliques, i.e., implicitly stored conflicts */
    const int nCliques = (cgraph->nodeCliques[node]) ? iset_n_elements(cgraph->nodeCliques[node]) : 0;
    for(currClique = 0; currClique < nCliques; currClique++)
    {
        const int idxClique = iset_element(cgraph->nodeCliques[node], currClique);
        const int *clqEl = clq_set_clique_elements(cgraph->clqSet, idxClique);
        clqSize = clq_set_clique_size(cgraph->clqSet, idxClique);

        for(j = 0; j < clqSize; j++)
            if(!iv[clqEl[j]] && clqEl[j] != node)
            {
                iv[clqEl[j]] = 1;
                neighs[nConfs++] = clqEl[j];
            }

        if(nConfs > maxSize)
            goto NO_SPACE;
    }

    free(iv);

    return nConfs;

NO_SPACE:
    fprintf( stderr, "ERROR: cgraph_get_all_conflicting:: Not enough space specified in maxSize.\n" );
    fprintf( stderr, "Working with node %d, which appears in %d cliques. When adding clique %d size %d. Result %d. MaxSize %d.\n", node, iset_n_elements(cgraph->nodeCliques[node]), currClique, clqSize, nConfs, maxSize );
    fprintf( stderr, "at: %s:%d\n", __FILE__, __LINE__ );
    abort();
    exit( EXIT_FAILURE );
}

int compare_int (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

void cgraph_save( CGraph *cgraph, const char *fileName )
{
    int *w = cgraph->w;

    FILE *f = fopen( fileName, "w" );
    if (!f)
    {
        fprintf( stderr, "Could not open file %s.", &(fileName[0]) );
        exit( EXIT_FAILURE );
    }

    /* counting edges */
    int nEdges = 0;
    int nodes = cgraph_size( cgraph );

    int *degree = xmalloc( sizeof(int)*nodes );
    memset( degree, 0, sizeof(int)*nodes );

    {
        int i,j;
        for ( i=0 ; (i<nodes-1) ; ++i )
            for ( j=i+1 ; (j<nodes) ; ++j )
                if ( cgraph_conflicting_nodes( cgraph, i, j ) )
                {
                    degree[i]++;
                    degree[j]++;
                }
        for ( i=0 ; (i<nodes-1) ; ++i )
            for ( j=i+1 ; (j<nodes) ; ++j )
                if ( ( cgraph_conflicting_nodes( cgraph, i, j ) )  )
                    nEdges++;

    }

    fprintf( f, "p edges %d %d\n", cgraph_size( cgraph ), nEdges );

    {
        int i, j;
        for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
        {
            const ISet *is = cgraph->nodeConflicts[i];
            const int nsize = (is) ? iset_n_elements(is) : 0;
            
            if(nsize == 0)
            	continue;

            int *confs = xmalloc(sizeof(int) * (nsize));
            for ( j=0 ; (j<nsize) ; ++j )
                confs[j] = iset_element(is, j);

            qsort( confs, nsize, sizeof(int), compare_int );
            
            for ( j=0 ; (j<nsize) ; ++j )
                if (confs[j] > i)
                    fprintf(f, "e %d %d\n", i + 1, confs[j] + 1);
            free(confs);
        }
    }

    if (w)
    {
        {
            int i;
            for ( i=0 ; (i<nodes) ; ++i )
                fprintf( f, "w %d %d\n", i+1, w[i] );
        }
    }

    {
        int i,j;
        for ( i=0 ; (i<clq_set_number_of_cliques( cgraph->clqSet)) ; ++i )
        {
            const int size = clq_set_clique_size( cgraph->clqSet, i );
            const int *elements = clq_set_clique_elements( cgraph->clqSet, i );
            fprintf( f, "c %d\n", size );

            int *confs = xmalloc(sizeof(int) * (size));
            for ( j=0 ; (j<size ) ; ++j )
            	confs[j] = elements[j];

            qsort( confs, size, sizeof(int), compare_int );

            for ( j=0 ; (j<size ) ; ++j )
                fprintf( f, "%d\n", elements[j]+1 );
            free(confs);
        }
    }

    if ( cgraph->nodeNames )
    {
        int i;
        for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
        {
            if ( cgraph_get_node_name( cgraph, i ) )
            {
                char name[256];
                strcpy( name, cgraph_get_node_name( cgraph, i ) );
                if ( name[ strlen(name)-1 ] == '\n' )
                    name[ strlen(name)-1 ] = '\0';
                fprintf( f, "n %d %s\n", i+1, name );
            }
        }
    }

    if ( cgraph->origIdx )
    {
        int i;
        for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
            fprintf( f, "o %d %d\n", i+1, cgraph->origIdx[i]+1 );
    }

    free( degree );

    fclose( f );
}

int compArc( const void *p1, const void *p2 )
{
    if ( ((const CGArc*)p1)->tail != ((const CGArc*)p2)->tail )
        return ((const CGArc*)p1)->tail - ((const CGArc*)p2)->tail;

    return ((const CGArc*)p1)->head - ((const CGArc*)p2)->head;
}

void cgraph_add_vector_arcs( CGraph *cgraph, const CGArc *arcs, const int nArcs,int neighs[] )
{
    if ( nArcs <= 0 )
        return;

    neighs[ 0 ] = arcs[0].head;
    int nNeighs = 1;
    int lastNode = arcs[0].tail;
    const int nArcsP1 = nArcs+1;

    {
        int i;
        for ( i=1 ; (i<(nArcsP1)) ; ++i )
        {
            if ( lastNode != arcs[i].tail )
            {
                assert( lastNode >= 0 );
                assert( lastNode < cgraph_size( cgraph ) );
                cgraph_add_node_conflicts_no_sim( cgraph, lastNode, neighs, nNeighs );
                lastNode = arcs[i].tail;
                nNeighs = 0;
                if ( i == nArcs )
                    break;
            }

            assert( arcs[ i ].head >= 0 );
            assert( arcs[ i ].head < cgraph_size( cgraph ) );
            neighs[nNeighs++] = arcs[ i ].head;
        }
    }
}

CGraph *cgraph_load( const char *fileName )
{
    int nodes = -1;
    int edges = -1;
    char recomputeDegree = 0;

    FILE *f = fopen( fileName, "r" );
    if (!f)
    {
        fprintf( stderr, "Could not open file %s.\n", fileName );
        exit( EXIT_FAILURE );
    }

    char line[LINE_SIZE];
    while ( fgets( line, LINE_SIZE, f ) )
    {
        if ( (line[0] == 'p') ||  (line[0] == 'P') )
        {
            char *ptr = line;
            char *end = ptr + LINE_SIZE;
            while ( (!isdigit(*ptr)) && (ptr<end) )
                ++ptr;
            sscanf( ptr, "%d %d", &nodes, &edges );
            break;
        }
    }

    if (nodes==-1)
    {
        fprintf( stderr, "Invalid file format.\n" );
        exit( EXIT_FAILURE );
    }

    CGArc *arcs = NULL;
    int arcCap = 8192;
    int nArcs = 0;
    arcs = xmalloc( sizeof(CGArc)*(arcCap+1) );
    int *clqElement = xmalloc( sizeof(int)*nodes );
    CliqueSet *clqSet = clq_set_create();
    int weightsLoaded = 0; /* lineNumber = 1; */

    char **names = NULL;
    int *tw = NULL;

    while ( fgets( line, LINE_SIZE, f ) )
    {
        if ( (line[0] == 'e') || (line[0] == 'E') )
        {
            int tail, head;
            sscanf( line+2, "%d %d", &tail, &head );
            --tail;
            --head;
            assert( tail >= 0 );
            assert( tail < nodes );
            assert( head >= 0 );
            assert( head < nodes );
            if ( tail > head )
            {
                int t = tail;
                tail = head;
                head = t;
            }

            /*         ADJUST_VECTOR_CAPACITY( arcs, arcCap, nArcs+2 );   */
            vmg_adjust_vector_capacity( (void**) (&arcs), &arcCap, nArcs+2, sizeof(CGArc) );

            arcs[nArcs].tail = tail;
            arcs[nArcs].head = head;
            nArcs++;
        }
        else
        {
            if ( ( (line[0] == 'w') || (line[0] == 'W') ) )
            {
                if (!tw)
                {
                    tw = xmalloc( sizeof(int)*nodes );
                }

                int nodeIdx, nodeWeight;
                sscanf( line+2, "%d %d", &nodeIdx, &nodeWeight );
                nodeIdx--;

                assert( ( nodeIdx >=0 ) && ( nodeIdx < nodes ) );
                tw[nodeIdx] = nodeWeight;

                weightsLoaded = 1;
            }
            else
            {
                if ( ( (line[0] == 'c') || (line[0] == 'C') ) )
                {
                    assert( weightsLoaded );
                    int clqSize, totalWeight = 0;
                    sscanf( line+2, "%d", &clqSize );
                    {
                        int i;
                        for ( i=0 ; (i<clqSize) ; ++i )
                        {
                            if ( fscanf( f, "%d", clqElement+i ) != 1 )
                            {
                                fprintf( stderr, "\nERROR: Expected clique element not found.\n" );
                                exit( EXIT_FAILURE );
                            }
                            clqElement[i]--;
                            assert( clqElement[i]>=0 );
                            assert( clqElement[i]<nodes );
                            totalWeight += tw[ clqElement[i] ];
                        }
                    }
                    clq_set_add( clqSet, clqSize, clqElement, totalWeight );
                    recomputeDegree = 1;
                } /* reading clique */
                else
                {
                    /* reading column names */
                    if ( ( (line[0] == 'n') || (line[0] == 'N') ) )
                    {
                        if ( names == NULL )
                        {
                            names = xmalloc( sizeof(char*)*nodes );
                            names[0] = xmalloc( sizeof(char)*MAX_NAME_SIZE*nodes );
                            int i;
                            for ( i=1 ; (i<nodes) ; ++i )
                                names[i] = names[i-1] + MAX_NAME_SIZE;
                        }

                        int idxNode;
                        sscanf( line+2, "%d", &idxNode );
                        int digits = 1;
                        idxNode--;
                        if ( idxNode != 0 )
                            digits = (int) ( log10( ((float)idxNode) ) + 1.0 );

                        strncpy( names[idxNode], line + 2 + digits + 1, MAX_NAME_SIZE );
                        int len = strlen( names[idxNode] );
                        if ( len > 0 )
                            if ( names[idxNode][len-1]=='\n' )
                                names[idxNode][len-1] = '\0';
                    }
                } /* not clique */
            } /* not weight */
        } /* not edge */
    }

    qsort( arcs, nArcs, sizeof(CGArc), compArc );

    /* removing repetitions */
    int nRep = 0;
    CGArc lastArc = { -1, -1 };
    {
        int i;
        for ( i=0 ; (i<nArcs) ; ++i )
            if ( compArc( arcs+i, &lastArc ) == 0 )
            {
                arcs[i].tail = INT_MAX;
                arcs[i].head = INT_MAX;
                nRep++;
            }
            else
                lastArc = arcs[i];
    }

    if (nRep)
        qsort( arcs, nArcs, sizeof(CGArc), compArc );

    nArcs -= nRep;

    CGraph *cgraph = cgraph_create( nodes );

    if (names)
        cgraph->nodeNames = names;
    if (tw)
    {
        cgraph->w = xmalloc( sizeof(int)*nodes );
        memcpy( cgraph->w, tw, sizeof(int)*nodes );
        free( tw );
    }

    int *neighs = malloc( sizeof(int)*nArcs );

    arcs[nArcs].tail = INT_MAX;
    arcs[nArcs].head = INT_MAX;

    cgraph_add_vector_arcs( cgraph, arcs, nArcs, neighs );

    /* inverting head and tail */
    {
        int i;
        for ( i=0 ; (i<nArcs) ; ++i )
        {
            const int t = arcs[i].tail;
            arcs[i].tail = arcs[i].head;
            arcs[i].head = t;
        }
    }

    qsort( arcs, nArcs, sizeof(CGArc), compArc );

    cgraph_add_vector_arcs( cgraph, arcs, nArcs, neighs );

    {
        int *t = xmalloc( sizeof(int)*nodes );
        int i;
        for ( i=0 ; (i<clq_set_number_of_cliques(clqSet)) ; ++i )
        {
            const int size = clq_set_clique_size(clqSet, i);
            const int *el = clq_set_clique_elements(clqSet, i);
            memcpy( t, el, sizeof(int)*size );
            cgraph_add_clique( cgraph, t, size );
        }
        free(t);
    }

    fclose( f );
    free( arcs );
    free( neighs );
    free(clqElement);
    clq_set_free( &clqSet );

    if(recomputeDegree)
        cgraph_recompute_degree(cgraph);
    else
        cgraph_update_min_max_degree( cgraph );

    return cgraph;
}

void cgraph_print( CGraph *cgraph, const int w[] )
{
    int *neighs = xmalloc( sizeof(int)*cgraph_size(cgraph)*100  );
    {
        int i;
        for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
        {
            printf("[%d] ", i+1 );
            int n = cgraph_get_all_conflicting( cgraph, i, neighs, cgraph_size(cgraph)*100 );
            {
                int j;
                for ( j=0 ; (j<n) ; ++j )
                    printf("%d ", neighs[j]+1 );
                printf("\n");
            }
        }
    }
    if (w)
    {
        int i;
        for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
            printf( "w[%d] %d\n", i+1, w[i] );
    }

    free( neighs );
}

void cgraph_update_min_max_degree( CGraph *cgraph )
{
    int i;

    cgraph->minDegree = INT_MAX;
    cgraph->maxDegree = 0;
    for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
    {
        cgraph->minDegree = MIN( cgraph->minDegree, cgraph->degree[i] );
    	cgraph->maxDegree = MAX( cgraph->maxDegree, cgraph->degree[i] );
    }
}

int cgraph_degree( const CGraph *cgraph, const int node )
{
    return cgraph->degree[ node ];
}

int cgraph_min_degree( const CGraph *cgraph )
{
    return cgraph->minDegree;
}

int cgraph_max_degree( const CGraph *cgraph )
{
    return cgraph->maxDegree;
}

void cgraph_free( CGraph **cgraph )
{
    int i;
    for ( i=0 ; (i<(*cgraph)->nodeSize) ; ++i ) {
        if((*cgraph)->nodeConflicts[i])
            iset_free(&(*cgraph)->nodeConflicts[i]);
        if((*cgraph)->nodeCliques[i])
            iset_free(&(*cgraph)->nodeCliques[i]);
    }

    free( (*cgraph)->nodeConflicts );
    free( (*cgraph)->nodeCliques );

    clq_set_free( &((*cgraph)->clqSet) );

    free( (*cgraph)->degree );

    if ((*cgraph)->origIdx)
        free( (*cgraph)->origIdx );

    if ((*cgraph)->nodeNames)
    {
        free((*cgraph)->nodeNames[0]);
        free((*cgraph)->nodeNames);
        (*cgraph)->nodeNames = NULL;
    }

    if ((*cgraph)->w)
        free((*cgraph)->w);

    free( (*cgraph) );

    (*cgraph) = NULL;
}

void cgraph_load_dimensions( const char *fileName, int *nodes, int *edges )
{
    FILE *f = fopen( fileName, "r" );
    if (!f)
    {
        fprintf( stderr, "Could not open file %s\n", &(fileName[0]) );
        exit( EXIT_FAILURE );
    }

    *nodes = -1;
    *edges = -1;

    char line[LINE_SIZE];
    while ( fgets( line, LINE_SIZE, f ) )
    {
        if ( (line[0] == 'p') || (line[0] == 'P') )
        {
            const char *p;
            if ( ((p = strstr( line, "edge" ))!=NULL) )
            {
                int shift = 7;
                if ( p[5] == 's' )
                    ++shift;
                if ( sscanf( line+shift, "%d %d", nodes, edges ) != 2 )
                {
                    fprintf( stderr, "Error reading dimensions.\n" );
                    exit( EXIT_FAILURE );
                }
            }

        }
    }

    fclose( f );
}

CGraph *cgraph_create_induced_subgraph( const CGraph *cgraph, const int nindexes[] )
{
    const int n = cgraph_size( cgraph );
    char recomputeDegree = 0;

    int newN = 0;

    int i;
    for ( i=0 ; (i<n) ; ++i )
        if ( nindexes[i] != -1 )
            ++newN;

    CGraph *result = cgraph_create( newN );
    if ( result->origIdx )
        free( result->origIdx );
    result->origIdx = xmalloc( sizeof(int)*newN );

    {
        int last = 0;
        for ( i=0 ; (i<n) ; ++i )
            if ( nindexes[i] != -1 )
                result->origIdx[ last++ ] = i;
    }

    int *neighs = xmalloc( sizeof(int)*newN );


    if (cgraph->w)
        result->w = xmalloc( sizeof(int)*newN );
    if (cgraph->nodeNames)
    {
        result->nodeNames = xmalloc( sizeof(char*)*cgraph_size(result) );
        result->nodeNames[0] = xmalloc( sizeof(char)*MAX_NAME_SIZE*cgraph_size(result) );
        int i;
        for ( i=1 ; (i<cgraph_size(result)) ; ++i )
            result->nodeNames[i] = result->nodeNames[i-1] + MAX_NAME_SIZE;
    }

    /* filling new conflicts */
    for ( i=0 ; (i<n) ; ++i )
    {
        if ( nindexes[i] == -1 )
            continue;

        int nNeighs = 0;

        const ISet *nodeConflicts = cgraph->nodeConflicts[i];
        const int nc = (nodeConflicts) ? iset_n_elements( nodeConflicts ) : 0;

        int j;
        for ( j=0 ; (j<nc) ; ++j )
        {
            const int c = iset_element(nodeConflicts, j);
            if(nindexes[c] != -1)
                neighs[nNeighs++] = nindexes[c];
        }

        if(nNeighs > 0)
            cgraph_add_node_conflicts_no_sim( result, nindexes[i], neighs, nNeighs );

        if ((cgraph->w)&&(result->w))
            result->w[nindexes[i]] = cgraph->w[i];
        if ((cgraph->nodeNames)&&(result->nodeNames))
            strncpy( result->nodeNames[nindexes[i]], cgraph->nodeNames[i], MAX_NAME_SIZE );
    }

    /* filling new cliques */
    const int nCliques = clq_set_number_of_cliques( cgraph->clqSet );
    for ( i=0 ; (i<nCliques) ; ++i )
    {
        const int nEl = clq_set_clique_size( cgraph->clqSet, i );
        const int *el = clq_set_clique_elements( cgraph->clqSet, i );

        int newClqSize = 0;
        int j;
        for ( j=0 ; (j<nEl) ; ++j )
        {
            if ( nindexes[ el[j] ] != -1 )
                neighs[ newClqSize++ ] = nindexes[ el[j] ];
        }

        if ( newClqSize > 1 )
        {
            if ( newClqSize < MIN_CLIQUE_ROW )
                /* not adding as clique anymore */
                cgraph_add_clique_as_normal_conflicts( result, newClqSize, neighs );
            else
                cgraph_add_clique( result, neighs, newClqSize );
            recomputeDegree = 1;
        }
    }

    free( neighs );

    if(recomputeDegree)
        cgraph_recompute_degree(result);
    else
        cgraph_update_min_max_degree( result );

    return result;
}

int cgraph_get_node_weight( const CGraph *cgraph, int node )
{
    if (cgraph->w)
        return cgraph->w[node];

    return 0;
}

void cgraph_set_node_weight( CGraph *cgraph, int node, int weight )
{
    if (!cgraph->w)
        cgraph->w = xmalloc( sizeof(int)*cgraph_size(cgraph));

    cgraph->w[node] = weight;
}

const char *cgraph_get_node_name( const CGraph *cgraph, int node )
{
    if (cgraph->nodeNames)
        return cgraph->nodeNames[node];

    return NULL;
}

void cgraph_set_node_name( CGraph *cgraph, int node, const char *name )
{
    if (!cgraph->nodeNames)
    {
        cgraph->nodeNames = xmalloc( sizeof(char*)*cgraph_size(cgraph) );
        cgraph->nodeNames[0] = xmalloc( sizeof(char)*MAX_NAME_SIZE*cgraph_size(cgraph) );
        int i;
        for ( i=1 ; (i<cgraph_size(cgraph)) ; ++i )
            cgraph->nodeNames[i] = cgraph->nodeNames[i-1] + MAX_NAME_SIZE;
    }

    strncpy( cgraph->nodeNames[node], name, MAX_NAME_SIZE );
}

const int *cgraph_get_node_weights( const CGraph *cgraph )
{
    return cgraph->w;
}

int cgraph_weight( const double w )
{
    return (int)(w*1000.0);
}

int cgraph_get_original_node_index( const CGraph *cgraph, const int node )
{
    if (cgraph->origIdx)
        return cgraph->origIdx[node];

    return -1;
}


#ifndef NDEBUG
void cgraph_check_node_cliques( const CGraph *cgraph )
{
    /* all nodes */
    int i,j;
    for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
    {
        const int nCliques = iset_n_elements(cgraph->nodeCliques[i]);
        const int *el = iset_elements( cgraph->nodeCliques[i] );
        for ( j=0 ; ( j<nCliques ) ; ++j )
        {
            const int idx = el[j];
            if (!(clq_set_clique_has_element( cgraph->clqSet, idx, i )))
            {
                printf("\nnode %d should appear on clique %d but does not\n", i, idx );
                fflush(stdout);
            }
            assert( clq_set_clique_has_element( cgraph->clqSet, idx, i ) );
        }
    }

    /* printf("Information in node cliques indicate cliques which really have the node.\n"); */

    const int nc = clq_set_number_of_cliques( cgraph->clqSet );
    for ( i=0 ; (i<nc) ; ++i )
    {
        const int n   = clq_set_clique_size( cgraph->clqSet, i );
        const int *el = clq_set_clique_elements( cgraph->clqSet, i );

        for ( j=0 ; (j<n) ; ++j )
        {
            const int currNode = el[j];
            /* this must appear in node clique */
            int l;
            int nn = iset_n_elements( cgraph->nodeCliques[currNode] );
            const int *eel = iset_elements( cgraph->nodeCliques[currNode] );

            for ( l=0 ; (l<nn) ; ++l )
                if (eel[l] == i)
                    break;

            if ( l == nn )
            {
                fprintf( stderr, "ERROR: in clique %d node %d appears but this clique does not appears in its node list.\n", i, currNode );
                fflush( stdout );
                fflush( stderr );
            }

            assert( l<nn );
        }
    }
}

void cgraph_check_neighs( const CGraph *cgraph )
{
    const int n = cgraph_size( cgraph );
    int nNeighs;
    int *neighs = xmalloc( sizeof(int)*n*100 );

    int i, j;
    for ( i=0 ; (i<n) ; ++i )
    {
        nNeighs = cgraph_get_all_conflicting( cgraph, i, neighs, n*100 );

        /* computed number of neighs */
        int cn = 0;
        for ( j=0 ; ( j<n ) ; ++j )
            if (cgraph_conflicting_nodes( cgraph, i, j ))
                cn++;

        assert( cn == nNeighs );

        for ( j=0 ; ( j<n ) ; ++j )
            if (cgraph_conflicting_nodes( cgraph, i, j ))
                assert ( bsearch( &j, neighs, nNeighs, sizeof(int), iset_cmp_int ) );
            else
                assert ( !bsearch( &j, neighs, nNeighs, sizeof(int), iset_cmp_int ) );
    }

    free( neighs );
}

void cgraph_check_preproc( const CGraph *ppgraph, const CGraph *cgraph )
{
    assert( ppgraph->origIdx );
    assert( cgraph_size(cgraph) >= cgraph_size(ppgraph) );
    const int n = cgraph_size( ppgraph );
    int i,j;

    int *ppIdx = xmalloc( sizeof(int)*cgraph_size(cgraph) );
    for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
        ppIdx[i] = -1;

    for ( i=0 ; (i<cgraph_size(ppgraph)) ; ++i )
        ppIdx[ ppgraph->origIdx[i] ] = i;

    for ( i=0 ; (i<n) ; ++i )
    {
        for ( j=i+1 ; (i<n) ; ++i )
        {
            const int cij = cgraph_conflicting_nodes( ppgraph, i, j );
            const int cji = cgraph_conflicting_nodes( ppgraph, j, i );
            assert( cij == cji );

            const int oi = ppgraph->origIdx[i];
            assert( oi != -1 );
            const int oj = ppgraph->origIdx[j];
            assert( oj != -1 );

            const int cijo = cgraph_conflicting_nodes( cgraph, oi, oj );
            const int cjio = cgraph_conflicting_nodes( cgraph, oj, oi );

            assert( cijo == cjio );

            if ( cijo != cij )
            {
                fprintf( stderr, "Error: Conflicts different in original and ppgraph.\n" );
            }
            assert( cijo == cij );
            assert( cjio == cji );
        }
    }

    for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
    {
        if ( ppIdx[i] == -1 )
            continue;
        const int pIdxI = ppIdx[i];
        for ( j=i+1 ; (j<cgraph_size(cgraph)) ; ++j )
        {
            if ( ppIdx[j] == -1 )
                continue;
            const int pIdxJ = ppIdx[j];
            const int co = cgraph_conflicting_nodes( cgraph, i, j );
            const int ci = cgraph_conflicting_nodes( ppgraph, pIdxI, pIdxJ );
            if (co!=ci)
                printf("bra");
            assert( co==ci );
        }
    }


    free( ppIdx );
}

#endif

const int *cgraph_get_original_node_indexes( const CGraph *cgraph )
{
    return cgraph->origIdx;
}

void cgraph_print_summary( CGraph *cgraph, const char *name )
{
    cgraph_update_min_max_degree( cgraph );
    int nodes = cgraph_size( cgraph );
    int nSmallDegree = 0;
    int edges = 0;
    {
        int i;
        for ( i=0 ; (i<nodes) ; ++i )
        {
            edges += cgraph_degree( cgraph,  i );
            if ( cgraph_degree( cgraph, i ) <= cgraph->lowDegree )
                nSmallDegree++;
        }
    }
    printf("%s : n %d e %d - minD %d maxD %d smallD %d\n", name, cgraph_size(cgraph), edges, cgraph_min_degree(cgraph), cgraph_max_degree(cgraph), nSmallDegree );
}

void cgraph_set_low_degree( CGraph *cgraph, const int lowDegree )
{
    cgraph->lowDegree = lowDegree;
}

void cgraph_add_clique_as_normal_conflicts( CGraph *cgraph, const int nEl, const int elClique[] )
{
    int i,j,count;
    const int nconflicts=nEl-1;
    int *iv = xmalloc(sizeof(int) * cgraph->nodeSize);

    for ( i=0 ; (i<nEl) ; ++i )
    {
        const int node = elClique[i];
        count = 0;
        for ( j=0 ; (j<nEl) ; ++j )
        {
            if ( elClique[j] == node )
                continue;

            iv[count++] = elClique[j];
        }
        cgraph_add_node_conflicts_no_sim( cgraph, node, iv, nconflicts );
    }

    free(iv);
}

int cgraph_get_n_conflicting( const CGraph *cgraph, int node, int neighs[], int n,
                              int v[], const int vcap )
{
    return n;
}


/* NeighIterator functions */

NeighIterator *nit_create( )
{
    NeighIterator *result = xmalloc( sizeof(NeighIterator) );

    result->bhv    = NULL;
    result->bhvCap = 0;

    result->previousNode = -1;

    return result;
}

void nit_check_capacity( NeighIterator *nit, int required )
{
    required++;

    if ( nit->bhvCap < required )
    {
        if ( nit->bhvCap )
        {
            nit->bhvCap = MAX( required, nit->bhvCap*2 );
            nit->bhv = xrealloc( nit->bhv, sizeof(NodeCost)*nit->bhvCap );
        }
        else
        {
            nit->bhvCap = MAX( required, 1024 );
            nit->bhv = xmalloc( sizeof(NodeCost)*nit->bhvCap );
        }
    }
}

/* Heap macros */
/* position of root node in vector */
#define rootPosBHV( node )  ( (((node+1)/2)-1) )
/* position of the first child node in vector */
#define childPosBHV( node ) ( (node*2)+1 )

#define bhv_swap( bhv,  p1, p2 ) { NodeCost tempBHV=bhv[p1]; bhv[p1]=bhv[p2]; bhv[p2]=tempBHV; }

void nit_bhv_up( NodeCost *bhv, int index )
{
    int root, child = index;

    while ( (root=rootPosBHV(child)) >=0 )
    {
        if ( bhv[root].cost > bhv[child].cost )
        {
            bhv_swap( bhv, child, root );
            child = root;
        }
        else
            return;
    }
}

void nit_bhv_down( NodeCost *bhv, int size, int index )
{
    int root = index, child;

#ifndef NDEBUG
    assert( index >= 0 );
#endif

    while ( (child=childPosBHV(root)) < size )
    {
#ifdef DEBUG
        assert( child >= 0 );
        assert( child < size );
#endif
        /* child with the smallest cost */
        if ( (child+1<size) && (bhv[child].cost>bhv[child+1].cost) )
            child++;

#ifdef DEBUG
        assert( child >= 0 );
        assert( child < size );
#endif
        if ( bhv[root].cost > bhv[child].cost )
        {
            bhv_swap( bhv, root, child );
            root = child;
        }
        else
            break;
    }
}

void nit_bhv_add( NeighIterator *nit, int node, int cost )
{
#ifdef DEBUG
    assert( node >= 0 );
    assert( (nit->bhvSize>=0) && (nit->bhvSize<nit->bhvCap) );
#endif


    nit->bhv[nit->bhvSize].node = node;
    nit->bhv[nit->bhvSize].cost = cost;

    nit_bhv_up( nit->bhv, nit->bhvSize );

    nit->bhvSize++;
}

void nit_fill_bhv( NeighIterator *nit, const CGraph *cgraph, int node, const int costs[] )
{
#ifdef DEBUG
    assert( node >= 0 );
    assert( node < cgraph_size(cgraph) );
#endif
    {
        const ISet *nodeConflicts = cgraph->nodeConflicts[node];
        const int nConf = (nodeConflicts) ? iset_n_elements( nodeConflicts ) : 0;
        int i;
        for ( i=0 ; (i<nConf) ; ++i )
        {
            const int conf = iset_element(nodeConflicts, i);
#ifdef DEBUG
            assert( (conf>=0) );
            assert( (conf<cgraph_size(cgraph)) );
#endif
            nit_bhv_add( nit, conf, costs[conf] );
        }
    }

    {
        /* conflicts stored in cliques */
        const int nCliques = (cgraph->nodeCliques[node]) ? iset_n_elements(cgraph->nodeCliques[node]) : 0;
        const CliqueSet *clqSet = cgraph->clqSet;
        int i;
        for ( i=0 ; (i<nCliques) ; ++i )
        {
            const int cliqueIdx = iset_element(cgraph->nodeCliques[node], i);
            const ISet *clique = clq_set_get_clique( clqSet, cliqueIdx );
            const int size = iset_n_elements(clique);
            const int *el = iset_elements(clique);
            int j;
            for ( j=0 ; (j<size) ; ++j )
                nit_bhv_add( nit, el[j], costs[el[j]] );
        }
    }
}

void nit_start( NeighIterator *nit, const CGraph *cgraph, int node, const int costs[] )
{
    /* nodespace : needed for store neighs and all cliques */
    int nodeSpace      = (cgraph->nodeConflicts[node]) ? iset_n_elements(cgraph->nodeConflicts[node]) : 0;
    int nCliques       = (cgraph->nodeCliques[node]) ? iset_n_elements(cgraph->nodeCliques[node]) : 0;
    const CliqueSet *clqSet = cgraph->clqSet;
    {
        int i;
        for ( i=0 ; (i<nCliques) ; ++i ) {
            const int clique = iset_element(cgraph->nodeCliques[node], i);
            nodeSpace += iset_n_elements(clq_set_get_clique(clqSet, clique));
        }
    }

    nit_check_capacity( nit, nodeSpace );

    nit->bhvSize      =  0;
    nit->previousNode = -1;

    nit_fill_bhv( nit, cgraph, node, costs );

    nit->node = node;
}

int nit_next( NeighIterator *nit )
{
    int result = -1;

NEXT:

    if (nit->bhvSize==0)
        return INT_MAX;

    NodeCost *bhv = nit->bhv;

    result = bhv[0].node;

    --nit->bhvSize;
    bhv[0] = bhv[nit->bhvSize];

    nit_bhv_down( bhv, nit->bhvSize, 0 );

    if ((result==nit->previousNode)||(result==nit->node))
        goto NEXT;

    nit->previousNode = result;
    return result;
}

void cgraph_recompute_degree( CGraph *cgraph )
{
    const int size = cgraph->nodeSize;
    const ISet *cliques;
    char *iv = xmalloc(sizeof(char) * size);

    cgraph->minDegree = INT_MAX;
    cgraph->maxDegree = 0;
    memset( cgraph->degree, 0, sizeof(int)*size );

    for(int i = 0; i < size; i++)
    {
        memset(iv, 0, sizeof(char)*size);

        const ISet *isnc = cgraph->nodeConflicts[i];
        const int nConfs = (isnc) ? iset_n_elements(isnc) : 0;

        //individual conflicts
        cgraph->degree[i] += nConfs;
        for(int j = 0; j < nConfs; j++) {
            const int conf = iset_element(isnc, j);
            iv[conf] = 1;
        }

        //conflicts stored as cliques
        cliques = cgraph->nodeCliques[i];
        const int nCliques = (cliques) ? iset_n_elements(cliques) : 0;
        for(int j = 0; j < nCliques; j++)
        {
            const int idxClique = iset_element(cliques, j);
            const int clqSize = clq_set_clique_size(cgraph->clqSet, idxClique);
            const int *clqEl = clq_set_clique_elements(cgraph->clqSet, idxClique);

            for(int k = 0; k < clqSize; k++)
                if(!iv[clqEl[k]] && clqEl[k] != i)
                {
                    iv[clqEl[k]] = 1;
                    cgraph->degree[i]++;
                }
        }

        //updating min and max degree
        cgraph->minDegree = MIN( cgraph->minDegree, cgraph->degree[i] );
        cgraph->maxDegree = MAX( cgraph->maxDegree, cgraph->degree[i] );
    }

    free(iv);
}

void nit_free( NeighIterator **nit )
{
    free( (*nit)->bhv );
    free( (*nit) );
    (*nit) = NULL;
}
