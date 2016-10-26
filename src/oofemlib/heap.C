#include "heap.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

namespace oofem {
  
#define MIN(a, b) ( ( a ) > ( b ) ? ( b ) : ( a ) )

Heap :: Heap(int N) {
    // Heuristic estimation of heap size
    Initial_Heap_Alloc_Size = N / 4 + 1;

    heapCount = 0;
    allocatedSize = Initial_Heap_Alloc_Size;

    // Initialize Heap and H2T which are the same size.
    Keys = ( double * )       calloc( Initial_Heap_Alloc_Size,
                                      sizeof( double ) );
    H2T      = ( int * )        calloc( Initial_Heap_Alloc_Size,
                                        sizeof( int ) );

    T2H      = ( int * )      calloc( N,       sizeof( int ) );
    // Initialize T2H by setting all values to -1, signifying that
    // there are no heap elements corresponding to the pixels.
    for ( int i = 0; i < N; i++ ) {
        T2H [ i ] = -1;
    }
}

Heap :: ~Heap() {
    // Destructor
    free(Keys);
    free(H2T);
    free(T2H);
}


int Heap :: parentInd(int inInd) {
    return ( inInd - 1 ) / 2;
}

int Heap :: leftChildInd(int inInd) {
    return inInd * 2 + 1;
}

int Heap :: rightChildInd(int inInd) {
    return inInd * 2 + 2;
}

int Heap :: lastParentInd() {
    return ( heapCount - 2 ) / 2;
}




void Heap :: swapElements(int Ind1, int Ind2) {
    double tmpKey;
    int tmpInd;

    if ( Ind1 == Ind2 ) {
        return;
    }

    // Swap keys
    tmpKey = Keys [ Ind1 ];
    Keys [ Ind1 ] = Keys [ Ind2 ];
    Keys [ Ind2 ] = tmpKey;

    // Swap T2H values
    // NB: Must be done before the H2T swaps
    T2H [ H2T [ Ind1 ] ] = Ind2;
    T2H [ H2T [ Ind2 ] ] = Ind1;

    // Swap H2T elems
    tmpInd = H2T [ Ind1 ];
    H2T [ Ind1 ] = H2T [ Ind2 ];
    H2T [ Ind2 ] = tmpInd;
}

void Heap :: upHeap(int Ind) {
    // Index of parent of Ind
    int parent;

    while ( Ind > 0 ) {
        parent = ( Ind - 1 ) / 2;
        if ( Keys [ Ind ] < Keys [ parent ] ) {
            swapElements(Ind, parent);
            Ind = parent;
        } else   {
            break;
        }
    }
}

void Heap :: downHeap(int Ind) {
    // Index of first child. Add one to get second child.
    int child1, child2, minChild;

    // Special case (not more than 1 element in the heap)
    if ( heapCount < 2 ) {
        return;
    }

    // Loop until the biggest parent node
    while ( Ind <= ( ( heapCount - 2 ) / 2 ) ) {
        child1 = 2 * Ind + 1;
        child2 = child1 + 1;
        minChild = Ind;

        // Determine the child with the minimal value
        if ( Keys [ child1 ] < Keys [ Ind ] ) {
            minChild = child1;
        }
        if ( ( child2 < heapCount ) && ( Keys [ child2 ] < Keys [ minChild ] ) ) {
            minChild = child2;
        }

        // If there was a smaller child, swap with Ind
        if ( minChild != Ind ) {
            swapElements(Ind, minChild);
            Ind = minChild;
        } else   {
            break;
        }
    }
}


bool Heap :: isInHeap(int Ind) {
    // Since T2H is initialized to -1 for all elements,
    // any element that is not -1 (i.e. >= 0) is already in the heap.
    return ( T2H [ Ind ] >= 0 );
}

int Heap :: nElems() {
    return heapCount;
}

void Heap :: insert(double time, int Ind) {
    // Reallocate more memory if necessary
    if ( heapCount == ( allocatedSize - 1 ) ) {
        Keys = ( double * ) realloc( Keys,
                                     ( allocatedSize + Initial_Heap_Alloc_Size ) * sizeof( double ) );
        H2T  = ( int * )    realloc( H2T,
                                     ( allocatedSize + Initial_Heap_Alloc_Size ) * sizeof( int ) );
    }

    // Insert element at the end of the heap
    Keys [ heapCount ] = time;
    H2T [ heapCount ] = Ind;
    T2H [ Ind ] = heapCount;

    heapCount++;

    // Ensure the heap is maintained by moving the
    // new element as necessary upwards in the heap.
    upHeap(heapCount - 1);
}

void Heap :: update(double time, int Ind) {
    Keys [ T2H [ Ind ] ] = time;

    // Must do one upHead run because maybe this new, lower, time
    // is lower than those of the parents. By design of Fast-Marching
    // it should never be larger though.
    upHeap(T2H [ Ind ]);
}

double Heap :: getSmallest(int *Ind) {
    // Smallest element is the first of the heap
    double heapRoot = Keys [ 0 ];

    // Set 'pointer' to the Ind in T of the min heap element
    * Ind = H2T [ 0 ];

    // Replace root by the last heap element
    swapElements(0, heapCount - 1);

    // Decrement heapCount.
    // NB: Important to do this before downHeap is run!
    heapCount--;

    // Run downHeap from root.
    downHeap(0);

    return heapRoot;
}

int Heap :: checkHeapProperty(int pInd) {
    int lChild = leftChildInd(pInd);
    int rChild = lChild + 1;

    if ( ( lChild < heapCount ) && ( Keys [ lChild ] < Keys [ pInd ] ) ) {
        return pInd;
    }
    if ( ( rChild < heapCount ) && ( Keys [ rChild ] < Keys [ pInd ] ) ) {
        return pInd;
    }

    if ( ( lChild <= lastParentInd() ) ) {
        return checkHeapProperty(lChild);
    }
    if ( ( rChild <= lastParentInd() ) ) {
        return checkHeapProperty(rChild);
    }

    return -1;
}

void Heap :: print() {
    printf("\nHeap:  ");
    for ( int iter = 0; iter < heapCount; iter++ ) {
        printf("%.3g ", Keys [ iter ]);
    }
}

void Heap :: recurse(int row, int pad, int spacing, int S) {
    // If not the first row, recursivly call itself
    if ( row > 1 ) {
        int newSpacing = ( ( int ) ceil( 2.0 * ( ( double ) spacing ) +
                                         ( ( double ) S ) ) );
        int newPad = ( ( int ) ceil( ( ( double ) pad ) + ( ( double ) spacing ) / 2.0
                                     + ( ( double ) S ) / 2.0 ) );
        recurse(row - 1, newPad, newSpacing, S);
    }

    // Calculate the first and last elements on the row
    int beg = ( pow( 2, ( row - 1 ) ) - 1 );
    int end = MIN( ( heapCount - 1 ), ( pow(2, row) - 2 ) );

    // Newline and padding
    printf("\n");
    for ( int i = 0; i < pad; i++ ) {
        printf(" ");
    }

    // Print elements
    for ( int elem = beg; elem <= end; elem++ ) {
        printf("%-*.*g", S + spacing, S, Keys [ elem ]);
    }
}

void Heap :: printTree() {
    if ( heapCount == 0 ) {
        return;
    }

    // Constants
    int S = 3;
    int B = 4;

    int nRows = 1 + floor( log2(heapCount) );

    // Call recurse from the last row
    printf("\n");
    recurse(nRows, 0, B, S);
    printf("\n\n");
}

// Insert narrowBand pixels into a m x n matrix.
// Unused?
double *Heap :: formMatrix(int m, int n) {
    double *mat = ( double * ) calloc( m * n, sizeof( double ) );
    for ( int iter = 0; iter < heapCount; iter++ ) {
        mat [ H2T [ iter ] ] = Keys [ iter ];
    }
    return mat;
}
} // end namespace oofem
