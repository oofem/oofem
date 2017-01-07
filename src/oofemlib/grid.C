#include "grid.h"
#include "error.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

namespace oofem {

#define MIN(a, b) ( ( a ) > ( b ) ? ( b ) : ( a ) )
#define MAX(a, b) ( ( a ) > ( b ) ? ( a ) : ( b ) )

// NB: When using matrix indices (i,j), they run from 1 to (m,n).
// When using 1-dim array index 'ind', it runs from 0 to (m*n-1).

    // Index offsets of Patch
    int iOffsets_full[] = {
        -1, -1, -1,  0,  0,  1,  1,  1
    };
    int jOffsets_full[] = {
        -1,  0,  1, -1,  1, -1,  0,  1
    };

    // Index offsets of cross-neighbours
    int iOffsets[] = {
        -1, 1, 0, 0
    };
    int jOffsets[] = {
        0, 0, -1, 1
    };

  // Offsets used by calcTime
    // NB: Don't change the ordering without also changing the if-else
    // statements below that set xmin, ymin.
    int icalcOffsets[] = {
        -2, -1, 1, 2,    0,  0, 0, 0
    };
    int jcalcOffsets[] = {
        0,  0, 0, 0,   -2, -1, 1, 2
    };

  
  Grid :: Grid(int N, int M)
  {
    n = N;
    m = M;
    solutionAvailable = false;
    T = NULL;
    F = NULL;
    // (calloc sets all bits of the allocated memory to zero, so the initial values are "false")
    Frozen = ( bool * ) calloc( m * n, sizeof( bool ) );
    narrowBand = new Heap(m * n);
  }
  
  Grid :: ~Grid()
  {
    if (T) free(T);
    if (F) free(F);
    if (Frozen) free(Frozen);
    delete narrowBand;
  }

int
Grid :: giveSize(int dir)
{
  switch (dir){
  case 1: return n;
  case 2: return m;
  default: return 0;
  }
}
  
void
Grid :: setPrescribedField(double* PF, int size)
{
  // For efficiency, F is passed by its pointer and from now on
  // will be taken care of (and deleted) by the Grid.
  if (F) delete F;
  F = PF;
}
    
void
Grid :: setZeroValues(FloatMatrix* gridCoords)
{
  // gridCoords contains coordinates of points at which the solution is set to zero.
  // These coordinates must be expressed the grid coordinate system,
  // with the first coordinate (j) from 1 to n and the second (i) from 1 to m.
  // However, they do not need to be integers - the points do not need
  // to coincide with grid nodes.

  if (gridCoords->giveNumberOfColumns()!=2){
    OOFEM_ERROR("Matrix gridCoords passed to Grid :: setZeroValues should have 2 columns.");
  }
  int nValues = gridCoords->giveNumberOfRows();
  for ( int ival = 0; ival < nValues; ival++ ) {
    // get coordinates of the specified point
    double CPx = gridCoords->at(ival,1);
    double CPy = gridCoords->at(ival,2);
    // find nearest grid node (i, j)
    double j = (int) ceil(CPx - 0.5);
    if (j<1) j=1; else if (j>n) j=n;   
    double i = (int) ceil(i - 0.5);
    if (i<1) i=1; else if (i>m) i=m;
    // determine the index of the nearest grid node
    int CPInd = ij2ind(i, j, m);

    // Calculate time-distance of nearest grid node and freeze the value
    T [ CPInd ] = sqrt( (i - CPy)*(i - CPy) + (j - CPx)*(j - CPx) ) / F [ CPInd ];
    Frozen [ CPInd ] = true;

    // For all eight neighbours of the nearest grid node, do the same
        for ( int neigh = 0; neigh < 8; neigh++ ) {
            int ni = i + iOffsets_full [ neigh ];
            int nj = j + jOffsets_full [ neigh ];
            if ( isInDomain(ni, nj, m, n) ) {
                int nInd = ij2ind(ni, nj, m);
                double time = sqrt( (ni - CPy)*(ni - CPy) + (nj - CPx)*(nj - CPx) ) / F [ nInd ];
                if ( Frozen [ nInd ] ) {
                    T [ nInd ] = MIN(time, T [ nInd ]);
                } else   {
                    T [ nInd ] = time;
                    Frozen [ nInd ] = true;
                }
            }
        }
  }
  // after each external change, the solution becomes invalid
  solutionAvailable = false;
}

void
Grid :: setSolutionValueAt(int i, int j, double value)
{
  int ind = ij2ind(i, j, m);
  T[ind] = value;
  // after each external change, the solution becomes invalid
  solutionAvailable = false;
}
 
double
Grid :: giveSolutionValueAt(int i, int j)
{
  if (!solutionAvailable){
    int eFlag;
    // for the moment, the order is fixed as 2, because the second-order approach seems to be more accurate and not much slower
    fastMarch(2, eFlag);
    // TBD: one should check eFlag and send an error message if it is > 0
    solutionAvailable = true;
  }
  int ind = ij2ind(i, j, m);
  return T[ind];
}
 
void Grid :: printSolutionAsMatrix()
{
  if (!solutionAvailable){
    OOFEM_WARNING("Grid is printing a solution which has not been computed yet:");
  }
  for(int ii=0; ii<m; ii++) {
    printf("\n");
    for(int jj=0; jj<n; jj++)
      printf("%g  ", T[ii+jj*m]);
  }
  printf("\n");
}
    
void
Grid :: printSolutionAsData()
{
  for(int jj=0; jj<n; jj++){
    printf("\n");
    for(int ii=0; ii<m; ii++) {
      printf("%g %g %g\n", (double) jj+1, (double) ii+1, T[ii+jj*m]);
    }
  }
}
  
/// Fast marching method, solving the eikonal equation
void
Grid :: fastMarch(int order, int &eFlag)
{
    double time;
    int i,j,ni,nj,nInd;
    int tmpFlag = 0;
    eFlag = 0;

    // Create the initial narrow band
    if (!narrowBand) narrowBand = new Heap(m * n);
    // Loop over all grid points (not efficient, but done only once)
    for ( int ind = 0; ind < ( m * n ); ind++ ) {
        if ( Frozen [ ind ] ) {
            i = ind2i(ind, m);
            j = ind2j(ind, m);

            for ( int neigh = 0; neigh < 4; neigh++ ) {
                ni = i + iOffsets [ neigh ];
                nj = j + jOffsets [ neigh ];
                nInd = ij2ind(ni, nj, m);

                if ( isInDomain(ni, nj, m, n) && !Frozen [ nInd ] ) {
                    if ( !narrowBand->isInHeap(nInd) ) {
                        time = calcTime(ni, nj, F [ nInd ], order, tmpFlag);
                        narrowBand->insert(time, nInd);
                        if ( tmpFlag > eFlag ) {
                            eFlag = tmpFlag;
                        }
                    }
                }
            }
        }
    }

    // Loop
    int lCount = 0;
    int CPInd;
    while ( narrowBand->nElems() > 0 ) {
        lCount++;

        // Get minimal time-distance and index of this narrow-band element
        time = narrowBand->getSmallest(& CPInd);
        i = ind2i(CPInd, m);
        j = ind2j(CPInd, m);

        // Freeze and set time
        Frozen [ CPInd ] = true;
        T [ CPInd ] = time;

        // For all neighbours
        for ( int neigh = 0; neigh < 4; neigh++ ) {
            ni      = i + iOffsets [ neigh ];
            nj      = j + jOffsets [ neigh ];
            nInd = ij2ind(ni, nj, m);

            // If valid for consideration
            if ( isInDomain(ni, nj, m, n) && !Frozen [ nInd ] ) {
                time = calcTime(ni, nj, F [ nInd ], order, tmpFlag);
                // If T(ni,nj) has not been previously calculated
                if ( !narrowBand->isInHeap(nInd) ) {
                    narrowBand->insert(time, nInd);
                } else   {
                    narrowBand->update(time, nInd);
                }

                if ( tmpFlag > eFlag ) {
                    eFlag = tmpFlag;
                }
            }
        }
    }
}

// Time-dist calculation (for a grid with unit spacing)
double
Grid :: calcTime(int i, int j, double Fij, int order, int &eFlag) {
    // Get infinity.
    // NB: Is this good practice? As in the matlab code, using Inf as
    // a flag simplifies the code, so that's why I use this. Maybe this
    // doesn't work with some compilers?
    //double Inf = 1.0 / 0.0;
    double Inf = std::numeric_limits<float>::infinity();
    
    // Temporary error flag
    int tmpFlag = 0;

    // Frozen values at neighbors 
    double CrossVals [ 8 ];

    // time calculated
    double time;

    // Indices of neighbouring nodes
    int ni, nj, nInd;

    // Variables used in the final formula 
    double xmin1, xmin2, ymin1, ymin2;

    // Get values at surrounding nodes (at those that are frozen)
    for ( int iter = 0; iter < 8; iter++ ) {
        ni = i + icalcOffsets [ iter ];
        nj = j + jcalcOffsets [ iter ];
        nInd = ij2ind(ni, nj, m);
        if ( isInDomain(ni, nj, m, n) && Frozen [ nInd ] ) {
            CrossVals [ iter ] = T [ nInd ];
        } else   {
            CrossVals [ iter ] = Inf;
        }
    }

    // If none of the immediate neighbors is frozen, we cannot proceed
    if ( ( CrossVals [ 1 ] == Inf ) && ( CrossVals [ 2 ] == Inf ) &&
         ( CrossVals [ 5 ] == Inf ) && ( CrossVals [ 6 ] == Inf ) ) {
        eFlag = 0;
        time = Inf;
        return time;
    }

    // Calculate coefficients of quadratic equation
    double a = 0.;
    double b = 0.;
    double c = -1. / ( Fij * Fij );

    switch ( order ) {
    // First-order algorithm
    case 1:
        // Contribution of y-direction
        if ( !( ( CrossVals [ 1 ] == Inf ) && ( CrossVals [ 2 ] == Inf ) ) ) {
            ymin1 = MIN(CrossVals [ 1 ], CrossVals [ 2 ]);
            a += 1.;
            b -= 2. * ymin1;
            c += ymin1 * ymin1;
        }

        // Contribution of x-direction
        if ( !( ( CrossVals [ 5 ] == Inf ) && ( CrossVals [ 6 ] == Inf ) ) ) {
            xmin1 = MIN(CrossVals [ 5 ], CrossVals [ 6 ]);
            a += 1.;
            b -= 2. * xmin1;
            c += xmin1 * xmin1;
        }
        break;

    // Second-order algorithm
    case 2:
        // Contribution of y-direction
        if ( !( ( CrossVals [ 1 ] == Inf ) && ( CrossVals [ 2 ] == Inf ) ) ) {
            if ( CrossVals [ 1 ] < CrossVals [ 2 ] ) {
                ymin1 = CrossVals [ 1 ];
                if ( CrossVals [ 0 ] <= CrossVals [ 1 ] ) { //second-order formula can be applied
                    ymin2 = CrossVals [ 0 ];
                    a += 2.25;
                    b += -6. * ymin1 + 1.5 * ymin2;
                    c += 4.0 * ymin1 * ymin1 + 0.25 * ymin2 * ymin2 - 2.0 * ymin1 * ymin2;
                } else {      //first-order formula
                    a += 1.;
                    b -= 2. * ymin1;
                    c += ymin1 * ymin1;
                }
            } else {
                ymin1 = CrossVals [ 2 ];
                if ( CrossVals [ 3 ] <= CrossVals [ 2 ] ) { //second-order formula can be applied
                    ymin2 = CrossVals [ 3 ];
                    a += 2.25;
                    b += -6. * ymin1 + 1.5 * ymin2;
                    c += 4.0 * ymin1 * ymin1 + 0.25 * ymin2 * ymin2 - 2.0 * ymin1 * ymin2;
                } else {      //first-order formula
                    a += 1.;
                    b -= 2. * ymin1;
                    c += ymin1 * ymin1;
                }
            }
        }
        // Contribution of x-direction
        if ( !( ( CrossVals [ 5 ] == Inf ) && ( CrossVals [ 6 ] == Inf ) ) ) {
            if ( CrossVals [ 5 ] < CrossVals [ 6 ] ) {
                xmin1 = CrossVals [ 5 ];
                if ( CrossVals [ 4 ] <= CrossVals [ 5 ] ) { //second-order formula can be applied
                    xmin2 = CrossVals [ 4 ];
                    a += 2.25;
                    b += -6. * xmin1 + 1.5 * xmin2;
                    c += 4.0 * xmin1 * xmin1 + 0.25 * xmin2 * xmin2 - 2.0 * xmin1 * xmin2;
                } else {      //first-order formula
                    a += 1.;
                    b -= 2. * xmin1;
                    c += xmin1 * xmin1;
                }
            } else {
                xmin1 = CrossVals [ 6 ];
                if ( CrossVals [ 7 ] <= CrossVals [ 6 ] ) { //second-order formula can be applied
                    xmin2 = CrossVals [ 7 ];
                    a += 2.25;
                    b += -6. * xmin1 + 1.5 * xmin2;
                    c += 4.0 * xmin1 * xmin1 + 0.25 * xmin2 * xmin2 - 2.0 * xmin1 * xmin2;
                } else {      //first-order formula
                    a += 1.;
                    b -= 2. * xmin1;
                    c += xmin1 * xmin1;
                }
            }
        }
    }

    // Solve quadratic equation
    double d = ( b * b ) - ( 4.0 * a * c );

    // Error treatment and appropriate time-dist calculation
    if ( ( d < 0.0 ) && ( order == 2 ) ) {
      // if second-order did not work, try first-order formula
        time = calcTime(i, j, Fij, 1, tmpFlag);
        eFlag = MAX(1, tmpFlag);
    } else if ( ( d < 0.0 ) && ( order == 1 ) )             {
        if ( ( CrossVals [ 1 ] == Inf ) && ( CrossVals [ 2 ] == Inf ) ) {
            ymin1 = Inf;
        }
        if ( ( CrossVals [ 5 ] == Inf ) && ( CrossVals [ 6 ] == Inf ) ) {
            xmin1 = Inf;
        }

        eFlag = 2;
        time = MIN(xmin1, ymin1) + 1.0 / Fij;
    } else   {
        eFlag = 0;
        time = ( -b + sqrt(d) ) / ( 2.0 * a );
    }

    return time;
}

}
