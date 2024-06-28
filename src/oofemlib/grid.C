#include "grid.h"
#include "error.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

namespace oofem {

// NB: When using matrix indices (i,j), they run from 1 to (m,n).
// When using 1-dim array index 'ind', it runs from 0 to (m*n-1).

// Index offsets of all neighbors (direct and diagonal ones)
int iOffsets_full[] = {
    -1, -1, -1,  0,  0,  1,  1,  1
};
int jOffsets_full[] = {
    -1,  0,  1, -1,  1, -1,  0,  1
};

// Indicator which neighbors are diagonal ones
bool is_diag[] = {
    true, false, true, false, false, true, false, true
};

// Index offsets of direct neighbours
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


Grid :: Grid(int N, int M) :
    n(N), m(M),
    T(m, n),
    F(m, n)
{
    solutionAvailable = false;

    // (calloc sets all bits of the allocated memory to zero, so the initial values are "false")
    Frozen = ( bool * ) calloc( m * n, sizeof( bool ) );
    narrowBand = new Heap(m *n);
}

Grid :: ~Grid()
{
    if ( Frozen ) {
        free(Frozen);
    }
    delete narrowBand;
}

int
Grid :: giveSize(int dir)
{
    switch ( dir ) {
    case 1: return n;

    case 2: return m;

    default: return 0;
    }
}

/*
 * void
 * Grid :: setPrescribedField(double *PF)
 * {
 * // For efficiency, F is passed by its pointer and from now on
 * // will be taken care of (and deleted) by the Grid.
 * if ( F ) {
 *    delete F;
 * }
 * F = PF;
 * }
 */

void
Grid :: unFreeze()
{
    for ( int ind = 0; ind < ( m * n ); ind++ ) {
        Frozen [ ind ] = false;
    }
    narrowBand->setToEmpty(m * n);
    solutionAvailable = false;
}

void
Grid :: setZeroValues(FloatMatrix *gridCoords)
{
    // gridCoords contains coordinates of points at which the solution is set to zero.
    // These coordinates must be expressed the grid coordinate system,
    // with the first coordinate (j) from 1 to n and the second (i) from 1 to m.
    // However, they do not need to be integers - the points do not need
    // to coincide with grid nodes.

    if ( gridCoords->giveNumberOfColumns() != 2 ) {
        OOFEM_ERROR("Matrix gridCoords passed to Grid :: setZeroValues should have 2 columns.");
    }
    int nValues = gridCoords->giveNumberOfRows();
    for ( int ival = 1; ival <= nValues; ival++ ) {
        // get coordinates of the specified point
        double CPx = gridCoords->at(ival, 1);
        double CPy = gridCoords->at(ival, 2);
        // find nearest grid node (i, j)
        double j = ( int ) ceil(CPx - 0.5);
        if ( j < 1 ) {
            j = 1;
        } else if ( j > n ) {
            j = n;
        }
        double i = ( int ) ceil(CPy - 0.5);
        if ( i < 1 ) {
            i = 1;
        } else if ( i > m ) {
            i = m;
        }
        // determine the index of the nearest grid node
        int CPInd = ij2ind((int)i, (int)j, (int)m);

        // Calculate time-distance of nearest grid node and freeze the value
        double Fij = F.at((size_t)i, (size_t)j);
        T.at((size_t)i, (size_t)j) = sqrt( ( i - CPy ) * ( i - CPy ) + ( j - CPx ) * ( j - CPx ) ) / Fij;
        Frozen [ CPInd ] = true;

        // For four direct neighbors or all eight neighbours of the nearest grid node, do the same
        // (depending on parameter initDiag)

        if ( initDiag <= 0. ) { // initialize direct neighbors only
            for ( int neigh = 0; neigh < 4; neigh++ ) {
                int ni = (int)(i + iOffsets [ neigh ]);
                int nj = (int)(j + jOffsets [ neigh ]);
                if ( isInDomain(ni, nj, m, n) ) {
                    int nInd = ij2ind(ni, nj, m);
                    double time;
                    if ( centDiff > 0 ) { // use the average speed
                        time = sqrt( ( ni - CPy ) * ( ni - CPy ) + ( nj - CPx ) * ( nj - CPx ) ) / ( 0.5 * ( Fij + F.at(ni, nj) ) );
                    } else { // use the speed at the arrival point
                        time = sqrt( ( ni - CPy ) * ( ni - CPy ) + ( nj - CPx ) * ( nj - CPx ) ) / F.at(ni, nj);
                    }
                    if ( Frozen [ nInd ] ) {
                        T.at(ni, nj) = std::min( time, T.at(ni, nj) );
                    } else {
                        T.at(ni, nj) = time;
                        Frozen [ nInd ] = true;
                    }
                }
            }
        } else { // initialize all neighbors
            for ( int neigh = 0; neigh < 8; neigh++ ) {
                int ni = (int) (i + iOffsets_full [ neigh ]);
                int nj = (int) (j + jOffsets_full [ neigh ]);
                if ( isInDomain(ni, nj, m, n) ) {
                    int nInd = ij2ind(ni, nj, m);
                    double time;
                    if ( centDiff > 0 ) { // use the average speed
                        time = sqrt( ( ni - CPy ) * ( ni - CPy ) + ( nj - CPx ) * ( nj - CPx ) ) / ( 0.5 * ( Fij + F.at(ni, nj) ) );
                    } else { // use the speed at the arrival point
                        time = sqrt( ( ni - CPy ) * ( ni - CPy ) + ( nj - CPx ) * ( nj - CPx ) ) / F.at(ni, nj);
                    }
                    // for diagonal neighbors, use initDiag as a scaling factor
                    if ( is_diag [ neigh ] ) {
                        time *= initDiag;
                    }
                    if ( Frozen [ nInd ] ) {
                        T.at(ni, nj) = std::min( time, T.at(ni, nj) );
                    } else {
                        T.at(ni, nj) = time;
                        Frozen [ nInd ] = true;
                    }
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
    T.at(i, j) = value;
    // after each external change, the solution becomes invalid
    solutionAvailable = false;
}

double
Grid :: giveSolutionValueAt(int i, int j)
{
    if ( !solutionAvailable ) {
        int eFlag;
        fastMarch(eFlag);
        // TBD: one should check eFlag and send an error message if it is > 0
        solutionAvailable = true;
    }
    return T.at(i, j);
}

void
Grid :: printSolutionAsData()
{
    for ( int j = 1; j <= n; j++ ) {
        printf("\n");
        for ( int i = 1; i <= m; i++ ) {
            printf( "%d %d %g\n", j, i, T.at(i, j) );
        }
    }
}

/// Fast marching method, solving the eikonal equation
void
Grid :: fastMarch(int &eFlag)
{
    eFlag = 0;

    // Create the initial narrow band
    if ( !narrowBand ) {
        narrowBand = new Heap(m *n);
    }

    // Loop over all grid points (not efficient, but done only once)
    for ( int ind = 0; ind < ( m * n ); ind++ ) {
        if ( Frozen [ ind ] ) {
            int i = ind2i(ind, m);
            int j = ind2j(ind, m);

            for ( int neigh = 0; neigh < 4; neigh++ ) {
                int ni = i + iOffsets [ neigh ];
                int nj = j + jOffsets [ neigh ];
                int nInd = ij2ind(ni, nj, m);

                if ( isInDomain(ni, nj, m, n) && !Frozen [ nInd ] ) {
                    if ( !narrowBand->isInHeap(nInd) ) {
                        int tmpFlag = 0;
                        double time = calcTime(ni, nj, F.at(ni, nj), order, tmpFlag);
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
    int CPInd;
    while ( narrowBand->nElems() > 0 ) {
        // Get minimal time-distance and index of this narrow-band element
        double time = narrowBand->getSmallest(& CPInd);
        int i = ind2i(CPInd, m);
        int j = ind2j(CPInd, m);

        // Freeze and set time
        Frozen [ CPInd ] = true;
        T.at(i, j) = time;

        // For all neighbours
        for ( int neigh = 0; neigh < 4; neigh++ ) {
            int ni = i + iOffsets [ neigh ];
            int nj = j + jOffsets [ neigh ];
            int nInd = ij2ind(ni, nj, m);

            // If valid for consideration
            if ( isInDomain(ni, nj, m, n) && !Frozen [ nInd ] ) {
                int tmpFlag;
                time = calcTime(ni, nj, F.at(ni, nj), order, tmpFlag);
                // If T(ni,nj) has not been previously calculated
                if ( !narrowBand->isInHeap(nInd) ) {
                    narrowBand->insert(time, nInd);
                } else {
                    narrowBand->update(time, nInd);
                }

                if ( tmpFlag > eFlag ) {
                    eFlag = tmpFlag;
                }
            }
        }
    }
    solutionAvailable = true;
}

// Time-dist calculation (for a grid with unit spacing)
double
Grid :: calcTime(int i, int j, double Fij, int ord, int &eFlag)
{
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
    //mj
    double Fx = -1.;
    double Fy = -1.;

    // time calculated
    double time;

    // Variables used in the final formula
    double xmin1 = 0., xmin2, ymin1 = 0., ymin2;

    // Get values at surrounding nodes (at those that are frozen)
    for ( int iter = 0; iter < 8; iter++ ) {
        int ni = i + icalcOffsets [ iter ];
        int nj = j + jcalcOffsets [ iter ];
        int nInd = ij2ind(ni, nj, m);
        if ( isInDomain(ni, nj, m, n) && Frozen [ nInd ] ) {
            CrossVals [ iter ] = T.at(ni, nj);
        } else {
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
    // mj double c = -1. / ( Fij * Fij );
    double c = 0.;

    switch ( ord ) {
    // First-order algorithm
    case 1:
        // Contribution of y-direction
        if ( !( ( CrossVals [ 1 ] == Inf ) && ( CrossVals [ 2 ] == Inf ) ) ) {
            ymin1 = std::min(CrossVals [ 1 ], CrossVals [ 2 ]);
            a += 1.;
            b -= 2. * ymin1;
            c += ymin1 * ymin1;
            //mj
            int ni, nj;
            if ( CrossVals [ 1 ] < CrossVals [ 2 ] ) {
                ni = i + icalcOffsets [ 1 ];
                nj = j + jcalcOffsets [ 1 ];
            } else {
                ni = i + icalcOffsets [ 2 ];
                nj = j + jcalcOffsets [ 2 ];
            }
            Fy = F.at(ni, nj);
            //end mj
        }

        // Contribution of x-direction
        if ( !( ( CrossVals [ 5 ] == Inf ) && ( CrossVals [ 6 ] == Inf ) ) ) {
            xmin1 = std::min(CrossVals [ 5 ], CrossVals [ 6 ]);
            a += 1.;
            b -= 2. * xmin1;
            c += xmin1 * xmin1;
            //mj
            int ni, nj;
            if ( CrossVals [ 5 ] < CrossVals [ 6 ] ) {
                ni = i + icalcOffsets [ 5 ];
                nj = j + jcalcOffsets [ 5 ];
            } else {
                ni = i + icalcOffsets [ 6 ];
                nj = j + jcalcOffsets [ 6 ];
            }
            Fx = F.at(ni, nj);
            //end mj
        }
        break;

    // Second-order algorithm
    case 2:
    default:
        // Contribution of y-direction
        if ( !( ( CrossVals [ 1 ] == Inf ) && ( CrossVals [ 2 ] == Inf ) ) ) {
            if ( CrossVals [ 1 ] < CrossVals [ 2 ] ) {
                ymin1 = CrossVals [ 1 ];
                //mj
                int ni = i + icalcOffsets [ 1 ];
                int nj = j + jcalcOffsets [ 1 ];
                Fy = F.at(ni, nj);
                //end mj
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
                //mj
                int ni = i + icalcOffsets [ 2 ];
                int nj = j + jcalcOffsets [ 2 ];
                Fy = F.at(ni, nj);
                //end mj
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
                //mj
                int ni = i + icalcOffsets [ 5 ];
                int nj = j + jcalcOffsets [ 5 ];
                Fx = F.at(ni, nj);
                //end mj
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
                //mj
                int ni = i + icalcOffsets [ 6 ];
                int nj = j + jcalcOffsets [ 6 ];
                Fx = F.at(ni, nj);
                //end mj
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

    // for centDiff=2, the average speed is used, otherwise the speed at arrival point is used
    double Faver;
    if ( centDiff != 2 ) {
        Faver = Fij;
    } else if ( Fx >= 0. && Fy >= 0. ) {
        Faver = 0.5 * Fij + 0.25 * ( Fx + Fy );
    } else if ( Fx >= 0. ) {
        Faver = 0.5 * ( Fij + Fx );
    } else if ( Fy >= 0. ) {
        Faver = 0.5 * ( Fij + Fy );
    } else {
        Faver = Fij;
    }

    c -= 1. / ( Faver * Faver );

    // Solve quadratic equation
    double d = ( b * b ) - ( 4.0 * a * c );

    // Error treatment and appropriate time-dist calculation
    if ( ( d < 0.0 ) && ( ord == 2 ) ) {
        // if second-order method did not work, try first-order formula
        time = calcTime(i, j, Fij, 1, tmpFlag);
        eFlag = std::max(1, tmpFlag);
    } else if ( ( d < 0.0 ) && ( ord == 1 ) ) {
        if ( ( CrossVals [ 1 ] == Inf ) && ( CrossVals [ 2 ] == Inf ) ) {
            ymin1 = Inf;
        }
        if ( ( CrossVals [ 5 ] == Inf ) && ( CrossVals [ 6 ] == Inf ) ) {
            xmin1 = Inf;
        }

        eFlag = 2;
        time = std::min(xmin1, ymin1) + 1.0 / Fij;
    } else {
        eFlag = 0;
        time = ( -b + sqrt(d) ) / ( 2.0 * a );
    }

    return time;
}
}
