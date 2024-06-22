#ifndef grid_h
#define grid_h

#include "floatmatrix.h"
#include "heap.h"

namespace oofem {
/**
 * Class that solves certain problems on a regular 2D grid, consisting of n x m nodes.
 * Currently it implements the fast marching algorithm that solves the eikonal equation.
 * This can be useful for instance for tracing evolving curves.
 * One could also implement finite differences or similar techniques.
 *
 * The indexing convention is a bit unusual but it follows the original code,
 * which was developed for image analysis. Subscripts "i" and "j" refer to rows
 * and columns of the matrix, but this means that "i" grows in the vertical
 * direction and corresponds to physical coordinate "y".
 * i ... from 1 to m, corresponds to coordinate y
 * j ... from 1 to n, corresponds to coordinate x
 * Actually, for the numerical method it does not matter which direction is horizontal
 * and which one is vertical, but the user should know how to interpret the data.
 * The constructor should get the horizontal dimension first and the vertical one second (n,m).
 * The coordinates of points at which zero value is prescribed should be passed as (x,y).
 * The prescribed field (e.g., propagating front velocity) should be stored by rows,
 * i.e., (1,1), (1,2), (1,3), ... (1,m), (2,1), (2,2,) ... (n,m).
 * If the solution is printed as a matrix, it is arranged on the screen in the natural way.
 * If the solution is printed as data, each output line corresponds to (x,y,value),
 * where x=j and y=i. Conversion to grids with non-unit spacing and shifter corner
 * must be done by the user.
 */

class Grid
{
public:
    /// Constructor (N = horizontal dimension, M = vertical dimension)
    Grid(int N, int M);
    /// Destructor
    ~Grid();

    /// Size information
    int giveSizeHorizontal() { return n; }
    int giveSizeVertical() { return m; }
    int giveSize(int dir);

    /// Set all flags to "unfrozen"
    void unFreeze();

    /// Set the details of the algorithm to be used
    void setMethod(int o, double i, int c) { order = o; initDiag = i; centDiff = c; }

    /// Method setting the values of the prescribed field (e.g., front velocities)
    //void setPrescribedField(double *PF);
    void setPrescribedFieldValue(int i, int j, double val) { F.at(i, j) = val; }
    FloatMatrix &givePrescribedField() { return F; }

    /// Methods setting the values of the unknown field at selected points (e.g., zero times)
    void setZeroValues(FloatMatrix *gridCoords);
    void setSolutionValueAt(int i, int j, double value);

    /// Output methods
    double giveSolutionValueAt(int i, int j);
    void printSolutionAsData();


private:
    /// Grid dimensions: number of grid nodes horizontally (n)  and vertically (m)
    int n, m;

    /// Algorithmic parameters
    int order;
    double initDiag;
    int centDiff;

    /// Flag indicating whether the solution has been computed
    bool solutionAvailable;

    /// Matrix storing the values of the unknown (computed) field (e.g., arrival times)
    FloatMatrix T;

    /// Matrix storing the values of the prescribed field (e.g., front velocities)
    FloatMatrix F;

    /// Array storing indicators of frozen points (at which the solution is already known)
    bool *Frozen;

    /// Heap used for efficient sorting and detection of smallest candidate
    Heap *narrowBand;

    /// Fast marching method, solving the eikonal equation
    void fastMarch(int &eFlag);

    /// Auxiliary function that evaluates the tentative value at one grid point, exploited by the fast marching method
    double calcTime(int i, int j, double Fij, int ord, int &eFlag);

    /// Utility methods
    int ij2ind(int _i, int _j, int _m) { return ( ( _i - 1 ) + ( _j - 1 ) * _m ); }
    int ind2i(int ind, int _m) { return 1 + ind - ( ind / _m ) * _m; }
    int ind2j(int ind, int _m) { return ( ( ind ) / _m ) + 1; }
    bool isInDomain(int i, int j, int _m, int _n) { return ( ( i >= 1 ) && ( j >= 1 ) && ( i <= _m ) && ( j <= _n ) ); }
};
} // end namespace oofem
#endif // grid_h
