/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef floatarray_h
#define floatarray_h

#include "oofemcfg.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#include <initializer_list>
#include <vector>
#include <iosfwd>

namespace oofem {
class IntArray;
class FloatMatrix;
class DataStream;

/**
 * Class representing vector of real numbers. This array can grow or shrink to
 * desired dimension. The lower value index of array is 1,
 * upper depends on array size.
 *
 * Tasks:
 * - Storing and returning a coefficient (method 'at') ;
 * - Expanding its size in order to store additional coefficients (method growTo )
 * - Performing basic operations : summation, product, rotation, etc.
 * - Assembling to itself another array, typically an elemental or nodal load vector (method 'assemble').
 * - Reading/writing its description on a given file.
 * - Introduced allocatedSize variable to allow dynamic resizing of array
 *   size possibly without memory reallocation. At startup array occupies space
 *   given by allocatedSpace = size. Then there can be
 *   - Further request for resizing array to smaller dimension
 *     then we only change size variable, but allocatedSize
 *     variable remain untouched - expecting possible array grow and then re-using
 *     previously allocated space.
 *   - If further request for growing then is necessary memory reallocation.
 *     This process is controlled in resize member function.
 *
 * Remarks:
 * - Method givePointer is an encapsulation crime. It is used only for
 *   speeding up method 'dot' of class RowColumn and for speeding method
 *   initialize.
 * 
 * @author Mikael Öhman
 * @author Erik Svenning
 * @author Jim Brouzoulis
 * @author many others (please add yourselves)
 */
class OOFEM_EXPORT FloatArray
{
protected:
    /// Stored values.
    std::vector< double > values;

public:
    /// @name Iterator for for-each loops:
    //@{
    std::vector< double > :: iterator begin() { return this->values.begin(); }
    std::vector< double > :: iterator end() { return this->values.end(); }
    std::vector< double > :: const_iterator begin() const { return this->values.begin(); }
    std::vector< double > :: const_iterator end() const { return this->values.end(); }
    //@}

    /// Constructor for sized array. Data is zeroed.
    FloatArray(int n = 0) : values(n) { }
    /// Disallow double parameter, which can otherwise give unexpected results.
    FloatArray(double)  = delete;
    /// Copy constructor. Creates the array from another array.
    FloatArray(const FloatArray &src) : values(src.values) { }
    /// Move constructor. Creates the array from another array.
    FloatArray(FloatArray &&src) : values(std::move(src.values)) { }
    /// Initializer list constructor.
    inline FloatArray(std :: initializer_list< double >list) : values(list) { }
    /// Destructor.
    virtual ~FloatArray() {};

    /// Assignment operator
    FloatArray &operator = (const FloatArray &src) { values = src.values; return *this; }
    /// Move operator
    FloatArray &operator = (FloatArray &&src) { values = std::move(src.values); return *this; }
    /// Assignment operator.
    inline FloatArray &operator = (std :: initializer_list< double >list) { values = list; return *this; }

    /// Add one element
    void push_back(const double &iVal) {values.push_back(iVal);}

    /// Returns true if no element is NAN or infinite
    bool isFinite() const;

    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Provides 1-based indexing access.
     * @param i Position of coefficient in array.
     */
#ifdef DEBUG
    double &at(int i);
#else
    inline double &at(int i) { return values [ i - 1 ]; }
#endif
    /**
     * Coefficient access function. Returns l-value of coefficient at given
     * position of the receiver. Provides 1-based indexing access.
     * @param i Position of coefficient in array.
     */
#ifdef DEBUG
    double at(int i) const;
#else
    inline double at(int i) const { return values [ i - 1 ]; }
#endif

    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Provides 0-based indexing access.
     * @param i Position of coefficient in array.
     */
#ifdef DEBUG
    double &operator() (int i);
    double &operator[] (int i);
#else
    inline double &operator() (int i) {
        return values [ i ];
    };
    inline double &operator[] (int i) {
        return values [ i ];
    };
#endif
    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Provides 0-based indexing access.
     * @param i Position of coefficient in array.
     */
#ifdef DEBUG
    const double &operator() (int i) const;
    const double &operator[] (int i) const;
#else
    inline const double &operator() (int i) const { return values [ i ]; }
    inline const double &operator[] (int i) const { return values [ i ]; }
#endif
    /** Checks size of receiver towards requested bounds.
     * Current implementation will call exit(1), if dimension
     * mismatch found.
     * @param i Required size of receiver.
     */
    void checkBounds(int i) const;
    /**
     * Checks size of receiver towards values stored in loc array.
     * Expands the receiver if loc points to coefficients beyond the size of receiver.
     * @param loc Array with indices.
     */
    void checkSizeTowards(const IntArray &loc);
    /**
     * Allocates enough size to fit s, and clears the array.
     * @param s New reserved size.
     */
    void reserve(int s);
    /**
     * Checks size of receiver towards requested bounds.
     * If dimension mismatch, size is adjusted accordingly.
     * Old values are copied over and new space is zeroed.
     * @param s New size.
     * @param allocChunk Additional space to allocate.
     */
    void resizeWithValues(int s, int allocChunk = 0);
    /**
     * Resizes receiver towards requested size. Array is zeroed.
     * @param s New size.
     */
    void resize(int s);
    /**
     * Clears receiver (zero size).
     * Same effect as resizing to zero, but has a clearer meaning and intent when used.
     */
    void clear() { this->values.clear(); }
    /**
     * Resizes the size of the receiver to requested bounds. Memory allocation always happens, more preferably use
     * resize() function instead. Array is zeroed.
     * @param s New size.
     */
    void hardResize(int s);
    /**
     * Returns nonzero if all coefficients of the receiver are 0, else returns zero.
     */
    bool containsOnlyZeroes() const;
    /// Returns the size of receiver.
    int giveSize() const { return (int)this->values.size(); }
    /// Returns true if receiver is not empty.
    bool isNotEmpty() const { return !this->values.empty(); }
    /// Returns true if receiver is empty.
    bool isEmpty() const { return this->values.empty(); }
    /**
     * Switches the sign of every coefficient of receiver.
     * @return receiver.
     */
    void negated();
    /**
     * Print receiver on stdout. Useful for debugging.
     */
    virtual void printYourself() const;
    /**
     * Print receiver on stdout with custom name.
     * @param name Display name of reciever.
     */    
    virtual void printYourself(const std::string &name) const;
    /**
     * Print receiver on stdout with high accuracy.
     */
    virtual void pY() const;
    /// Zeroes all coefficients of receiver.
    void zero();
    /**
     * Appends array to reciever.
     * @param a Values to be appended.
     */
    void append(const FloatArray &a);
    /**
     * Appends value to reciever.
     * @param a Value to be appended.
     */
    void append(double a);
    /**
     * Receiver becomes the result of the product of aMatrix and anArray.
     * Adjusts the size of receiver if necessary.
     */
    void beProductOf(const FloatMatrix &aMatrix, const FloatArray &anArray);
    /**
     * Receiver becomes the result of the product of aMatrix^T and anArray.
     * Adjusts the size of receiver if necessary.
     */
    void beTProductOf(const FloatMatrix &aMatrix, const FloatArray &anArray);
    /**
     * Computes vector product (or cross product) of vectors given as parameters,
     * @f$ v_1 \times v_2 @f$, and stores the result into receiver.
     * @param v1 First vector in the product.
     * @param v2 Second vector in the product.
     */
    void beVectorProductOf(const FloatArray &v1, const FloatArray &v2);
    /**
     * Sets receiver to be @f$ a = s b @f$.
     * @param s Scaling factor.
     * @param b Vector to be scaled.
     */
    void beScaled(double s, const FloatArray &b);
    /**
     * Adds array src to receiver. If the receiver's size is zero, it adjusts its size
     * to size of src array. If receiver's size is nonzero and different from src
     * array size an error is generated.
     * @param src Array to add to receiver.
     */
    void add(const FloatArray &src);
    /**
     * Adds array times factor to receiver.
     * If the receiver's size is zero, it adjusts its size
     * to size of b array. If receiver's size is nonzero and different from b
     * array size an error is generated.
     * @param factor Factor to be multiplied with b.
     * @param b Will be multiplied by factor and added to receiver.
     */
    void add(double factor, const FloatArray &b);
    /**
     * Adds the product @f$ b^T . s dV @f$.
     * If the receiver's size is zero, it adjusts its size.
     * @param b Matrix b in the equation.
     * @param s Array s in the equation.
     * @param dV Scalar dV in the equation.
     */
    void plusProduct(const FloatMatrix &b, const FloatArray &s, double dV);
    /**
     * Adds scalar to receiver.
     * @param offset Scalar to add
     */
    void add(double offset);
    /**
     * Subtracts array src to receiver. If the receiver's size is zero, it adjusts its size
     * to size of src array. If recever's size is nonzero and different from src
     * array size an error is generated.
     */
    void subtract(const FloatArray &src);
    /**
     * Multiplies receiver with scalar.
     * @param s Scalar to multiply by.
     */
    void times(double s);
    /**
     * Sets receiver to maximum of a or b's respective elements.
     * @param a Array of size n.
     * @param b Array of size n.
     */
    void beMaxOf(const FloatArray &a, const FloatArray &b);
    /**
     * Sets receiver to be minimum of a or b's respective elements.
     * @param a Array of size n.
     * @param b Array of size n.
     */
    void beMinOf(const FloatArray &a, const FloatArray &b);
    /**
     * Sets receiver to be a - b.
     * @param a Array which receiver goes to.
     * @param b Array which receiver comes from.
     */
    void beDifferenceOf(const FloatArray &a, const FloatArray &b);
    /**
     * Sets receiver to be a - b, using only the first n entries.
     * @param a Array which receiver goes to.
     * @param b Array which receiver comes from.
     * @param n Only first n entries are taken.
     */
    void beDifferenceOf(const FloatArray &a, const FloatArray &b, int n);
    /**
     * Extract sub vector form src array and stores the result into receiver.
     * @param src source vector for sub vector
     * @param indx Determines sub vector. Receiver size will be indx max value,
     * and on i-th position of subVector will be src(indx->at(i)) value.
     */
    void beSubArrayOf(const FloatArray &src, const IntArray &indx);
    /**
     * Adds the given vector as sub-vector to receiver. The sub-vector values will be added to receivers
     * corresponding receiver's values at positions (si,...,si+src.size). The size of receiver will be
     * adjusted, if necessary.
     * @param src Sub-vector to be added.
     * @param si Determines the position (receiver's 1-based index) of first src value to be added.
     */
    void addSubVector(const FloatArray &src, int si);
    /**
     * Assembles the array fe (typically, the load vector of a finite
     * element) into the receiver, using loc as location array.
     * @param fe Array to be assembled.
     * @param loc Array of code numbers. src(i) value will
     * be added to receiver value at position loc(i)
     * (if this loc(i) value is nonzero).
     */
    void assemble(const FloatArray &fe, const IntArray &loc);
    /**
     * Assembles the array fe with each component squared.
     * @param fe Array to be assembled (with each component squared)
     * @param loc Location array.
     */
    void assembleSquared(const FloatArray &fe, const IntArray &loc);
    /**
     * Copy the given vector as sub-vector to receiver. The sub-vector values will be set to receivers
     * values starting at at positions (si,...,si+src.size). The size of receiver will be
     * adjusted, if necessary.
     * @param src Sub-vector to be added
     * @param si Determines the position (receiver's 1-based index) of first src value to be added.
     */
    void copySubVector(const FloatArray &src, int si);
    /**
     * Computes the distance between position represented by receiver and position given as parameter.
     * @param x Coordinate to calculate distance from.
     */
    double distance(const FloatArray &x) const;
    /// @see distance
    double distance(const FloatArray *x) const { return this->distance(* x); }

    /**
     * Computes distance between the position represented by the reciever and a line segment represented by it's start
     * point iP1 and it's end point iP2.
     * The local coordinate oXi is in the range [0,1]
     * @author Erik Svenning, August 2013.
     */
    double distance(const FloatArray &iP1, const FloatArray &iP2, double &oXi, double &oXiUnbounded) const;
    double distance_square(const FloatArray &iP1, const FloatArray &iP2, double &oXi, double &oXiUnbounded) const;

    /**
     * Computes the square of distance between position represented by receiver and position given as parameter.
     * @param x Coordinate to calculate squared distance from.
     */
    double distance_square(const FloatArray &x) const;

    /**Returns index (between 1 and Size) of minimum element in the array
     */
    int giveIndexMinElem(void);

    /**Returns index (between 1 and Size) of maximum element in the array
     */
    int giveIndexMaxElem(void);

    /**
     * Computes the dot product (or inner product) of receiver and argument.
     * @param x Vector to contract to receiver.
     */
    double dotProduct(const FloatArray &x) const;

    /**
     * Computes the dot product (or inner product) of receiver and argument.
     * @param x Vector to contract to receiver.
     * @param size Number of elements to contract. May not be larger than
     */
    double dotProduct(const FloatArray &x, int size) const;

    /**
     * Normalizes receiver. Euclidean norm is used, after operation receiver
     * will have this norm equal to 1.0.
     * @return modified receiver
     */
    double normalize();
    /**
     * Computes the norm (or length) of the vector.
     * @return The Euclidean norm of the vector.
     */
    double computeNorm() const;

    /**
     * Computes the square of the norm.
     * @return Squared norm.
     */
    double computeSquaredNorm() const;
    /**
     * Computes the sum of receiver values.
     * @return Sum of receiver.
     */
    double sum() const;
    /**
     * Computes the product of receiver values.
     * @return Product of receiver.
     */
    double product() const;

    /**
     * Returns the receiver a rotated according the change-of-base matrix r.
     * @param r Rotation matrix.
     * @param mode If mode == 't' the method performs the operation  @f$ a = t^{\mathrm{T}} \cdot r @f$,
     * else if mode = 'n' then the method performs the operation  @f$ a = t \cdot r @f$.
     * @return modified receiver.
     */
    void rotatedWith(FloatMatrix &r, char mode);
    /**
     * Gives the pointer to the raw data, breaking encapsulation.
     * @return Pointer to values of array
     */
    inline const double *givePointer() const { return values.data(); }
    inline double *givePointer() { return values.data(); }
    /**
     * Reciever will be a vector with 9 components formed from a 3x3 matrix.
     * Order of matrix components in vector: 11, 22, 33, 23, 13, 12, 32, 31, 21
     * @param aMatrix Matrix to transform.
     */
    void beVectorForm(const FloatMatrix &aMatrix);
    /**
     * Reciever will be a vector with 6 components formed from a 3x3 matrix.
     * Off-diagonals of the matrix are symmetrized.
     * Order of matrix components in vector: 11, 22, 33, 23, 13, 12
     * @param aMatrix Matrix to transform.
     */
    void beSymVectorForm(const FloatMatrix &aMatrix);
    void beSymVectorFormOfStrain(const FloatMatrix &aMatrix);

    /**
     * Swaps the fourth and sixth index in the array. This is to reorder the indices
     * from OOFEM's order to Abaqus' and vice versa.
     */
    void changeComponentOrder();


    /**
     * Reciever will be set to a given column in a matrix
     */
    void beColumnOf(const FloatMatrix &mat, int col);

    contextIOResultType storeYourself(DataStream &stream) const;
    contextIOResultType restoreYourself(DataStream &stream);
    int givePackSize(DataStream &buff) const;

    friend std :: ostream &operator << ( std :: ostream & out, const FloatArray & x );

    friend class FloatMatrix;

    ///@name IML compatibility
    //@{
    /// Assignment of scalar to all components of receiver
    FloatArray &operator = ( const double & );
    //@}

#ifdef BOOST_PYTHON
    void __setitem__(int i, double val) { this->values[i] = val; }
    double __getitem__(int i) { return this->values[i]; }
    void beCopyOf(const FloatArray &src) { this->operator = ( src ); }
#endif
};

const FloatArray ZeroVector = {0.0,0.0,0.0};

///@name IML compatibility
//@{
/// Vector multiplication by scalar
FloatArray &operator *= ( FloatArray & x, const double & a );
FloatArray operator *( const double & a, const FloatArray & x );
FloatArray operator *( const FloatArray & x, const double & a );
FloatArray operator + ( const FloatArray & x, const FloatArray & y );
FloatArray operator - ( const FloatArray & x, const FloatArray & y );
FloatArray &operator += ( FloatArray & x, const FloatArray & y );
FloatArray &operator -= ( FloatArray & x, const FloatArray & y );

double norm(const FloatArray &x);
double dot(const FloatArray &x, const FloatArray &y);
//@}
} // end namespace oofem
#endif // floatarray_h
