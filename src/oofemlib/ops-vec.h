#pragma once

#ifndef __MATH_OPS
    #error __MATH_OPS must be defined to use this header (only oofem internal headers should include it)
#endif

#include<string>
#include<initializer_list>

#ifndef _MATHOPS_STAGE
    #define _MATHOPS_STAGE 100
#endif


namespace oofem{
    #ifndef _USE_EIGEN
        class FloatMatrix;
        class FloatArray;
    #endif
    class IntArray;
}

namespace oofem::vec{

    /***** FloatArray *****/

    /**
     * Returns the receiver a rotated according the change-of-base matrix r.
     * @param r Rotation matrix.
     * @param mode If mode == 't' the method performs the operation  @f$ a = t^{\mathrm{T}} \cdot r @f$,
     * else if mode = 'n' then the method performs the operation  @f$ a = t \cdot r @f$.
     * @return modified receiver.
     */
    void rotatedWith(FloatArray& a, FloatMatrix &r, char mode);

    /**
     * Print receiver on stdout. Useful for debugging.
     */
    void printYourself(const FloatArray& a);
    /**
     * Print receiver on stdout with custom name.
     * @param name Display name of reciever.
     */    
    void printYourself(const FloatArray& a, const std::string &name);
    /**
     * Print receiver to file
     * @param filename Name of recieving file.
     * @param showDimensions Determins if dimesions should be included in output
     */
    void printYourselfToFile(const FloatArray& a, const std::string filename, const bool showDimensions=true);
    /**
     * Print receiver on stdout with high accuracy.
     */
    void pY(const FloatArray& a);

    /**
     * Computes distance between the position represented by the reciever and a line segment represented by it's start
     * point iP1 and it's end point iP2.
     * The local coordinate oXi is in the range [0,1]
     * @author Erik Svenning, August 2013.
     */
    double distance(const FloatArray& pos, const FloatArray &iP1, const FloatArray &iP2, double &oXi, double &oXiUnbounded);
    double distance_square(const FloatArray& pos, const FloatArray &iP1, const FloatArray &iP2, double &oXi, double &oXiUnbounded);

    /**
     * Reciever will be a vector with 9 components formed from a 3x3 matrix.
     * Order of matrix components in vector: 11, 22, 33, 23, 13, 12, 32, 31, 21
     * @param aMatrix Matrix to transform.
     */
    FloatArray beVectorForm(const FloatMatrix &aMatrix);
    /**
     * Reciever will be a vector with 6 components formed from a 3x3 matrix.
     * Off-diagonals of the matrix are symmetrized.
     * Order of matrix components in vector: 11, 22, 33, 23, 13, 12
     * @param aMatrix Matrix to transform.
     */
    FloatArray beSymVectorForm(const FloatMatrix &aMatrix);
    FloatArray beSymVectorFormOfStrain(const FloatMatrix &aMatrix);

    /**
     * Swaps the fourth and sixth index in the array. This is to reorder the indices
     * from OOFEM's order to Abaqus' and vice versa.
     */
    void changeComponentOrder(FloatArray& a);


    /**
     * Assembles the array fe (typically, the load vector of a finite
     * element) into the receiver, using loc as location array.
     * @param fe Array to be assembled.
     * @param loc Array of code numbers. src(i) value will
     * be added to receiver value at position loc(i)
     * (if this loc(i) value is nonzero).
     */
    void assemble(FloatArray& dst, const FloatArray &fe, const IntArray &loc);
    /**
     * Assembles the array fe with each component squared.
     * @param fe Array to be assembled (with each component squared)
     * @param loc Location array.
     */
    void assembleSquared(FloatArray& dst, const FloatArray &fe, const IntArray &loc);



   /**** free functions from floatarray.h *****/
   double norm(const FloatArray &x);
   double norm_square(const FloatArray &x);
   double dot(const FloatArray &x, const FloatArray &y);
   double distance(const FloatArray &x, const FloatArray &y);
   double distance_square(const FloatArray &x, const FloatArray &y);
   bool isfinite(const FloatArray &x);
   bool iszero(const FloatArray &x);
   double sum(const FloatArray & x);
   double product(const FloatArray & x);

} /* namespace oofem::vec */



#if 0
namespace oofem{
    using vec::norm;
    using vec::norm_square;
    using vec::dot;
    using vec::distance;
    using vec::distance_square;
    using vec::isfinite;
    using vec::iszero;
    using vec::sum;
    using vec::product;
}; /* namespace oofem */
#endif
