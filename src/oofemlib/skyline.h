/* $Header: /home/cvs/bp/oofem/oofemlib/src/skyline.h,v 1.11 2003/04/06 14:08:25 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


//   *************************************************************************
//   *** CLASS SPARSE MATRIX SKYLINE CLASS                        ************
//   *************************************************************************

#ifndef skyline_h
#define skyline_h

#include "sparsemtrx.h"

/**
 * Class implementing sparse matrix stored in skyline form. This class
 * assumes symmetric form of matrix to be stored, i.e., only upper half
 * part is therefore stored. Coefficients are stored in one dimensional
 * float array, containing for each column the values from the diagonal 
 * to the last nonzero entry. 
 * This requires to remember the adresses of diagonal members in such
 * array (addr attribute)
 * @see SparseMtrx class
 */
class Skyline : public SparseMtrx
{
    /*
     * This class implements a symmetric matrix stored in a compact (skyline)
     * form. A skyline is usually an attribute of the linear system.
     * DESCRIPTION :
     *
     * Attribute 'isFac-torized' is True if the skyline is already in U(T).D.U
     *         factorized form, else it is False.
     * Attribute adr  is array of diagonal members, it's size is size+1 (adr[0]=neni, adr[1]=1)
     * Attribute mtrx is double pointer to skyline stored in a array form
     *         (but we start from index 1)
     * TASKS :
     * - building its internal storage structure (method 'buildInternalStructure')
     * - store and localize local mtrices (method 'localize')
     * - performing standard operations : multiplication by array (method 'times')
     * - possible factorization and backSubstitution (recognized by nonzero result of
     *  canBeFactorized) (methods 'factorize' and 'bacSobstitution')
     * - setting all coefficients to zero (method 'zero')
     */

protected:
    /// Total number of nonzero coefficients stored.
    int nwk;
    /// Integer array holding adresses of diagonal memebers.
    IntArray *adr;
    /// Vector of stored coefficients.
    double *mtrx;
    /// Flag indicating whether factorized.
    int isFactorized;


public:

    /** Constructor. Before any operation an internal profile must be built.
     * @see builInternalStructure
     */
    Skyline(int n);
    /** Constructor. Before any operation an internal profile must be built.
     * @see builInternalStructure
     */
    Skyline();
    /// Destructor
    ~Skyline();

    /** Returns {\bf newly allocated} copy of receiver. Programmer must take
     * care about proper deallocation of allocated space.
     * @return newly allocated copy of receiver */
    SparseMtrx *GiveCopy() const;

    /** Evaluates a product of receiver with vector.
     * @param x array to be multiplied with receiver
     * @param answer result of product of receiver and x parameter
     */
    void times(const FloatArray &x, FloatArray &answer) const;
    /** Multiplies receiver by scalar value.
     * @param x value to multiply receiver
     */
    virtual void times(double x);
    /// Builds internal structure of receiver
    int buildInternalStructure(EngngModel *, int, EquationID);
    /**
     * Allocates and builds internal structure according to given
     * array holding adresses of diagonal members values (addr).
     */
    int setInternalStructure(IntArray *a);

    // int assemble (FloatMatrix*, IntArray*) ;
    /** Assembles receiver from local element contributions.
     * @param loc location array. The values corresponding to zero loc array value are not assembled.
     * @param mat contribution to be assembled using loc array.
     */
    int assemble(const IntArray &loc, const FloatMatrix &mat);
    /** Assembles receiver from local element contributions.
     * @param rloc row location array. The values corresponding to zero loc array value are not assembled.
     * @param cloc column location array. The values corresponding to zero loc array value are not assembled.
     * @param mat contribution to be assembled using loc array.
     */
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);

    /// Determines, whether receiver can be factorized.
    int canBeFactorized() const { return 1; }
    /// Performs factorization of receiver
    SparseMtrx *factorized();
    /**
     * Computes the solution of linear system \f$A x = \f$. A is receiver.
     * solution vector x overwrites the right hand side vector y.
     * Receiver must be in factorized form.
     * @param y right hand side on input, solution on output.
     * @return pointer to y array
     * @see factorized method
     */
    FloatArray *backSubstitutionWith(FloatArray &) const;
    /// Zeroes the receiver.
    SparseMtrx *zero();
    /**
     * Splits the receiver to LDLT form,
     * and computes the rigid body motions.
     * @param r matrix containing the rigid body motions base vectors.
     * @param nse number of rigid body motions
     * @param se array containing indexes of singular eqs.
     * @param limit - determines linear dependence or independence
     * @param   tc - typ vypoctu
     * = 1 - rozklad matice na tvar LDL
     * = 2 - konstruovani baze prostoru Ker A
     * = 3 - rozklad matice na tvar LDL a konstruovani baze prostoru Ker A
     */
    void rbmodes(FloatMatrix &r, int &nse, IntArray &se,
                 double limit, int tc);
    /**
     * Solves the singular system of equations, the receiver should be factorized
     * using rbmodes service.
     * @param x solution vector
     * @param y right hand side
     * @param nse number of rigid body motions
     * @param limit determines the liner dependency or independency
     * @param se indexes of singular equations
     */
    void ldl_feti_sky(FloatArray &x, FloatArray &y,
                      int nse, double limit, IntArray &se);
    /// Returns coefficient at position (i,j).
    double &at(int, int);
    /// Returns coefficient at position (i,j).
    double at(int i, int j) const;

    void toFloatMatrix(FloatMatrix &answer) const;
    /// Prints receiver to stdout.
    void printYourself() const;
    int  giveAllocatedSize() { return nwk; }

    SparseMtrxType  giveType() const { return SMT_Skyline; }
    int isAntisymmetric() const { return 0; }

#ifdef IML_COMPAT
    // /***********************************/
    // /*  Matrix/Vector multiply         */
    // /***********************************/

    virtual FloatArray trans_mult(const FloatArray &x) const
    { FloatArray answer;
      this->times(x, answer);
      return answer; }

#endif



protected:
    Skyline(int, int, double *, IntArray *);
};

#endif // skyline_h

