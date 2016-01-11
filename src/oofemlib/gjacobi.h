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

#ifndef gjacobi_h
#define gjacobi_h

#include "nummet.h"
#include "floatarray.h"
#include "nmstatus.h"

namespace oofem {
class Domain;
class EngngModel;

/**
 * This class implements the Generalized Jacobi Eigenvalue Problem Solver.
 * Perform solution of eigenvalue problem in the form
 * @f[
 *  K\cdot x = \omega^2 M\cdot x
 * @f]
 */
class OOFEM_EXPORT GJacobi : public NumericalMethod
{
private:
    int n, nsmax;
    double rtol;
    int solved;

public:
    GJacobi(Domain * d, EngngModel * m);
    virtual ~GJacobi();

    /**
     * Solves the given sparse generalized eigenvalue system of equations @f$ K\cdot x = w^2 M\cdot x @f$.
     * @param K Coefficient matrix.
     * @param M Coefficient matrix.
     * @param x Eigenvector(s).
     * @param w Eigenvalue(s).
     * @return Status.
     */
    virtual NM_Status solve(FloatMatrix &K, FloatMatrix &M, FloatArray &w, FloatMatrix &x);

    virtual const char *giveClassName() const { return "GeneralizedJacobiSolver"; }
    std :: string errorInfo(const char *func) { return std :: string(this->giveClassName()) + func; }
};
} // end namespace oofem
#endif // gjacobi_h
