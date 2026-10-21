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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
#ifndef termlibrary5_h
#define termlibrary5_h

/**
 * Hardwired terms of MPM weak forms, built from the state operators of the cells (MPElement) and
 * constitutive coefficients queried from the material by MatResponseMode.
 *
 * Each term is the hardwired counterpart of an expression of the symbolic formulation; the
 * documentation of each class gives the expression.
 */

#include "mpm.h"
#include "matresponsemode.h"
#include "valuemodetype.h"
#include "stateoperator.h"

#include <vector>

#define _IFT_OperatorTerm_Name "OperatorTerm"
#define _IFT_OperatorTerm_testop "testop"
#define _IFT_OperatorTerm_op "op"
#define _IFT_OperatorTerm_lmode "lmode"
#define _IFT_OperatorTerm_rmode "rmode"
#define _IFT_OperatorTerm_coefficient "coefficient"
#define _IFT_OperatorTerm_sign "sign"
#define _IFT_OperatorTerm_vm "vm"
#define _IFT_OperatorTerm_nolhs "nolhs"
#define _IFT_OperatorTerm_norhs "norhs"

namespace oofem {

/**
 * Bilinear term sign * A(w)^T D B(u), where A and B are operators (value, gradient, symmetric
 * gradient, divergence) applied to the test field w and the unknown u, and D is a constitutive
 * coefficient queried from the material (or a constant). The residual is
 * sign * A^T D B r, r being the nodal values of the unknown -- the total ones, or their rates for
 * the terms of the rate part of the weak form.
 *
 * The residual may use a different coefficient than the tangent (a stabilized tangent), and either
 * part may be absent.
 *
 * Symbolic equivalent: "sign*A(w,gp).T * MDer(gp,ts,lmode) * B(u,gp)" and
 * "sign*A(w,gp).T * MDer(gp,ts,rmode) * B(u,gp) * ru(u,cell,ts)" (rv for rates).
 *
 * Input record: OperatorTerm n variable "u" testvariable "w" mmode m testop o op o
 * [lmode m] [rmode m] [coefficient c] [sign s] [vm v] [nolhs] [norhs] [nip n].
 * The operators are StateOperator values; without lmode (rmode), the coefficient is the scalar c
 * (default 1).
 */
class OperatorTerm : public GenericCellTerm
{
protected:
    StateOperator testOp = SO_Value, op = SO_Value;
    MatResponseMode lhsMode = MatResponseMode::MRM_ScalarOne, rhsMode = MatResponseMode::MRM_ScalarOne;
    bool lhsFromMaterial = false, rhsFromMaterial = false;
    double coefficient = 1.0;
    double sign = 1.0;
    ValueModeType vm = VM_TotalIntrinsic;
    bool hasLhs = true, hasRhs = true;

public:
    OperatorTerm() : GenericCellTerm() { }
    /// Term with material coefficients.
    OperatorTerm(const Variable *testField, StateOperator testOp, const Variable *unknownField, StateOperator op,
                 MatResponseMode lhsMode, MatResponseMode rhsMode, double sign = 1.0, ValueModeType vm = VM_TotalIntrinsic,
                 bool hasLhs = true, bool hasRhs = true);
    /// Term with a constant coefficient.
    OperatorTerm(const Variable *testField, StateOperator testOp, const Variable *unknownField, StateOperator op,
                 double coefficient, double sign = 1.0, ValueModeType vm = VM_TotalIntrinsic);

    void evaluate_lin(FloatMatrix &answer, MPElement &cell, GaussPoint *gp, TimeStep *tStep) const override;
    void evaluate(FloatArray &answer, MPElement &cell, GaussPoint *gp, TimeStep *tStep) const override;
    void getDimensions(Element &cell) const override { }
    void initializeFrom(const std::shared_ptr< InputRecord > &ir, EngngModel *problem) override;

protected:
    /// Computes A^T D, D being the coefficient: a material matrix, or the constant (a scalar factor).
    void giveTestCoefficientProduct(FloatMatrix &answer, const FloatMatrix &a, bool fromMaterial, MatResponseMode mode,
                                    MPElement &cell, GaussPoint *gp, TimeStep *tStep) const;
};

} // end namespace oofem
#endif // termlibrary5_h
