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

#include "termlibrary5.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "engngm.h"

namespace oofem {

REGISTER_Term(OperatorTerm)

namespace {
Material *giveMaterial(MPElement &cell, GaussPoint *gp)
{
    return cell.giveCrossSection()->giveMaterial(gp);
}
} // end anonymous namespace


OperatorTerm::OperatorTerm(const Variable *testField, StateOperator testOp, const Variable *unknownField, StateOperator op,
                           MatResponseMode lhsMode, MatResponseMode rhsMode, double sign, ValueModeType vm,
                           bool hasLhs, bool hasRhs) :
    GenericCellTerm(testField, unknownField, _Unknown), testOp(testOp), op(op), lhsMode(lhsMode), rhsMode(rhsMode),
    lhsFromMaterial(true), rhsFromMaterial(true), sign(sign), vm(vm), hasLhs(hasLhs), hasRhs(hasRhs)
{ }


OperatorTerm::OperatorTerm(const Variable *testField, StateOperator testOp, const Variable *unknownField, StateOperator op,
                           double coefficient, double sign, ValueModeType vm) :
    GenericCellTerm(testField, unknownField, _Unknown), testOp(testOp), op(op), coefficient(coefficient), sign(sign), vm(vm)
{ }


void
OperatorTerm::giveTestCoefficientProduct(FloatMatrix &answer, const FloatMatrix &a, bool fromMaterial, MatResponseMode mode,
                                         MPElement &cell, GaussPoint *gp, TimeStep *tStep) const
{
    if ( fromMaterial ) {
        FloatMatrix d;
        giveMaterial(cell, gp)->giveCharacteristicMatrix(d, mode, gp, tStep);
        answer.beTProductOf(a, d);
    } else {
        answer.beTranspositionOf(a);
        answer.times(this->coefficient);
    }
}


void
OperatorTerm::evaluate_lin(FloatMatrix &answer, MPElement &cell, GaussPoint *gp, TimeStep *tStep) const
{
    if ( !this->hasLhs ) {
        answer.clear();
        return;
    }
    FloatMatrix a, b, atd;
    cell.computeStateOperatorAt(a, this->testOp, this->testField, gp);
    cell.computeStateOperatorAt(b, this->op, this->field, gp);
    this->giveTestCoefficientProduct(atd, a, this->lhsFromMaterial, this->lhsMode, cell, gp, tStep);
    answer.beProductOf(atd, b);
    answer.times(this->sign);
}


void
OperatorTerm::evaluate(FloatArray &answer, MPElement &cell, GaussPoint *gp, TimeStep *tStep) const
{
    if ( !this->hasRhs ) {
        answer.clear();
        return;
    }
    FloatMatrix a, b, atd, k;
    FloatArray r;
    cell.computeStateOperatorAt(a, this->testOp, this->testField, gp);
    cell.computeStateOperatorAt(b, this->op, this->field, gp);
    this->giveTestCoefficientProduct(atd, a, this->rhsFromMaterial, this->rhsMode, cell, gp, tStep);
    k.beProductOf(atd, b);
    cell.getUnknownVector(r, this->field, this->vm, tStep);
    answer.beProductOf(k, r);
    answer.times(this->sign);
}


void
OperatorTerm::initializeFrom(const std::shared_ptr< InputRecord > &ir, EngngModel *problem)
{
    GenericCellTerm::initializeFrom(ir, problem);
    int value;
    IR_GIVE_FIELD(ir, value, _IFT_OperatorTerm_testop);
    this->testOp = static_cast< StateOperator >(value);
    IR_GIVE_FIELD(ir, value, _IFT_OperatorTerm_op);
    this->op = static_cast< StateOperator >(value);

    this->lhsFromMaterial = ir->hasField(_IFT_OperatorTerm_lmode);
    if ( this->lhsFromMaterial ) {
        IR_GIVE_FIELD(ir, value, _IFT_OperatorTerm_lmode);
        this->lhsMode = static_cast< MatResponseMode >(value);
    }
    this->rhsFromMaterial = ir->hasField(_IFT_OperatorTerm_rmode) || this->lhsFromMaterial;
    this->rhsMode = this->lhsMode;
    if ( ir->hasField(_IFT_OperatorTerm_rmode) ) {
        IR_GIVE_FIELD(ir, value, _IFT_OperatorTerm_rmode);
        this->rhsMode = static_cast< MatResponseMode >(value);
    }
    IR_GIVE_OPTIONAL_FIELD(ir, this->coefficient, _IFT_OperatorTerm_coefficient);
    IR_GIVE_OPTIONAL_FIELD(ir, this->sign, _IFT_OperatorTerm_sign);
    value = this->vm;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_OperatorTerm_vm);
    this->vm = static_cast< ValueModeType >(value);
    this->hasLhs = !ir->hasField(_IFT_OperatorTerm_nolhs);
    this->hasRhs = !ir->hasField(_IFT_OperatorTerm_norhs);
}

} // end namespace oofem
