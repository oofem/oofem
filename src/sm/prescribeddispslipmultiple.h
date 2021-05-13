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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#ifndef PRESCRIBEDDISPSLIPBCMULTIPLE_H_
#define PRESCRIBEDDISPSLIPBCMULTIPLE_H_

#include "prescribeddispsliphomogenization.h"
#include "activebc.h"

#include <memory>

#define _IFT_PrescribedDispSlipMultiple_Name   "prescribeddispslipmultiple"
#define _IFT_PrescribedDispSlipMultiple_BCs    "bcs"

namespace oofem {
class Node;
class Element;
/**
 * Allows for imposing multiple boundary conditions of type PrescribedDispSlip on
 * the RVE
 *
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT PrescribedDispSlipMultiple : public ActiveBoundaryCondition, public PrescribedDispSlipHomogenization
{
public:
    PrescribedDispSlipMultiple(int n, Domain *d);
    virtual ~PrescribedDispSlipMultiple();

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    bcType giveType() const override { return UnknownBT; }
    DofManager *giveInternalDofManager(int i) override;

    void setSlipField(const FloatArray &t) override;
    void setDispGradient(const FloatArray &t) override;
    void setSlipGradient(const FloatArray &t) override;
    void setCenterCoordinate(FloatArray &x) override;
    void scale(double s) override;

    const char *giveClassName() const override { return "PrescribedDispSlipMultiple"; }
    const char *giveInputRecordName() const override { return _IFT_PrescribedDispSlipMultiple_Name; }

    void computeStress(FloatArray &stress, TimeStep *tStep) override;
    void computeTransferStress(FloatArray &bStress, TimeStep *tStep) override;
    void computeReinfStress(FloatArray &rStress, TimeStep *tStep) override;

    void computeTangent(FloatMatrix &tangent, TimeStep *tStep) override;

protected:
    IntArray bcs;
};
} /* namespace oofem */

#endif /* PRESCRIBEDDISPSLIPBCMULTIPLE_H_ */
