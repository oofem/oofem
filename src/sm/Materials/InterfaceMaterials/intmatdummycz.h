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

#ifndef INTMATDUMMYCZ_H_
#define INTMATDUMMYCZ_H_

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for IntMatDummyCZ
//@{
#define _IFT_IntMatDummyCZ_Name "intmatdummycz"
//@}

namespace oofem {

/**
 * Dummy cohesive zone model. The purpose of the model is
 * to store, and thus allow export of, the displacement jump.
 *
 * @date Jan 24, 2014
 * @author Erik Svenning
 */
class IntMatDummyCZ : public StructuralInterfaceMaterial
{
public:
    IntMatDummyCZ(int n, Domain *d);

    const char *giveClassName() const override { return "IntMatBilinearCZ"; }
    const char *giveInputRecordName() const override { return _IFT_IntMatDummyCZ_Name; }


    FloatArrayF<3> giveFirstPKTraction_3d(const FloatArrayF<3> &jump, const FloatMatrixF<3,3> &F, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<3,3> give3dStiffnessMatrix_dTdj(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    bool hasAnalyticalTangentStiffness() const override { return true; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<StructuralInterfaceMaterialStatus>(gp); }
    void printYourself() override;
};

} /* namespace oofem */
#endif /* INTMATDUMMYCZ_H_ */
