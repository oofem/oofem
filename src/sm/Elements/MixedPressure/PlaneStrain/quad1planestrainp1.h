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
 *               Copyright (C) 1993 - 2015   Borek Patzak
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

#ifndef quad1planestrainp1_h
#define quad1planestrainp1_h

#include "../sm/Elements/PlaneStrain/quad1planestrain.h"
#include "../sm/Elements/MixedPressure/basemixedpressureelement.h"


#define _IFT_Quad1PlaneStrainP1_Name "quad1planestrainp1"

namespace oofem {
class Quad1PlaneStrainP1 : public Quad1PlaneStrain, public BaseMixedPressureElement
{
protected:

public:
    Quad1PlaneStrainP1(int n, Domain *d);
    virtual ~Quad1PlaneStrainP1() { }

protected:
    void computePressureNMatrixAt(GaussPoint *gp, FloatArray &Np) override;
    void computeVolumetricBmatrixAt(GaussPoint *gp, FloatArray &answer, NLStructuralElement *element) override;
    NLStructuralElement *giveElement() override { return this; }

public:
    const char *giveInputRecordName() const override { return _IFT_Quad1PlaneStrainP1_Name; }
    const char *giveClassName() const override { return "Quad1PlaneStrainP0"; }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    void giveDofManDofIDMask_u(IntArray &answer) override;
    void giveDofManDofIDMask_p(IntArray &answer) override;

    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep) override { BaseMixedPressureElement :: computeStiffnessMatrix(answer, mode, tStep); }
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override { BaseMixedPressureElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }


    int giveNumberOfPressureDofs() override { return 4; }
    int giveNumberOfDisplacementDofs() override { return 8; }
    int giveNumberOfDofs() override { return 12; }
    void postInitialize() override;
};
} // end namespace oofem
#endif // quad1planestrainp1_h
