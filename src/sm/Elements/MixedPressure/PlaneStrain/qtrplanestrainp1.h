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

#ifndef qtrplanestrainp1_h
#define qtrplanestrainp1_h

#include "../sm/Elements/PlaneStrain/qtrplanestrain.h"
#include "../sm/Elements/MixedPressure/basemixedpressureelement.h"

#define _IFT_QTrPlaneStrainP1_Name "qtrplanestrainp1"


namespace oofem {
class FEI2dTrLin;

class QTrPlaneStrainP1 : public QTrPlaneStrain, public BaseMixedPressureElement
{
protected:
    static FEI2dTrLin interpolation_lin;

public:
    QTrPlaneStrainP1(int n, Domain *d);
    virtual ~QTrPlaneStrainP1() { }

protected:
    void computePressureNMatrixAt(GaussPoint *, FloatArray &) override;
    void computeVolumetricBmatrixAt(GaussPoint *gp, FloatArray &Bvol, NLStructuralElement *element) override;
    NLStructuralElement *giveElement() override { return this; }
public:
    const char *giveClassName() const override { return "QTrPlaneStrainP1"; }
    const char *giveInputRecordName() const override { return _IFT_QTrPlaneStrainP1_Name; }


    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    void giveDofManDofIDMask_u(IntArray &answer) override;
    void giveDofManDofIDMask_p(IntArray &answer) override;


    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep) override { BaseMixedPressureElement :: computeStiffnessMatrix(answer, mode, tStep); }
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override { BaseMixedPressureElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }


    int giveNumberOfPressureDofs() override { return 3; }
    int giveNumberOfDisplacementDofs() override { return 12; }
    int giveNumberOfDofs() override { return 15; }
    void postInitialize() override;
};
} // end namespace oofem
#endif // qtrplanestrainp1_h
