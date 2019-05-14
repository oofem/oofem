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
    virtual void computePressureNMatrixAt(GaussPoint *, FloatArray &);
    virtual void computeVolumetricBmatrixAt(GaussPoint *gp, FloatArray &Bvol, NLStructuralElement *element);
    virtual NLStructuralElement *giveElement() { return this; }
public:
    virtual const char *giveClassName() const { return "QTrPlaneStrainP1"; }
    virtual const char *giveInputRecordName() const { return _IFT_QTrPlaneStrainP1_Name; }


    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual void giveDofManDofIDMask_u(IntArray &answer);
    virtual void giveDofManDofIDMask_p(IntArray &answer);


    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep) { BaseMixedPressureElement :: computeStiffnessMatrix(answer, mode, tStep); }
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) { BaseMixedPressureElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }


    virtual int giveNumberOfPressureDofs() { return 3; }
    virtual int giveNumberOfDisplacementDofs() { return 12; }
    virtual int giveNumberOfDofs() { return 15; }
    virtual void postInitialize();
};
} // end namespace oofem
#endif // qtrplanestrainp1_h
