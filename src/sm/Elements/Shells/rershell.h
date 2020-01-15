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

#ifndef rershell_h
#define rershell_h

#include "sm/Elements/Shells/cct3d.h"
#include "sm/CrossSections/layeredcrosssection.h"

#define _IFT_RerShell_Name "rershell"

namespace oofem {
#ifndef __CHARTENSOR // termitovo
 #define __CHARTENSOR
enum CharTensor {
    LocalStrainTensor,
    GlobalStrainTensor,
    LocalCurvatureTensor,
    GlobalCurvatureTensor,

    LocalForceTensor,
    GlobalForceTensor,
    LocalMomentTensor,
    GlobalMomentTensor
};
#endif

/**
 * This class implements an triangular three-node  shell (CCT+linear plan stress)
 * curved finite element. Each node has 5 degrees of freedom.
 */
class RerShell : public CCTPlate3d
{
protected:
    double Rx, Ry, Rxy;

public:
    RerShell(int n, Domain * d);
    virtual ~RerShell() { }

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { computeLumpedMassMatrix(answer, tStep); }

    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    void giveLocalCoordinates(FloatArray &answer, const FloatArray &global);

    bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords) override;
    //
    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);
    void printOutputAt(FILE *file, TimeStep *tStep) override;

    const FloatMatrix *computeGtoLRotationMatrix() override { return CCTPlate3d::computeGtoLRotationMatrix(); }
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;

    // layered cross section support functions
    void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                    GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep) override;

    Interface *giveInterface(InterfaceType it) override;

    int computeNumberOfDofs() override { return 18; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep) override;

    // io routines
#ifdef __OOFEG
    //void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    //void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_RerShell_Name; }
    const char *giveClassName() const override { return "RerShell"; }

    void initializeFrom(InputRecord &ir) override;

    integrationDomain giveIntegrationDomain() const override { return _Triangle; }
    MaterialMode giveMaterialMode() override { return _3dShell; }

protected:
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode) override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeGaussPoints() override;
    double giveArea();
};
} // end namespace oofem
#endif // rershell_h
