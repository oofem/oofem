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

#include "../sm/Elements/Shells/cct3d.h"
#include "../sm/CrossSections/layeredcrosssection.h"

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
    LocalMomentumTensor,
    GlobalMomentumTensor
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
    virtual ~RerShell() {
    }

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }

    virtual int giveLocalCoordinateSystem(FloatMatrix &answer);

    void giveLocalCoordinates(FloatArray &answer, const FloatArray &global);

    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    //
    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual const FloatMatrix *computeGtoLRotationMatrix() {return CCTPlate3d::computeGtoLRotationMatrix();}
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    


    // layered cross section support functions
    virtual void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                            GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType it);

    virtual int computeNumberOfDofs() { return 18; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep);


    // io routines
#ifdef __OOFEG
    //void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    //virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
#endif

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_RerShell_Name; }
    virtual const char *giveClassName() const { return "RerShell"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual integrationDomain giveIntegrationDomain() const { return _Triangle; }
    virtual MaterialMode giveMaterialMode() { return _3dShell; }

protected:
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void computeGaussPoints();
    virtual double giveArea();
};
} // end namespace oofem
#endif // rershell_h
