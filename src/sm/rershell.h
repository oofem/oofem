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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef rershell_h
#define rershell_h

#include "cct.h"
#include "layeredcrosssection.h"

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
class RerShell : public CCTPlate
{
protected:
    double Rx, Ry, Rxy;
    FloatMatrix *GtoLRotationMatrix;

    // Transformation Matrix form GtoL(3,3) is stored
    // at the element level for computation
    // efficiency

public:
    RerShell(int n, Domain *d);
    ~RerShell() { delete GtoLRotationMatrix; }

    // FloatMatrix* ComputeConstitutiveMatrixAt (GaussPoint *) ;
    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
    //void printOutputAt (TimeStep*) ;
    FloatMatrix *computeGtoLRotationMatrix();
    int giveLocalCoordinateSystem(FloatMatrix &answer);
    void giveLocalCoordinates(FloatArray &answer, const FloatArray &global);

    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    //
    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);
    void printOutputAt(FILE *file, TimeStep *tStep);

    // layered cross section support functions
    void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                                  GaussPoint *slaveGp, TimeStep *tStep);

    Interface *giveInterface(InterfaceType it);

    virtual int computeNumberOfDofs(EquationID ut) { return 18; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

    int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    Element *ZZNodalRecoveryMI_giveElement() { return this; }
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint,
                                                             InternalStateType type);

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }


    // io routines
#ifdef __OOFEG
    //void drawRawGeometry(oofegGraphicContext &);
    //void drawDeformedGeometry(oofegGraphicContext &);
    //virtual void drawScalar(oofegGraphicContext &context);
    //void drawInternalState(oofegGraphicContext &);
#endif

    // definition & identification
    const char *giveClassName() const { return "RerShell"; }
    classType giveClassID() const { return RerShellClass; }

    IRResultType initializeFrom(InputRecord *ir);

    integrationDomain giveIntegrationDomain() { return _Triangle; }
    MaterialMode giveMaterialMode() { return _3dShell; }

protected:
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);
    //void computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    int computeGtoLRotationMatrix(FloatMatrix &answer);
    void computeGaussPoints();
    double giveArea();
};
} // end namespace oofem
#endif // rershell_h
