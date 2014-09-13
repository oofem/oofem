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

#ifndef libeam3d_h
#define libeam3d_h

#include "../sm/Elements/structuralelement.h"
#include "../sm/CrossSections/fiberedcs.h"

///@name Input fields for LIBeam3d
//@{
#define _IFT_LIBeam3d_Name "libeam3d"
#define _IFT_LIBeam3d_refnode "refnode"
//@}

namespace oofem {
/**
 * This class implements a 3-dimensional mindlin theory Linear Isoparametric
 * beam element, with reduced integration.
 */
class LIBeam3d : public StructuralElement, public FiberedCrossSectionInterface
{
private:
    double length;
    int referenceNode;

public:
    LIBeam3d(int n, Domain * d);
    virtual ~LIBeam3d() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);

    virtual int testElementExtension(ElementExtension ext);

    virtual int computeNumberOfDofs() { return 12; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual int giveLocalCoordinateSystem(FloatMatrix &answer);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    // Fibered cross section support functions
    virtual void FiberedCrossSectionInterface_computeStrainVectorInFiber(FloatArray &answer, const FloatArray &masterGpStrain,
                                                                         GaussPoint *slaveGp, TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType it);

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_LIBeam3d_Name; }
    virtual const char *giveClassName() const { return "LIBeam3d"; }

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
#endif

    virtual integrationDomain giveIntegrationDomain() const { return _Line; }
    virtual MaterialMode giveMaterialMode() { return _3dBeam; }

protected:
    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp);
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer);
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);

    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    virtual void computeGaussPoints();
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    virtual double computeLength();
};
} // end namespace oofem
#endif // libeam3d_h
