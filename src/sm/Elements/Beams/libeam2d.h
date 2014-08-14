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

#ifndef libeam2d_h
#define libeam2d_h

#include "../sm/Elements/structuralelement.h"
#include "../sm/CrossSections/layeredcrosssection.h"

#define _IFT_LIBeam2d_Name "libeam2d"

namespace oofem {
class FEI2dLineLin;

/**
 * A 2-dimensional Linear Isoparametric
 * Mindlin theory beam element, with reduced integration.
 */
class LIBeam2d : public StructuralElement, public LayeredCrossSectionInterface
{
protected:
    /// Interpolation
    static FEI2dLineLin interpolation;

public:
    double pitch, length;

    LIBeam2d(int n, Domain * aDomain);
    virtual ~LIBeam2d() { }

    virtual FEInterpolation *giveInterpolation() const;

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) { computeLumpedMassMatrix(answer, tStep); }
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);

    // layered cross section support functions
    virtual void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                            GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    virtual Interface *giveInterface(InterfaceType it);

    virtual int computeNumberOfDofs() { return 6; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;
    virtual double computeVolumeAround(GaussPoint *gp);


    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_LIBeam2d_Name; }
    virtual const char *giveClassName() const { return "LIBeam2d"; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialMode giveMaterialMode() { return _2dBeam; }

protected:
    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *);
    virtual void giveEdgeDofMapping(IntArray &answer, int) const;
    virtual double computeEdgeVolumeAround(GaussPoint *, int);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &, int, GaussPoint *);
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer);
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);
    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &);
    virtual void computeGaussPoints();
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    virtual double computeLength();
    double givePitch();
};
} // end namespace oofem
#endif // libeam2d_h
