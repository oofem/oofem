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

#ifndef libeam2d_h
#define libeam2d_h

#include "structuralelement.h"
#include "layeredcrosssection.h"
#include "fei2dlinelin.h"

namespace oofem {
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

    LIBeam2d(int n, Domain *aDomain);
    virtual ~LIBeam2d() { }

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) { computeLumpedMassMatrix(answer, tStep); }
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);

    // layered cross section support functions
    virtual void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                    GaussPoint *slaveGp, TimeStep *tStep);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    virtual Interface *giveInterface(InterfaceType it);

    virtual int computeNumberOfDofs(EquationID ut) { return 6; }
    virtual void giveDofManDofIDMask(int inode, EquationID eid, IntArray &) const;
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual FEInterpolation *giveInterpolation() { return &interpolation; }

    // definition & identification
    virtual const char *giveClassName() const { return "LIBeam2d"; }
    virtual classType giveClassID() const { return LIBeam2dClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual integrationDomain giveIntegrationDomain() { return _Line; }
    virtual MaterialMode giveMaterialMode() { return _2dBeam; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

protected:
    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *);
    virtual void giveEdgeDofMapping(IntArray &answer, int) const;
    virtual double computeEdgeVolumeAround(GaussPoint *, int);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
    { computeGlobalCoordinates( answer, * ( gp->giveCoordinates() ) ); }
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &, int, GaussPoint *);
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer);
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);
    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(GaussPoint *, FloatMatrix &);
    virtual void computeGaussPoints();
    double giveLength();
    double givePitch();
};
} // end namespace oofem
#endif // libeam2d_h
