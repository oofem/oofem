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

#ifndef libeam2d_h
#define libeam2d_h

#include "structuralelement.h"
#include "layeredcrosssection.h"
#include "fei2dlinelin.h"

namespace oofem {
/**
 * A 2-dimensional Linear Isoparametric.
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
    ~LIBeam2d() { }

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) { computeLumpedMassMatrix(answer, tStep); }
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    int computeGtoLRotationMatrix(FloatMatrix &);

    // layered cross section support functions
    void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                    GaussPoint *slaveGp, TimeStep *tStep);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    Interface *giveInterface(InterfaceType it);

    virtual int computeNumberOfDofs(EquationID ut) { return 6; }
    virtual void giveDofManDofIDMask(int inode, EquationID eid, IntArray &) const;
    double computeVolumeAround(GaussPoint *gp);

    virtual FEInterpolation *giveInterpolation() { return &interpolation; }

    // definition & identification
    const char *giveClassName() const { return "LIBeam2d"; }
    classType giveClassID() const { return LIBeam2dClass; }
    IRResultType initializeFrom(InputRecord *ir);

    integrationDomain giveIntegrationDomain() { return _Line; }
    MaterialMode giveMaterialMode() { return _2dBeam; }

protected:
    // edge load support
    void computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *);
    void giveEdgeDofMapping(IntArray &answer, int) const;
    double computeEdgeVolumeAround(GaussPoint *, int);
    void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
    { computeGlobalCoordinates( answer, * ( gp->giveCoordinates() ) ); }
    int computeLoadLEToLRotationMatrix(FloatMatrix &, int, GaussPoint *);
    int computeLoadGToLRotationMtrx(FloatMatrix &answer);
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *, FloatMatrix &);
    void computeGaussPoints();
    double giveLength();
    double givePitch();
};
} // end namespace oofem
#endif // libeam2d_h
