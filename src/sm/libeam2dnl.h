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

#ifndef libeam2dnl_h
#define libeam2dnl_h

#include "nlstructuralelement.h"
#include "layeredcrosssection.h"

namespace oofem {
/**
 * This class implements a 2-dimensional Linear Isoparametric
 * Mindlin theory beam element, with reduced integration.
 * Geometric nonlinearities are taken into account.
 */
class LIBeam2dNL : public NLStructuralElement, public LayeredCrossSectionInterface
{
protected:
    double pitch, length;

public:
    LIBeam2dNL(int n, Domain *d);
    ~LIBeam2dNL() { }

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
    int computeGtoLRotationMatrix(FloatMatrix &answer);
    void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    // layered cross section support functions
    void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                            GaussPoint *slaveGp, TimeStep *tStep);

    Interface *giveInterface(InterfaceType it);

    virtual int computeNumberOfDofs(EquationID ut) { return 6; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    double computeVolumeAround(GaussPoint *gp);

    // definition & identification
    const char *giveClassName() const { return "LIBeam2dNL"; }
    classType giveClassID() const { return LIBeam2dClass; }
    IRResultType initializeFrom(InputRecord *ir);

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
#endif

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

    //void computeTemperatureStrainVectorAt(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode);
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    // nonlinear part of geometrical eqs. for i-th component of strain vector.
    void computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, int);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeGaussPoints();
    double giveLength();
    double givePitch();
};
} // end namespace oofem
#endif // libeam2dnl_h
