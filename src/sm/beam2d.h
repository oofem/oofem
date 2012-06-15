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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef beam2d_h
#define beam2d_h

#include "structuralelement.h"
#include "layeredcrosssection.h"
#include "fei2dlinelin.h"
#include "fei2dlinehermite.h"

namespace oofem {
/**
 * This class implements a 2-dimensional beam element
 * with cubic lateral displacement, quadratic rotations,
 * and linear longitudinal displacements and geometry.
 * This is an exact displacement approximation for beam with no
 * nonnodal loading.
 *
 * This class is not derived from linear beam or truss element, because it does not support
 * any material nonlinearities (if should, stiffness must be integrated)
 */
class Beam2d : public StructuralElement, public LayeredCrossSectionInterface
{
protected:
    double kappa, pitch, length;
    IntArray *dofsToCondense;

    static FEI2dLineLin interp_geom;
    static FEI2dLineHermite interp_beam;

public:
    Beam2d(int n, Domain *aDomain);
    virtual ~Beam2d();

    virtual void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass);
    virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual int giveLocalCoordinateSystem(FloatMatrix &answer);
    virtual void computeLocalForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual void giveEndForcesVector(FloatArray &answer, TimeStep *tStep);

    virtual int testElementExtension(ElementExtension ext) { return ( ext == Element_EdgeLoadSupport ); }

    virtual Interface *giveInterface(InterfaceType);

    virtual FEInterpolation *giveInterpolation() { return &interp_geom; }
    virtual FEInterpolation *giveInterpolation(DofIDItem id) { return NULL; }

    virtual int computeNumberOfDofs(EquationID ut) { return 6; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual void  printOutputAt(FILE *file, TimeStep *tStep);

    virtual const char *giveClassName() const { return "Beam2d"; }
    virtual classType giveClassID() const { return Beam2dClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
#endif

    virtual void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                     GaussPoint *slaveGp, TimeStep *tStep);

protected:
    virtual void computeEdgeLoadVectorAt(FloatArray &answer, Load *, int, TimeStep *, ValueModeType mode);
    virtual void computePrescribedStrainLocalLoadVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(GaussPoint *, FloatMatrix &);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);

    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);

    double giveKappaCoeff();
    double giveLength();
    double givePitch();
    virtual void computeClampedStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode, TimeStep *tStep);
    virtual void computeLocalStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseMode rMode, TimeStep *tStep);
    virtual void computeGaussPoints();
    virtual integrationDomain giveIntegrationDomain() { return _Line; }
    virtual MaterialMode giveMaterialMode() { return _2dBeam; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // beam2d_h
