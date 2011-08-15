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

#ifndef beam2d_h
#define beam2d_h

#include "structuralelement.h"
#include "layeredcrosssection.h"
#include "fei2dlinelin.h"
//#include "fei2dlinehermite.h"

namespace oofem {
/**
 * This class implements a 2-dimensional beam element
 * with cubic lateral displacement and geometry, quadratic rotations,
 * and linear longitudinal displacements.
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
    //static FEI2dLineHermite interp_beam;

public:
    Beam2d(int n, Domain *aDomain);
    ~Beam2d();

    void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass);
    void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    int giveLocalCoordinateSystem(FloatMatrix &answer);
    void computeLocalForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    void giveEndForcesVector(FloatArray &answer, TimeStep *tStep);
    /**
     * Computes the global coordinates from given element's local coordinates.
     * @returns Nonzero if successful.
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }

    Interface *giveInterface(InterfaceType);

    virtual FEInterpolation *giveInterpolation() { return &interp_geom; }
    virtual FEInterpolation *giveInterpolation(DofIDItem id) { return NULL; }

    virtual int computeNumberOfDofs(EquationID ut) { return 6; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    double computeVolumeAround(GaussPoint *gp);
    void  printOutputAt(FILE *file, TimeStep *tStep);

    const char *giveClassName() const { return "Beam2d"; }
    classType giveClassID() const { return Beam2dClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
#endif

    void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                     GaussPoint *slaveGp, TimeStep *tStep);

protected:
    void computeEdgeLoadVectorAt(FloatArray &answer, Load *, int, TimeStep *, ValueModeType mode);
    void computePrescribedStrainLocalLoadVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    //void computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *, FloatMatrix &);
    int computeGtoLRotationMatrix(FloatMatrix &);
    //int computeGtoNRotationMatrix (FloatMatrix&);

    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);

    double giveKappaCoeff();
    double giveLength();
    double givePitch();
    virtual void computeClampedStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode, TimeStep *tStep);
    virtual void computeLocalStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseMode rMode, TimeStep *tStep);
    void computeGaussPoints();
    integrationDomain giveIntegrationDomain() { return _Line; }
    MaterialMode giveMaterialMode() { return _2dBeam; }
    virtual int  giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // beam2d_h
