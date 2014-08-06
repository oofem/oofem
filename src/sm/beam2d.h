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

#ifndef beam2d_h
#define beam2d_h

#include "structuralelement.h"
#include "layeredcrosssection.h"
#include "fei2dlinelin.h"
#include "fei2dlinehermite.h"

///@name Input fields for Beam2d
//@{
#define _IFT_Beam2d_Name "beam2d"
#define _IFT_Beam2d_dofstocondense "dofstocondense"
//@}

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
    Beam2d(int n, Domain * aDomain);
    virtual ~Beam2d();

    virtual void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL);
    virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual int giveLocalCoordinateSystem(FloatMatrix &answer);
    virtual void computeLocalForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual void giveEndForcesVector(FloatArray &answer, TimeStep *tStep);

    virtual int testElementExtension(ElementExtension ext) { return ( ext == Element_EdgeLoadSupport ); }

    virtual Interface *giveInterface(InterfaceType);

    virtual FEInterpolation *giveInterpolation() const { return & interp_geom; }
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const { return NULL; }

    virtual int computeNumberOfDofs() { return 6; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual void  printOutputAt(FILE *file, TimeStep *tStep);

    virtual const char *giveClassName() const { return "Beam2d"; }
    virtual const char *giveInputRecordName() const { return _IFT_Beam2d_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
#endif

    virtual void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                            GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

protected:
    virtual void computeEdgeLoadVectorAt(FloatArray &answer, Load *, int, TimeStep *, ValueModeType mode);
    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);

    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);

    double giveKappaCoeff(TimeStep *tStep);
    virtual double computeLength();
    double givePitch();
    virtual void computeClampedStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode, TimeStep *tStep);
    virtual void computeLocalStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseMode rMode, TimeStep *tStep);
    virtual void computeGaussPoints();
    virtual MaterialMode giveMaterialMode() { return _2dBeam; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // beam2d_h
