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

#ifndef beam3d_h
#define beam3d_h

#include "../sm/Elements/structuralelement.h"
#include "../sm/CrossSections/fiberedcs.h"

///@name Input fields for Beam3d
//@{
#define _IFT_Beam3d_Name "beam3d"
#define _IFT_Beam3d_dofstocondense "dofstocondense"
#define _IFT_Beam3d_refnode "refnode"
//@}

namespace oofem {
/**
 * This class implements a 2-dimensional beam element
 * with cubic lateral displacement interpolation (rotations are quadratic)
 * and longitudial displacements are linear.
 * This is an exact displacement approximation for beam with no
 * nonnodal loading.
 *
 * This class is not derived from liBeam3d or truss element, because it does not support
 * any material nonlinearities (if should, stiffness must be integrated)
 */
class Beam3d : public StructuralElement, public FiberedCrossSectionInterface
{
protected:
    double kappay, kappaz, length;
    int referenceNode;
    IntArray *dofsToCondense;

public:
    Beam3d(int n, Domain * d);
    virtual ~Beam3d();

    virtual void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL);
    virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual int giveLocalCoordinateSystem(FloatMatrix &answer);
    virtual void computeLocalForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *, int useUpdatedGpRecord = 0);
    void giveEndForcesVector(FloatArray &answer, TimeStep *tStep);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);


    virtual int testElementExtension(ElementExtension ext) {
        return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 );
    }
    //int hasLayeredSupport () {return 1;}

    virtual int computeNumberOfDofs() { return 12; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    //
    // fibered cross section support functions
    //
    virtual void FiberedCrossSectionInterface_computeStrainVectorInFiber(FloatArray &answer, const FloatArray &masterGpStrain,
                                                                         GaussPoint *slaveGp, TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType it);

    // definition & identification
    virtual const char *giveClassName() const { return "Beam3d"; }
    virtual const char *giveInputRecordName() const { return _IFT_Beam3d_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    ///@todo Introduce interpolator and remove these two:
    virtual integrationDomain giveIntegrationDomain() const { return _Line; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
#endif

protected:
    virtual void computeEdgeLoadVectorAt(FloatArray &answer, Load *, int, TimeStep *, ValueModeType mode);
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer);
    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);

    double giveKappayCoeff(TimeStep *tStep);
    double giveKappazCoeff(TimeStep *tStep);
    void computeKappaCoeffs(TimeStep *tStep);
    virtual double computeLength();
    virtual void computeClampedStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode, TimeStep *tStep);
    virtual void computeLocalStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseMode rMode, TimeStep *tStep);
    virtual void computeGaussPoints();
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    virtual MaterialMode giveMaterialMode() { return _3dBeam; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // beam3d_h
