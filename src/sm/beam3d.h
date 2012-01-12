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

#ifndef beam3d_h
#define beam3d_h

#include "structuralelement.h"
#include "fiberedcs.h"

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
    Beam3d(int n, Domain *d);
    ~Beam3d();

    void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass);
    void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeStiffnessMatrix(FloatMatrix &answer,
                                         MatResponseMode rMode, TimeStep *tStep);
    int giveLocalCoordinateSystem(FloatMatrix &answer);
    void computeLocalForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    void giveInternalForcesVector(FloatArray &answer, TimeStep *, int useUpdatedGpRecord = 0);
    void giveEndForcesVector(FloatArray &answer, TimeStep *tStep);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);


    virtual int testElementExtension(ElementExtension ext) {
        return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 );
    }
    //int hasLayeredSupport () {return 1;}

    virtual int computeNumberOfDofs(EquationID ut) { return 12; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    double computeVolumeAround(GaussPoint *gp);

    void printOutputAt(FILE *file, TimeStep *tStep);

    //
    // fibered cross section support functions
    //
    void FiberedCrossSectionInterface_computeStrainVectorInFiber(FloatArray &answer, GaussPoint *masterGp,
                                                                 GaussPoint *slaveGp, TimeStep *tStep);

    Interface *giveInterface(InterfaceType it);

    // definition & identification
    const char *giveClassName() const { return "Beam3d"; }
    classType giveClassID() const { return Beam3dClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
#endif

protected:
    void computeEdgeLoadVectorAt(FloatArray &answer, Load *, int, TimeStep *, ValueModeType mode);
    int computeLoadGToLRotationMtrx(FloatMatrix &answer);
    void computePrescribedStrainLocalLoadVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    //void computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *, FloatMatrix &);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);

    double giveKappayCoeff();
    double giveKappazCoeff();
    void computeKappaCoeffs();
    double giveLength();
    virtual void computeClampedStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode, TimeStep *tStep);
    virtual void computeLocalStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseMode rMode, TimeStep *tStep);
    void computeGaussPoints();
    integrationDomain giveIntegrationDomain() { return _Line; }
    MaterialMode giveMaterialMode() { return _3dBeam; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // beam3d_h
