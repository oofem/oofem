/* $Header: /home/cvs/bp/oofem/sm/src/beam3d.h,v 1.5.4.1 2004/04/05 15:19:46 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ********************
//   *** CLASS Beam3d ***
//   ********************

#ifndef beam3d_h
#define beam3d_h

#include "structuralelement.h"
#include "fiberedcs.h"

class Beam3d : public StructuralElement, public FiberedCrossSectionInterface
{
    /*
     * This class implements a 2-dimensional beam element
     * with cubic lateral displace,ent interpolation (rotations are quadratic)
     * and longitudial displacements are linear.
     * This is an exact displacement approximation for beam with no
     * nonnodal loading.
     *
     * This class is not derived from liBeam3d or truss element, because it does not support
     * any material nonlinearities (if shoul, stiffness must be integrated)
     */
protected:
    double kappay, kappaz, length;
    int referenceNode;
    IntArray *dofsToCondense;

public:
    Beam3d(int, Domain *);                     // constructor
    ~Beam3d();                                 // destructor

    void          computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass);
    void          computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);
    void          computeStiffnessMatrix(FloatMatrix &answer,
                                         MatResponseMode rMode, TimeStep *tStep);
    int           giveLocalCoordinateSystem(FloatMatrix &answer);
    void          computeLocalForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    void          giveInternalForcesVector(FloatArray &answer,
                                           TimeStep *, int useUpdatedGpRecord = 0);
    void          giveEndForcesVector(FloatArray &answer, TimeStep *tStep);
    /**
     * Computes the global coordinates from given element's local coordinates.
     * @returns nonzero if successful
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);


    virtual int testElementExtension(ElementExtension ext) {
        return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 );
    }
    //int hasLayeredSupport () {return 1;}

    virtual int            computeNumberOfDofs(EquationID ut) { return 12; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    double        computeVolumeAround(GaussPoint *);

    void  printOutputAt(FILE *, TimeStep *);

    //
    // fibered cross section support functions
    //
    void FiberedCrossSectionInterface_computeStrainVectorInFiber(FloatArray &answer, GaussPoint *masterGp,
                                                                 GaussPoint *slaveGp, TimeStep *tStep);


    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);
    //
    // definition & identification
    //
    const char *giveClassName() const { return "Beam3d"; }
    classType            giveClassID() const { return Beam3dClass; }
    IRResultType initializeFrom(InputRecord *ir);

#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
#endif

protected:
    void          computeEdgeLoadVectorAt(FloatArray &answer, Load *, int, TimeStep *, ValueModeType mode);
    int           computeLoadGToLRotationMtrx(FloatMatrix &answer);
    void          computePrescribedStrainLocalLoadVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    //void          computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);
    void          computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void          computeNmatrixAt(GaussPoint *, FloatMatrix &);
    int           computeGtoLRotationMatrix(FloatMatrix &); // giveRotationMatrix () ;
    //  int           computeGtoNRotationMatrix (FloatMatrix&);

    double giveKappayCoeff();
    double giveKappazCoeff();
    void   computeKappaCoeffs();
    double        giveLength();
    virtual void computeClampedStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode, TimeStep *tStep);
    virtual void computeLocalStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseMode rMode, TimeStep *tStep);
    void          computeGaussPoints();
    integrationDomain  giveIntegrationDomain() { return _Line; }
    virtual int  giveNumberOfIPForMassMtrxIntegration() { return 4; }
};

#endif // beam3d_h
