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

#ifndef tr1bubblestokes_h
#define tr1bubblestokes_h

#include "fmelement.h"
#include "domain.h"
#include "zznodalrecoverymodel.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"
#include "elementinternaldofman.h"

#define _IFT_Tr1BubbleStokes_Name "tr1bubblestokes"

namespace oofem {

class FEI2dTrLin;

/**
 * Triangular element for Stokes flow using Bubble basis function.
 * Linear+Bubble interpolation of velocity, and linear interpolation of pressures.
 * The element is exported as a linear triangle (bubble dofs are not exported).
 * It can deal with nonlinear material models, but it is assumed that the fluid is without memory (which is usually the case).
 * 
 * @author Mikael Öhman
 */
class Tr1BubbleStokes : public FMElement,
    public ZZNodalRecoveryModelInterface,
    public SpatialLocalizerInterface,
    public EIPrimaryUnknownMapperInterface
{
protected:
    /// Number of gauss points. Same for pressure and velocity.
    int numberOfGaussPoints;
    /// Interpolation for pressure
    static FEI2dTrLin interp;
    /// Ordering of dofs in element. Used to assemble the element stiffness
    static IntArray momentum_ordering, conservation_ordering;
    /// Ordering of dofs on edges. Used to assemble edge loads
    static IntArray edge_ordering [ 3 ];
    /// Dummy variable
    static bool __initialized;
    /// Defines the ordering of the dofs in the local stiffness matrix.
    static bool initOrdering() {
        momentum_ordering.setValues(8,  1, 2, 4, 5, 7, 8, 10, 11);
        conservation_ordering.setValues(3,  3, 6, 9);
        edge_ordering [ 0 ].setValues(4,  1, 2, 4, 5);
        edge_ordering [ 1 ].setValues(4,  4, 5, 7, 8);
        edge_ordering [ 2 ].setValues(4,  7, 8, 1, 2);
        return true;
    }

    /// The extra dofs from the bubble
    ElementDofManager *bubble;
    // Coordinates associated with the bubble dofs.
    //FloatArray bubbleCoord; // Assumed fixed at 0 for now (i.e. only linear geometry)

public:
    Tr1BubbleStokes(int n, Domain *d);
    virtual ~Tr1BubbleStokes();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual void computeGaussPoints();

    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);

    void computeInternalForcesVector(FloatArray &answer, TimeStep *tStep);
    void computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep);

    void computeExternalForcesVector(FloatArray &answer, TimeStep *tStep);
    virtual void computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void computeBoundaryLoadVector(FloatArray &answer, Load *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep);

    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }
    virtual const char *giveClassName() const { return "Tr1BubbleStokes"; }
    virtual const char *giveInputRecordName() const { return _IFT_Tr1BubbleStokes_Name; }
    virtual classType giveClassID() const { return Tr1BubbleStokesElementClass; }
    virtual MaterialMode giveMaterialMode() { return _2dFlow; }

    virtual int computeNumberOfDofs(EquationID ut);

    virtual int giveNumberOfInternalDofManagers() const { return 1; }
    virtual DofManager *giveInternalDofManager(int i) const { return bubble; }
    virtual void giveInternalDofManDofIDMask(int i, EquationID eid, IntArray &answer) const;

    virtual FEInterpolation *giveInterpolation();
    virtual FEInterpolation *giveInterpolation(DofIDItem id);

    virtual void giveDofManDofIDMask(int inode, EquationID eid, IntArray &answer) const;

    virtual void updateYourself(TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType it);

    // Spatial localizer interface:
    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    // Element interpolation interface:
    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType u,
            TimeStep *stepN, const FloatArray &coords, FloatArray &answer);
    virtual void EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType u,
            TimeStep *stepN, const FloatArray &coords, FloatArray &answer);
    virtual void EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer);

    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }
};
} // end namespace oofem
#endif // tr1bubblestokes_h
