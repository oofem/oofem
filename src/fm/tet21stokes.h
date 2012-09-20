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

#ifndef tet21stokes_h
#define tet21stokes_h

#include "fmelement.h"
#include "domain.h"

#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"

namespace oofem {

class FEI3dTrLin;
class FEI3dTetQuad;

/**
 * Tetrahedral Taylor-Hood element for Stokes flow.
 * Quadratic interpolation of geometry and velocity, and linear interpolation of pressures.
 * @author Mikael Öhman
 */
class Tet21Stokes : public FMElement,
    public NodalAveragingRecoveryModelInterface,
    public SpatialLocalizerInterface,
    public EIPrimaryUnknownMapperInterface
{
protected:
    /// Number of gauss points. Same for pressure and velocity.
    int numberOfGaussPoints;
    /// Interpolation for pressure
    static FEI3dTrLin interpolation_lin;
    /// Interpolation for geometry and velocity
    static FEI3dTetQuad interpolation_quad;
    /// Ordering of momentum balance dofs in element. Used to assemble the element stiffness
    static IntArray momentum_ordering;
    /// Ordering of conservation dofs in element. Used to assemble the element stiffness
    static IntArray conservation_ordering;
    /// Ordering of dofs on surfaces. Used to assemble edge loads (only momentum balance)
    static IntArray surf_ordering [ 4 ];

    /// Dummy variable
    static bool __initialized;

    /// Defines the ordering of the dofs in the local stiffness matrix.
    static bool initOrdering() {
        for (int i = 0, j = 1; i < 10; ++i) {
            momentum_ordering(i*3+0) = j++;
            momentum_ordering(i*3+1) = j++;
            momentum_ordering(i*3+2) = j++;
            j++;
        }
        conservation_ordering.setValues(4, 4, 8, 12, 16);

        surf_ordering [ 0 ].setValues(18, 1,  2,  3,  // node 1
                                          9, 10, 11,  // node 3
                                          5,  6,  7,  // node 2
                                         23, 24, 25,  // node 7
                                         20, 21, 22,  // node 6
                                         17, 18, 19); // node 5

        surf_ordering [ 1 ].setValues(18, 1,  2,  3,  // node 1
                                          5,  6,  7,  // node 2
                                         13, 14, 15,  // node 4
                                         17, 18, 19,  // node 5
                                         29, 30, 31,  // node 9
                                         26, 27, 28); // node 8

        surf_ordering [ 2 ].setValues(18, 5,  6,  7,  // node 2
                                          9, 10, 11,  // node 3
                                         13, 14, 15,  // node 4
                                         20, 21, 22,  // node 6
                                         32, 33, 34,  // node 10
                                         29, 30, 31); // node 9

        surf_ordering [ 2 ].setValues(18, 1,  2,  3,  // node 1
                                         13, 14, 15,  // node 4
                                          9, 10, 11,  // node 3
                                         26, 27, 28,  // node 8
                                         32, 33, 34,  // node 10
                                         23, 24, 25); // node 7
        return true;
    }

public:
    Tet21Stokes(int n, Domain *d);
    virtual ~Tet21Stokes();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void computeGaussPoints();
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);

    virtual void computeInternalForcesVector(FloatArray &answer, TimeStep *tStep);
    virtual void computeLoadVector(FloatArray &answer, TimeStep *tStep);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeSurfaceBCSubVectorAt(FloatArray &answer, Load *load, int iSurf, TimeStep *tStep);
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep);

    virtual Element_Geometry_Type giveGeometryType() const { return EGT_tetra_2; }
    virtual const char *giveClassName() const { return "Tet21Stokes"; }
    virtual classType giveClassID() const { return Tet21StokesElementClass; }
    virtual MaterialMode giveMaterialMode() { return _3dFlow; }

    virtual int computeNumberOfDofs(EquationID ut);

    virtual FEInterpolation *giveInterpolation();
    virtual FEInterpolation *giveInterpolation(DofIDItem id);

    /**
     * Gives the dof ID mask for the element.
     * This element (Taylor-Hood) has V_u, V_v, P_f in node corner and V_u, V_v in edge nodes.
     * @param inode Node to check.
     * @param ut Equation ID to check.
     * @param answer List of dof IDs.
     */
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;

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

    // Nodal averaging interface:
    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type);
};
} // end namespace oofem
#endif // tet21stokes_h
