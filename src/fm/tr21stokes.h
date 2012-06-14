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

#ifndef tr21stokes_h
#define tr21stokes_h

#include "fmelement.h"
#include "domain.h"

#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"

namespace oofem {

class FEI2dTrLin;
class FEI2dTrQuad;

/**
 * Triangular Taylor-Hood element for Stokes flow.
 * Quadratic interpolation of geometry and velocity, and linear interpolation of pressures.
 *
 * @author Carl Sandström
 * @author Mikael Öhman
 */
class Tr21Stokes : public FMElement,
    public NodalAveragingRecoveryModelInterface,
    public SpatialLocalizerInterface,
    public EIPrimaryUnknownMapperInterface
{
protected:
    /// Number of gauss points. Same for pressure and velocity.
    int numberOfGaussPoints;
    /// Interpolation for pressure
    static FEI2dTrLin interpolation_lin;
    /// Interpolation for geometry and velocity
    static FEI2dTrQuad interpolation_quad;
    /// Ordering of dofs in element. Used to assemble the element stiffness
    static IntArray ordering;
    /// Ordering of dofs on edges. Used to assemble edge loads
    static IntArray edge_ordering [ 3 ];
    /// Dummy variable
    static bool __initialized;
    /// Defines the ordering of the dofs in the local stiffness matrix.
    static bool initOrdering() {
        ordering.setValues(15,  1, 2, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 3, 6, 9);
        edge_ordering [ 0 ].setValues(6,  1, 2, 4, 5, 10, 11);
        edge_ordering [ 1 ].setValues(6,  4, 5, 7, 8, 12, 13);
        edge_ordering [ 2 ].setValues(6,  7, 8, 1, 2, 14, 15);
        return true;
    }

public:
    Tr21Stokes(int n, Domain *d);
    virtual ~Tr21Stokes();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void computeGaussPoints();
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);

    void computeInternalForcesVector(FloatArray &answer, TimeStep *tStep);
    void computeLoadVector(FloatArray &answer, TimeStep *tStep);
    void computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep);
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep);

    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_2; }
    virtual const char *giveClassName() const { return "Tr21Stokes"; }
    virtual classType giveClassID() const { return Tr21StokesElementClass; }
    virtual MaterialMode giveMaterialMode() { return _2dFlow; }

    virtual int computeNumberOfDofs(EquationID ut);

    virtual FEInterpolation *giveInterpolation();
    virtual FEInterpolation *giveInterpolation(DofIDItem id);

    /**
     * Gives the dof ID mask for the element.
     * This element (Taylor-Hood) has V_u, V_v, P_f in node 1,2,3 and V_u, V_v in node 4,5,6.
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
    virtual double SpatialLocalizerI_giveClosestPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &gcoords);

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
#endif // tr21stokes_h
