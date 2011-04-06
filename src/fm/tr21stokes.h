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

#ifndef tr21stokes_h
#define tr21stokes_h

#include "fmelement.h"
#include "domain.h"
#include "fei2dtrlin.h"
#include "fei2dtrquad.h"

#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"

/**
 * Triangular Taylor-Hood element for Stokes flow.
 * Quadratic interpolation of geometry and velocity, and linear interpolation of pressures.
 * @author Carl Sandström
 * @author Mikael Öhman
 */
namespace oofem {
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
        //1 2 4 5 7 8 10 11 12 13 14 15 3 6 9
        ordering.at(1)  = 1;
        ordering.at(2)  = 2;
        ordering.at(3) = 4;
        ordering.at(4)  = 5;
        ordering.at(5)  = 7;
        ordering.at(6) = 8;
        ordering.at(7)  = 10;
        ordering.at(8)  = 11;
        ordering.at(9) = 12;
        ordering.at(10) = 13;
        ordering.at(11) = 14;
        ordering.at(12) = 15;
        ordering.at(13) = 3;
        ordering.at(14) = 6;
        ordering.at(15) = 9;

        //1, 2, 4, 5, 10, 11
        //4, 5, 7, 8, 12, 13
        //7, 8, 1, 2, 14, 15
        edge_ordering [ 0 ](0) = 1;
        edge_ordering [ 0 ](1) = 2;
        edge_ordering [ 0 ](2) = 4;
        edge_ordering [ 0 ](3) = 5;
        edge_ordering [ 0 ](4) = 10;
        edge_ordering [ 0 ](5) = 11;

        edge_ordering [ 1 ](0) = 4;
        edge_ordering [ 1 ](1) = 5;
        edge_ordering [ 1 ](2) = 7;
        edge_ordering [ 1 ](3) = 8;
        edge_ordering [ 1 ](4) = 12;
        edge_ordering [ 1 ](5) = 13;

        edge_ordering [ 2 ](0) = 7;
        edge_ordering [ 2 ](1) = 8;
        edge_ordering [ 2 ](2) = 1;
        edge_ordering [ 2 ](3) = 2;
        edge_ordering [ 2 ](4) = 14;
        edge_ordering [ 2 ](5) = 15;
        return true;
    }

public:
    Tr21Stokes(int, Domain *);
    ~Tr21Stokes();

    IRResultType initializeFrom(InputRecord *);

    void computeGaussPoints();
    void giveCharacteristicVector(FloatArray & answer, CharType, ValueModeType, TimeStep *);
    void giveCharacteristicMatrix(FloatMatrix & answer, CharType, TimeStep *);

    void computeInternalForcesVector(FloatArray &answer, TimeStep *tStep);
    void computeLoadVector(FloatArray &answer, TimeStep *atTime);
    void computeStiffnessMatrix(FloatMatrix &answer, TimeStep *atTime);
    void computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep);
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep);

    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_2; }
    const char *giveClassName() const { return "Tr21Stokes"; }
    classType giveClassID() const { return Tr21StokesElementClass; }

    /// @see Element::computeNumberOfDofs
    virtual int computeNumberOfDofs(EquationID ut);

    /**
     * Gives the dof ID mask for the element.
     * This element (Taylor-Hood) has V_u, V_v, P_f in node 1,2,3 and V_u, V_v in node 4,5,6.
     * @param inode Node to check.
     * @param ut Equation ID to check.
     * @param answer List of dof IDs.
     */
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;

    void updateYourself(TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType it);

    virtual double computeArea() const;

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual int computeLocalCoordinates(FloatArray &lcoords, const FloatArray &coords);

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
