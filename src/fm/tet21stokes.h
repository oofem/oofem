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
        conservation_ordering(0) = 4;
        conservation_ordering(1) = 8;
        conservation_ordering(2) = 12;
        conservation_ordering(3) = 16;

        surf_ordering [ 0 ](0)  = 1;  surf_ordering [ 0 ](1)  = 2;  surf_ordering [ 0 ](2)  = 3;  // node 1
        surf_ordering [ 0 ](3)  = 9;  surf_ordering [ 0 ](4)  = 10; surf_ordering [ 0 ](5)  = 11; // node 3
        surf_ordering [ 0 ](6)  = 5;  surf_ordering [ 0 ](7)  = 6;  surf_ordering [ 0 ](8)  = 7;  // node 2
        surf_ordering [ 0 ](9)  = 23; surf_ordering [ 0 ](10) = 24; surf_ordering [ 0 ](11) = 25; // node 7
        surf_ordering [ 0 ](12) = 20; surf_ordering [ 0 ](13) = 21; surf_ordering [ 0 ](14) = 22; // node 6
        surf_ordering [ 0 ](15) = 17; surf_ordering [ 0 ](16) = 18; surf_ordering [ 0 ](17) = 19; // node 5

        surf_ordering [ 1 ](0)  = 1;  surf_ordering [ 1 ](1)  = 2;  surf_ordering [ 1 ](2)  = 3;  // node 1
        surf_ordering [ 1 ](3)  = 5;  surf_ordering [ 1 ](4)  = 6;  surf_ordering [ 1 ](5)  = 7;  // node 2
        surf_ordering [ 1 ](6)  = 13; surf_ordering [ 1 ](7)  = 14; surf_ordering [ 1 ](8)  = 15; // node 4
        surf_ordering [ 1 ](9)  = 17; surf_ordering [ 1 ](10) = 18; surf_ordering [ 1 ](11) = 19; // node 5
        surf_ordering [ 1 ](12) = 29; surf_ordering [ 1 ](13) = 30; surf_ordering [ 1 ](14) = 31; // node 9
        surf_ordering [ 1 ](15) = 26; surf_ordering [ 1 ](16) = 27; surf_ordering [ 1 ](17) = 28; // node 8

        surf_ordering [ 2 ](0)  = 5;  surf_ordering [ 2 ](1)  = 6;  surf_ordering [ 2 ](2)  = 7;  // node 2
        surf_ordering [ 2 ](3)  = 9;  surf_ordering [ 2 ](4)  = 10; surf_ordering [ 2 ](5)  = 11; // node 3
        surf_ordering [ 2 ](6)  = 13; surf_ordering [ 2 ](7)  = 14; surf_ordering [ 2 ](8)  = 15; // node 4
        surf_ordering [ 2 ](9)  = 20; surf_ordering [ 2 ](10) = 21; surf_ordering [ 2 ](11) = 22; // node 6
        surf_ordering [ 2 ](12) = 32; surf_ordering [ 2 ](13) = 33; surf_ordering [ 2 ](14) = 34; // node 10
        surf_ordering [ 2 ](15) = 29; surf_ordering [ 2 ](16) = 30; surf_ordering [ 2 ](17) = 31; // node 9

        surf_ordering [ 3 ](0)  = 1;  surf_ordering [ 3 ](1)  = 2;  surf_ordering [ 3 ](2)  = 3;  // node 1
        surf_ordering [ 3 ](3)  = 13; surf_ordering [ 3 ](4)  = 14; surf_ordering [ 3 ](5)  = 15; // node 4
        surf_ordering [ 3 ](6)  = 9;  surf_ordering [ 3 ](7)  = 10; surf_ordering [ 3 ](8)  = 11; // node 3
        surf_ordering [ 3 ](9)  = 26; surf_ordering [ 3 ](10) = 27; surf_ordering [ 3 ](11) = 28; // node 8
        surf_ordering [ 3 ](12) = 32; surf_ordering [ 3 ](13) = 33; surf_ordering [ 3 ](14) = 34; // node 10
        surf_ordering [ 3 ](15) = 23; surf_ordering [ 3 ](16) = 24; surf_ordering [ 3 ](17) = 25; // node 7
        return true;
    }

public:
    Tet21Stokes(int n, Domain *d);
    ~Tet21Stokes();

    IRResultType initializeFrom(InputRecord *ir);

    void computeGaussPoints();
    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);

    void computeInternalForcesVector(FloatArray &answer, TimeStep *tStep);
    void computeLoadVector(FloatArray &answer, TimeStep *tStep);
    void computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeSurfaceBCSubVectorAt(FloatArray &answer, Load *load, int iSurf, TimeStep *tStep);
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep);

    Element_Geometry_Type giveGeometryType() const { return EGT_tetra_2; }
    const char *giveClassName() const { return "Tet21Stokes"; }
    classType giveClassID() const { return Tet21StokesElementClass; }
    MaterialMode giveMaterialMode() { return _3dFlow; }

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

    void updateYourself(TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType it);

    virtual double computeVolume() const;

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
#endif // tet21stokes_h
