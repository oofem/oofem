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

#ifndef tet21stokes_h
#define tet21stokes_h

#include "fmelement.h"
#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"
#include "matresponsemode.h"

#define _IFT_Tet21Stokes_Name "tet21stokes"

namespace oofem {
class FEI3dTetLin;
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
    /// Interpolation for pressure
    static FEI3dTetLin interpolation_lin;
    /// Interpolation for geometry and velocity
    static FEI3dTetQuad interpolation_quad;
    /// Ordering of momentum balance dofs in element. Used to assemble the element stiffness
    static IntArray momentum_ordering;
    /// Ordering of conservation dofs in element. Used to assemble the element stiffness
    static IntArray conservation_ordering;
    /// Ordering of dofs on surfaces. Used to assemble edge loads (only momentum balance)
    static IntArray surf_ordering [ 4 ];

public:
    Tet21Stokes(int n, Domain * d);
    virtual ~Tet21Stokes();

    virtual void computeGaussPoints();
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);

    void computeInternalForcesVector(FloatArray &answer, TimeStep *tStep);
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep);

    void computeExternalForcesVector(FloatArray &answer, TimeStep *tStep);
    virtual void computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep);

    virtual const char *giveClassName() const { return "Tet21Stokes"; }
    virtual const char *giveInputRecordName() const { return _IFT_Tet21Stokes_Name; }
    virtual MaterialMode giveMaterialMode() { return _3dFlow; }

    virtual int computeNumberOfDofs();

    virtual FEInterpolation *giveInterpolation() const;
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const;

    /**
     * Gives the dof ID mask for the element.
     * This element (Taylor-Hood) has V_u, V_v, P_f in node corner and V_u, V_v in edge nodes.
     * @param inode Node to check.
     * @param answer List of dof IDs.
     */
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;

    virtual void updateYourself(TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType it);

    // Element interpolation interface:
    virtual void EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType u,
                                                                       TimeStep *tStep, const FloatArray &coords, FloatArray &answer);

    // Nodal averaging interface:
    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
};
} // end namespace oofem
#endif // tet21stokes_h
