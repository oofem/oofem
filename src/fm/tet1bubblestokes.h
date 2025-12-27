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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef tet1bubblestokes_h
#define tet1bubblestokes_h

#include "fmelement.h"
#include "zznodalrecoverymodel.h"
#include "spatiallocalizer.h"
#include "elementinternaldofman.h"
#include "matresponsemode.h"

#include <memory>

#define _IFT_Tet1BubbleStokes_Name "tet1bubblestokes"

namespace oofem {
class FEI3dTetLin;

/**
 * Tetrahedral element for Stokes flow using Bubble basis function for stabilization.
 * Linear+Bubble interpolation of velocity, and linear interpolation of pressures.
 * The element is exported as a linear tetrahedron (bubble dofs are not exported).
 * It can deal with nonlinear material models, but it is assumed that the fluid is without memory (which is usually the case).
 *
 * @author Mikael Öhman
 */
class Tet1BubbleStokes : public FMElement,
public ZZNodalRecoveryModelInterface,
public SpatialLocalizerInterface
{
protected:
    /// Interpolation for pressure
    static FEI3dTetLin interp;
    /// Ordering of dofs in element. Used to assemble the element stiffness
    static IntArray momentum_ordering, conservation_ordering;
    /// Ordering of dofs on edges. Used to assemble edge loads
    static IntArray edge_ordering [ 6 ];
    /// Ordering of dofs on surfaces. Used to assemble surface loads
    static IntArray surf_ordering [ 4 ];

    /// The extra dofs from the bubble
    std :: unique_ptr< ElementDofManager > bubble;
    // Coordinates associated with the bubble dofs.
    //FloatArray bubbleCoord; // Assumed fixed at 0 for now (i.e. only linear geometry)

public:
    Tet1BubbleStokes(int n, Domain * d);

    double computeVolumeAround(GaussPoint *gp) override;

    void computeGaussPoints() override;

    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) override;
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) override;

    void computeInternalForcesVector(FloatArray &answer, TimeStep *tStep);
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep);

    void computeExternalForcesVector(FloatArray &answer, TimeStep *tStep);
    void computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep) override;
    void computeBoundarySurfaceLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global=true) override;

    const char *giveClassName() const override { return "Tet1BubbleStokes"; }
    const char *giveInputRecordName() const override { return _IFT_Tet1BubbleStokes_Name; }
    MaterialMode giveMaterialMode() override { return _3dFlow; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_tetra_1;}


    int computeNumberOfDofs() override;

    int giveNumberOfInternalDofManagers() const override { return 1; }
    DofManager *giveInternalDofManager(int i) const override { return bubble.get(); }
    void giveInternalDofManDofIDMask(int i, IntArray &answer) const override;

    FEInterpolation *giveInterpolation() const override;
    FEInterpolation *giveInterpolation(DofIDItem id) const override;

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    void updateYourself(TimeStep *tStep) override;

    Interface *giveInterface(InterfaceType it) override;

    void computeField(ValueModeType u, TimeStep *tStep, const FloatArray &coords, FloatArray &answer) override;
};
} // end namespace oofem
#endif // tet1bubblestokes_h
