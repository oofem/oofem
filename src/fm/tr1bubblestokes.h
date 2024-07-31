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

#ifndef tr1bubblestokes_h
#define tr1bubblestokes_h

#include "fmelement.h"
#include "zznodalrecoverymodel.h"
#include "spatiallocalizer.h"
#include "matresponsemode.h"

#define _IFT_Tr1BubbleStokes_Name "tr1bubblestokes"

namespace oofem {
class FEI2dTrLin;
class ElementDofManager;

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
public SpatialLocalizerInterface
{
protected:
    /// Interpolation for pressure
    static FEI2dTrLin interp;
    /// Ordering of dofs in element. Used to assemble the element stiffness
    static IntArray momentum_ordering, conservation_ordering;
    /// Ordering of dofs on edges. Used to assemble edge loads
    static IntArray edge_ordering [ 3 ];

    /// The extra dofs from the bubble
    std :: unique_ptr< ElementDofManager > bubble;
    // Coordinates associated with the bubble dofs.
    //FloatArray bubbleCoord; // Assumed fixed at 0 for now (i.e. only linear geometry)

public:
    Tr1BubbleStokes(int n, Domain * d);

    double computeVolumeAround(GaussPoint *gp) override;

    void computeGaussPoints() override;

    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) override;
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) override;

    void computeInternalForcesVector(FloatArray &answer, TimeStep *tStep);
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep);

    void computeExternalForcesVector(FloatArray &answer, TimeStep *tStep);
    void computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep) override;
    void computeBoundarySurfaceLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global=true) override;

    const char *giveClassName() const override { return "Tr1BubbleStokes"; }
    const char *giveInputRecordName() const override { return _IFT_Tr1BubbleStokes_Name; }
    MaterialMode giveMaterialMode() override { return _2dFlow; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_triangle_1;}

    int computeNumberOfDofs() override;

    int giveNumberOfInternalDofManagers() const override;
    DofManager *giveInternalDofManager(int i) const override;
    void giveInternalDofManDofIDMask(int i, IntArray &answer) const override;

    FEInterpolation *giveInterpolation() const override;
    FEInterpolation *giveInterpolation(DofIDItem id) const override;

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    void updateYourself(TimeStep *tStep) override;

    Interface *giveInterface(InterfaceType it) override;

    void computeField(ValueModeType u, TimeStep *tStep, const FloatArray &coords, FloatArray &answer) override;
};
} // end namespace oofem
#endif // tr1bubblestokes_h
