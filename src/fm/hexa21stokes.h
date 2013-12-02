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

#ifndef hexa21stokes_h
#define hexa21stokes_h

#include "fmelement.h"
#include "domain.h"
#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"

#define _IFT_Hexa21Stokes_Name "hexa21stokes"

namespace oofem {

class FEI3dHexaLin;
class FEI3dHexaTriQuad;

/**
 * Hexahedral Taylor-Hood element for Stokes flow.
 *(Tri)Quadratic interpolation of geometry and velocity, and linear interpolation of pressures.
 * @author Mikael Öhman
 */
class Hexa21Stokes : public FMElement,
    public NodalAveragingRecoveryModelInterface,
    public SpatialLocalizerInterface,
    public EIPrimaryUnknownMapperInterface
{
protected:
    /// Interpolation for pressure
    static FEI3dHexaLin interpolation_lin;
    /// Interpolation for geometry and velocity
    static FEI3dHexaTriQuad interpolation_quad;
    /// Ordering of momentum balance dofs in element. Used to assemble the element stiffness
    static IntArray momentum_ordering;
    /// Ordering of conservation dofs in element. Used to assemble the element stiffness
    static IntArray conservation_ordering;
    /// Ordering of dofs on surfaces. Used to assemble edge loads (only momentum balance)
    static IntArray surf_ordering [ 6 ];

    /// Dummy variable
    static bool __initialized;

    /// Defines the ordering of the dofs in the local stiffness matrix.
    static bool initOrdering() {
        for (int i = 0, j = 1; i < 27; ++i) {
            momentum_ordering(i*3+0) = j++;
            momentum_ordering(i*3+1) = j++;
            momentum_ordering(i*3+2) = j++;
            if ( i < 8 ) j++;
        }
        conservation_ordering.setValues(8, 4, 8, 12, 16, 20, 24, 28, 32);

        surf_ordering [ 0 ].setValues(27,
                                     5,   6,   7, // node 2
                                     1,   2,   3, // node 1
                                    13,  14,  15, // node 4
                                     9,  10,  11, // node 3
                                    33,  34,  35, // node 9
                                    42,  43,  44, // node 12
                                    39,  40,  41, // node 11
                                    36,  37,  38, // node 10
                                    69,  70,  71);// node 21

        surf_ordering [ 1 ].setValues(27,
                                    17,  18,  19, // node 5
                                    21,  22,  23, // node 6
                                    25,  26,  27, // node 7
                                    29,  30,  31, // node 8
                                    45,  46,  47, // node 13
                                    48,  49,  50, // node 14
                                    51,  52,  53, // node 15
                                    54,  55,  56, // node 16
                                    72,  73,  74);// node 22

        surf_ordering [ 2 ].setValues(27,
                                     1,   2,   3, // node 1
                                    17,  18,  19, // node 5
                                    21,  22,  23, // node 6
                                     5,   6,   7, // node 2
                                    57,  58,  59, // node 17
                                    45,  46,  47, // node 13
                                    60,  61,  62, // node 18
                                    33,  34,  35, // node 9
                                    75,  76,  77);// node 23

        surf_ordering [ 3 ].setValues(27,
                                     5,   6,   7, // node 2
                                     9,  10,  11, // node 3
                                    25,  26,  27, // node 7
                                    21,  22,  23, // node 6
                                    36,  37,  38, // node 10
                                    63,  64,  65, // node 19
                                    48,  49,  50, // node 14
                                    60,  61,  62, // node 18
                                    78,  79,  80);// node 24

        surf_ordering [ 4 ].setValues(27,
                                     9,  10,  11, // node 3
                                    13,  14,  15, // node 4
                                    29,  30,  31, // node 8
                                    25,  26,  27, // node 7
                                    39,  40,  41, // node 11
                                    66,  67,  68, // node 20
                                    51,  52,  53, // node 15
                                    63,  64,  65, // node 19
                                    81,  82,  83);// node 25

        surf_ordering [ 5 ].setValues(27,
                                    13,  14,  15, // node 4
                                     1,   2,   3, // node 1
                                    17,  18,  19, // node 5
                                    29,  30,  31, // node 8
                                    42,  43,  44, // node 12
                                    57,  58,  59, // node 17
                                    54,  55,  56, // node 16
                                    66,  67,  68, // node 20
                                    84,  85,  86);// node 26
        return true;
    }

public:
    Hexa21Stokes(int n, Domain *d);
    virtual ~Hexa21Stokes();

    virtual void computeGaussPoints();
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);

    void computeInternalForcesVector(FloatArray &answer, TimeStep *tStep);
    void computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep);

    void computeExternalForcesVector(FloatArray &answer, TimeStep *tStep);
    virtual void computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep);

    virtual const char *giveClassName() const { return "Hexa21Stokes"; }
    virtual const char *giveInputRecordName() const { return _IFT_Hexa21Stokes_Name; }
    virtual classType giveClassID() const { return Hexa21StokesElementClass; }
    virtual MaterialMode giveMaterialMode() { return _3dFlow; }

    virtual int computeNumberOfDofs();

    virtual FEInterpolation *giveInterpolation() const;
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const;

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
};
} // end namespace oofem
#endif // hexa21stokes_h
