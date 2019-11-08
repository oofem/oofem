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

//   ***************************************************************************************
//   *** 2D LINEAR TRIANGULAR ELEMENT FOR FLUID DYNAMIC PROBLEMS SOLVED WITH PFEM METHOD ***
//   ***************************************************************************************

#ifndef tr1_2d_pfem_h
#define tr1_2d_pfem_h


#include "pfemelement2d.h"
#include "femcmpnn.h"
#include "domain.h"
#include "floatmatrix.h"
#include "fei2dtrlin.h"

///@name Input fields for TR1PFEM
//@{
#define _IFT_TR1_2D_PFEM_Name "tr1pfem"
//@}

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * This class is the implementation of triangular PFEM element with linear (and equal order) interpolation of velocity and pressure fields.
 *
 * @author David Krybus
 */
class TR1_2D_PFEM : public PFEMElement2d
{
protected:
    double b [ 3 ];
    double c [ 3 ];
    double area;

    /// Interpolation for velocity unknowns
    static FEI2dTrLin velocityInterpolation;
    /// Interpolation for pressure unknowns
    static FEI2dTrLin pressureInterpolation;

    static IntArray edge_ordering [ 3 ];

    /// Mask of velocity Dofs
    static IntArray velocityDofMask;
    /// Mask of pressure Dofs
    static IntArray pressureDofMask;

public:
    /**
     * Constructor of TR1_2D_PFEM Element. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     * @param particle1 number of first PFEMParticle
     * @param particle2 number of second PFEMParticle
     * @param particle3 number of third PFEMParticle
     * @param mat number of associated Material
     * @param cs number of associated CrossSection
     */
    TR1_2D_PFEM(int n, Domain *aDomain, int particle1, int particle2, int particle3, int mat, int cs);
    /// Destructor
    virtual ~TR1_2D_PFEM();

    void computeDiagonalMassMtrx(FloatArray &answer, TimeStep *tStep) override;
    void computeDiagonalMassMtrx(FloatMatrix &answer, TimeStep *tStep) override;

    double computeCriticalTimeStep(TimeStep *tStep) override;

    const char *giveClassName() const override { return "TR1_2D_PFEM"; }
    const char *giveInputRecordName() const override { return _IFT_TR1_2D_PFEM_Name; }

    Element_Geometry_Type giveGeometryType() const override { return EGT_triangle_1; }

    void giveElementDofIDMask(IntArray &answer) const override;

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    int computeNumberOfDofs() override;
    void initializeFrom(InputRecord &ir) override;
    int checkConsistency() override;

    double computeVolumeAround(GaussPoint *gp) override;

    Interface *giveInterface(InterfaceType) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    Element *giveElement() override { return this; }

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *atTime) override;
    // Graphics output
    //void drawYourself(oofegGraphicContext&) override;
    void drawRawGeometry(oofegGraphicContext &) override;
    void drawScalar(oofegGraphicContext &context) override;
    //void drawDeformedGeometry(oofegGraphicContext&, UnknownType) override {}
#endif

    FEInterpolation *giveVelocityInterpolation() override { return & velocityInterpolation; }
    FEInterpolation *givePressureInterpolation() override { return & pressureInterpolation; }

    FEInterpolation *giveInterpolation() const override { return & velocityInterpolation; }
    FEInterpolation *giveInterpolation(DofIDItem id) const override { return id == P_f ? & pressureInterpolation : & velocityInterpolation; }

    const IntArray &giveVelocityDofMask() const override
    {
        return this->velocityDofMask;
    }
    const IntArray &givePressureDofMask() const override
    {
        return this->pressureDofMask;
    }

protected:
    void computeGaussPoints() override;
    void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;

    void computeDeviatoricStressDivergence(FloatArray &answer, TimeStep *tStep) override;

    void computeBodyLoadVectorAt(FloatArray &answer, BodyLoad *load, TimeStep *tStep, ValueModeType mode) override;
    void computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // tr1_2d_pfem_h
