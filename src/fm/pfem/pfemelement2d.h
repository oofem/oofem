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

//   *****************************************************************************
//   *** GENERAL 2D ELEMENT FOR FLUID DYNAMIC PROBLEMS SOLVED WITH PFEM METHOD ***
//   *****************************************************************************

#ifndef pfemelement2d_h
#define pfemelement2d_h


#include "pfemelement.h"
#include "femcmpnn.h"
#include "domain.h"
#include "floatmatrix.h"
#include "material.h"



#include "fei2dtrlin.h"

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * This class is the implementation of general 2d element with arbitrary interpolation of velocity and pressure fields.
 * Should be used with PFEM solution algorithm.
 *
 * @author David Krybus
 */
class PFEMElement2d : public PFEMElement
{
protected:

public:
    /// Constructor
    PFEMElement2d(int n, Domain *d);
    /// Destructor
    virtual ~PFEMElement2d();

    double computeCriticalTimeStep(TimeStep *tStep) override = 0;

    const char *giveClassName() const override { return "PFEMElement2d"; }
    Element_Geometry_Type giveGeometryType() const override { return EGT_triangle_1; }

    void giveElementDofIDMask(IntArray &answer) const override = 0;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override = 0;
    int computeNumberOfDofs() override = 0;
    void initializeFrom(InputRecord &ir, int priority) override;
    int checkConsistency() override;

    Interface *giveInterface(InterfaceType) override = 0;

    virtual Element *giveElement() { return this; }

#ifdef __OOFEG
    virtual int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int node, TimeStep *tStep) = 0;
    // Graphics output
    //virtual void drawYourself (oofegGraphicContext&);
    //virtual void drawRawGeometry(oofegGraphicContext &);
    //virtual void drawScalar(oofegGraphicContext &context);
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

    FEInterpolation *giveVelocityInterpolation() override = 0;
    FEInterpolation *givePressureInterpolation() override = 0;

protected:
    void computeGaussPoints() override = 0;
    void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override = 0;
    void computeDeviatoricStressDivergence(FloatArray &answer, TimeStep *tStep) override = 0;

    void computeBMatrix(FloatMatrix &answer, GaussPoint *gp) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep) override; //K
    void computePressureLaplacianMatrix(FloatMatrix &answer, TimeStep *tStep) override; //L
    void computeDivergenceMatrix(FloatMatrix &answerx, TimeStep *tStep) override; //D
    void computeGradientMatrix(FloatMatrix &answer, TimeStep *tStep) override; //G

    void computePrescribedRhsVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode) override;

    /// Calculates the shape function matrix on an edge
    void computeEdgeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    /// Calculates the shape function vector on an edge
    void computeEgdeNVectorAt(FloatArray &answer, int iedge, GaussPoint *gp);
    /// Calculates the volume around an edge
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    /// Gives the mapping for degrees of freedom on an edge
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
};
} // end namespace oofem
#endif // pfemelement_2d_h
