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
    ~TR1_2D_PFEM();

    virtual void computeDiagonalMassMtrx(FloatArray &answer, TimeStep *);
    virtual void computeDiagonalMassMtrx(FloatMatrix &answer, TimeStep *);

    virtual double computeCriticalTimeStep(TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "PFEMElement"; }
    virtual const char *giveInputRecordName() const { return _IFT_TR1_2D_PFEM_Name; }

    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

    virtual void giveElementDofIDMask(IntArray &answer) const;

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual int computeNumberOfDofs();
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void updateYourself(TimeStep *tStep);
    virtual int checkConsistency();

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual Interface *giveInterface(InterfaceType);

    virtual Element *giveElement() { return this; }

#ifdef __OOFEG
    virtual int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int node, TimeStep *atTime);
    //
    // Graphics output
    //
    //virtual void drawYourself (oofegGraphicContext&);
    virtual void  drawRawGeometry(oofegGraphicContext &);
    virtual void  drawScalar(oofegGraphicContext &context);
    //virtual void  drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

    /** Prints output of receiver to stream, for given time step */
    virtual void printOutputAt(FILE *, TimeStep *);

    virtual FEInterpolation *giveVelocityInterpolation() { return & velocityInterpolation; }
    virtual FEInterpolation *givePressureInterpolation() { return & pressureInterpolation; }

    virtual FEInterpolation *giveInterpolation() const { return & velocityInterpolation; }
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const { return id == P_f ? & pressureInterpolation : & velocityInterpolation; }

    virtual const IntArray &giveVelocityDofMask() const
    {
        return this->velocityDofMask;
    }
    virtual const IntArray &givePressureDofMask() const
    {
        return this->pressureDofMask;
    }

protected:
    virtual void computeGaussPoints();
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *);

    virtual void computeDeviatoricStressDivergence(FloatArray &answer, TimeStep *atTime);

    virtual void computeForceVector(FloatArray &answer, TimeStep *atTime); //F

    /// Calculates the body load vector
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *atTime);
    /// Calculates the boundary condition sub-vector on an edge
    virtual void computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *atTime);
};
} // end namespace oofem
#endif // tr1_2d_pfem_h
