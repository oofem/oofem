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

#ifndef springelement_h
#define springelement_h

#include "../sm/Elements/structuralelement.h"

///@name Input fields for spring element
//@{
#define _IFT_SpringElement_Name "spring"
#define _IFT_SpringElement_mode "mode"
#define _IFT_SpringElement_orientation "orientation"
#define _IFT_SpringElement_springConstant "k"
#define _IFT_SpringElement_mass "m"

//@}

namespace oofem {
/**
 * This class implements a simple spring element. Its purpose is to introduce
 * a longitudinal or torsional spring between two nodes. The spring element is defined
 * by its spring constant and orientation.
 */
class SpringElement : public StructuralElement
{
public:
    /// Defines type of spring element (longitudinal/rotational) spring.
    enum SpringElementType {
        SE_1D_SPRING = 0,             ///< 1D spring element along x-axis.
        SE_2D_SPRING_XY = 1,          ///< 2D spring element in xy plane, requires D_u and D_v DOFs in each node (orientation vector should be in this plane). 
        SE_2D_SPRING_XZ = 2,          /// < 2D spring element in xz plane, requires D_u and D_w DOFs in each node (orientation vector should be in this plane). 
        SE_2D_TORSIONALSPRING_XZ = 3, ///< 2D torsional spring element in xz plane, requires R_v DOFs in each node.
        SE_3D_SPRING = 4,             ///< 3D spring element in space, requires D_u, D_v, and D_w DOFs in each node.
        SE_3D_TORSIONALSPRING = 5     ///< 3D torsional spring in space, requires R_u, R_v, and R_w DOFs in each node.
    };

protected:
    /// The longitudinal spring constant [Force/Length], torsional spring constant [Force*Length/Radians].
    double springConstant;
    /// total mass of the spring; to be distributed to nodes
    double mass;
    /**
     * Orientation vector. Defines orientation of spring element- for spring it defines the direction of spring,
     * for torsional spring it defines the axis of rotation.
     */
    FloatArray dir;
    /// Mode.
    SpringElementType mode;

public:
    SpringElement(int n, Domain * d);
    virtual ~SpringElement() { }

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
    { answer.clear(); }
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType)
    { answer.clear(); }
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    virtual int computeNumberOfDofs() { return 2; }
    virtual int computeNumberOfGlobalDofs();

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;

    virtual void updateInternalState(TimeStep *tStep) { }
    virtual void updateYourself(TimeStep *tStep) { }
    virtual int checkConsistency() { return 1; }
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

#ifdef __OOFEG
    //void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    //void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
#endif

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_SpringElement_Name; }
    virtual const char *giveClassName() const { return "SpringElement"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_point; }

protected:
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
    { answer.clear(); }

    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                                  int lowerIndx = 1, int upperIndx = ALL_STRAINS)
    { answer.clear(); }
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) { answer.clear(); }
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    double computeSpringInternalForce(TimeStep *tStep);
};
} // end namespace oofem
#endif // springelement_h
