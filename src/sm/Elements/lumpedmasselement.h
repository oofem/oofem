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

#ifndef lumpedmasselement_h
#define lumpedmasselement_h

#include "../sm/Elements/structuralelement.h"

///@name Input fields for lumped mass element
//@{
#define _IFT_LumpedMassElement_Name "lumpedmass"
#define _IFT_LumpedMassElement_components "components"
//@}

namespace oofem {
/**
 * This class implements a simple lumped mass element. Its purpose is to introduce
 * an additional mass (mass components or rotary inertias) into a node.
 * The mass element is defined by a single node.
 * At present, mass is defined in the nodal coordinate system.
 * The same element can be used to add an additional stiffness if needed (Not yet implemented).
 */
class LumpedMassElement : public StructuralElement
{
protected:
    ///Mass and moments of inertia corresponding to nodal DOFs
    FloatArray components;

public:
    LumpedMassElement(int n, Domain * d);
    virtual ~LumpedMassElement() { }

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
    { answer.clear(); }
    virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
    { answer.clear(); }
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType)
    { answer.clear(); }
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0)
    { answer.clear(); }

    virtual int computeNumberOfDofs();
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;

    virtual void updateInternalState(TimeStep *tStep) { }
    virtual void updateYourself(TimeStep *tStep) { }
    virtual int checkConsistency();

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
#endif

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_LumpedMassElement_Name; }
    virtual const char *giveClassName() const { return "LumpedMassElement"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_point; }

protected:
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
    { answer.clear(); }

    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                                  int lowerIndx = 1, int upperIndx = ALL_STRAINS)
    { answer.clear(); }
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) { answer.clear(); }
};
} // end namespace oofem
#endif // lumpedmasselement_h
