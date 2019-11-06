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

#include "sm/Elements/structuralelement.h"

///@name Input fields for lumped mass element
//@{
#define _IFT_LumpedMassElement_Name "lumpedmass"
#define _IFT_LumpedMassElement_components "components"
#define _IFT_LumpedMassElement_dofs "dofs"
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
    /// Mass and moments of inertia corresponding to the specified DOFs
    FloatArray components;
    /// List of DOFs
    IntArray dofs;

public:
    LumpedMassElement(int n, Domain * d);
    virtual ~LumpedMassElement() { }

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { computeLumpedMassMatrix(answer, tStep); }
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override
    { answer.clear(); }
    void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { answer.clear(); }
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override
    { answer.clear(); }

    int computeNumberOfDofs() override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    void updateInternalState(TimeStep *tStep) override { }
    void updateYourself(TimeStep *tStep) override { }
    int checkConsistency() override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_LumpedMassElement_Name; }
    const char *giveClassName() const override { return "LumpedMassElement"; }
    void initializeFrom(InputRecord &ir) override;
    Element_Geometry_Type giveGeometryType() const override { return EGT_point; }

protected:
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override
    { answer.clear(); }
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override
    { answer.clear(); }

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                          int lowerIndx = 1, int upperIndx = ALL_STRAINS) override
    { answer.clear(); }
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override { answer.clear(); }
};
} // end namespace oofem
#endif // lumpedmasselement_h
