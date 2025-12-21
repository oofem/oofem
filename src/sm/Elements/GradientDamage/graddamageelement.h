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

#ifndef graddamageelement_h
#define graddamageelement_h

#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/nlstructuralelement.h"

#define _IFT_GradientDamageElement_penalty "penalty"


namespace oofem {
/**
 * Abstract class for gradient damage models
 * It can be also simply combined with damage-plastic models with yield function formulated in the effective stress space and damage driven by the nonlocal(over-nonlocal) cumulated plastic strain.
 * The new nonlocal degrees of freedom (they can, for example, have meaning of the nonlocal equivalent strain) are introduced. Usually, they are approximated using lower order approximation functions than the displacement field to avoid spurious stress oscillations.
 ***@author Martin Horak
 */




class GradientDamageElement
{
protected:
    int nPrimNodes, nPrimVars, nSecNodes, nSecVars;
    int totalSize, nlSize, locSize;
    IntArray locationArray_u;
    IntArray locationArray_d;
    double penalty;

public:
    GradientDamageElement();
    virtual ~GradientDamageElement() { }

    virtual void initializeFrom(InputRecord &ir);

protected:
    virtual StructuralElement *giveStructuralElement() = 0;
    virtual NLStructuralElement *giveNLStructuralElement() = 0;

    virtual void computeNdMatrixAt(GaussPoint *gp, FloatArray &answer) = 0;
    virtual void computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer) = 0;

    virtual void giveDofManDofIDMask_u(IntArray &answer) const = 0;
    virtual void giveDofManDofIDMask_d(IntArray &answer) const = 0;


    void computeStiffnessMatrix(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_uu(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_ud(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_dd(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_du(FloatMatrix &, MatResponseMode, TimeStep *);

    void computeDisplacementDegreesOfFreedom(FloatArray &answer, TimeStep *tStep);
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    void computeNonlocalDegreesOfFreedom(FloatArray &answer, TimeStep *tStep, ValueModeType vmt = VM_Total);
    void computeNonlocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep);
    void computeNonlocalDamageDrivingVariableGradient(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);


    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveInternalForcesVector_d(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);

    void computeInternalForces_dN(double &answer, double localDamageDrivingVariable, double damage, GaussPoint *gp, TimeStep *tStep);
    void computeInternalForces_dB(FloatArray &answer, double localDamageDrivingVariable, FloatArray damage_grad, GaussPoint *gp, TimeStep *tStep);

    void computeStressVector_and_localDamageDrivingVariable(FloatArray &answer, double &localCumulatedPlasticStrain, GaussPoint *gp, TimeStep *tStep);


    void giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_d, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u, const IntArray &dofIdArray_m);

    virtual void giveLocationArray_u(IntArray &answer) = 0;
    virtual void giveLocationArray_d(IntArray &answer) = 0;


    virtual const char *giveClassName() const { return "GradientDamageElement"; }

    virtual void postInitialize();
};
} // end namespace oofem

#endif // end gradient damage element
