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
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef graddpelement_h
#define graddpelement_h

#include "structuralelement.h"
#include "nlstructuralelement.h"

namespace oofem {
/**
 * Abstract class for gradient formulation of coupled damage-plasticity model(GradDp).
 * Yield function is formulated in the effective stress space and damage is driven by the nonlocal(over-nonlocal) cumulated plastic strain.
 * The new nonlocal degrees of freedom (with the meaning of the nonlocal cumulated plastic strain) are introduced with lower order of approximation functions than the displacement field to avoid spurious stress oscillations
 */
class GradDpElement
{
protected:
    int numberOfGaussPoints;
    int nPrimNodes, nPrimVars, nSecNodes, nSecVars;
    IntArray locU, locK;
    int totalSize, nlSize, locSize;
    int averType, nlGeo;

public:
    GradDpElement();
    virtual ~GradDpElement() {}

    virtual IRResultType initializeFrom(InputRecord *ir);
    //virtual void giveDofManDofIDMask (int inode, EquationID ut, IntArray& answer) const;

    // definition & identification
    virtual const char *giveClassName() const { return "GradDpElement"; }
    virtual classType giveClassID() const { return GradDpElementClass; }
    /***********************Predelat************************************/
    virtual int getNprimNodes() { return 0; }
    virtual int getNprimVars() { return 0; }
    virtual int getNsecNodes() { return 0; }
    virtual int getNsecVars() { return 0; }

    /************************************************************/
    //  virtual int  computeNumberOfDofs (EquationID ut){return (nPrimNodes*nPrimVars+nSecNodes*nSecVars);};
    //virtual int    checkConsistency();

protected:
    virtual StructuralElement *giveStructuralElement() = 0;
    virtual NLStructuralElement *giveNLStructuralElement() = 0;

    virtual void computeNkappaMatrixAt(GaussPoint *, FloatMatrix &) = 0;
    virtual void computeBkappaMatrixAt(GaussPoint *, FloatMatrix &) = 0;
    //void initialize();
    void setDisplacementLocationArray(IntArray &answer, int nPrimNodes, int nPrimVars, int nSecNodes, int nSecVars);
    void setNonlocalLocationArray(IntArray &answer, int nPrimNodes, int nPrimVars, int nSecNodes, int nSecVars);

    void computeStiffnessMatrix(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_uu(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_uk(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_kk(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_ku(FloatMatrix &, MatResponseMode, TimeStep *);

    void computeDisplacementDegreesOfFreedom(FloatArray &answer, TimeStep *stepN);
    void computeNonlocalDegreesOfFreedom(FloatArray &answer, TimeStep *stepN);
    void computeNonlocalGradient(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);

    void computeLocalStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
    void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
    void computeNonlocalCumulatedStrain(double &answer, GaussPoint *gp, TimeStep *stepN);

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveLocalInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveNonlocalInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);

    void computeStressVectorAndLocalCumulatedStrain(FloatArray &answer, double localCumulatedPlasticStrain, GaussPoint *gp, TimeStep *stepN);
    void computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    void computeLocForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    void computeLocNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    void computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    void computeDistanceToBoundary();

    //@todo this method has been removed from NLStructuralElement and has simply been put here to
    // not affect these derived elements
    virtual void computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, int i)
    { OOFEM_ERROR("computeNLBMatrixAt is not implemented"); }
};
} // end namespace oofem

#endif
