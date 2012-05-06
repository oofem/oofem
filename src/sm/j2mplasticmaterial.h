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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef j2mplasticmaterial_h
#define j2mplasticmaterial_h

#include "mplasticmaterial.h"

namespace oofem {
class Domain;

/**
 * This class implements a isotropic  plastic linear material (J2 plasticity condition is used)
 * in a finite element problem.
 * Both kinematic and isotropic hardening is supported.
 */
class J2MPlasticMaterial : public MPlasticMaterial
{
protected:
    int kinematicHardeningFlag, isotropicHardeningFlag;
    double kinematicModuli, isotropicModuli;
    double k;

public:
    J2MPlasticMaterial(int n, Domain *d);
    virtual ~J2MPlasticMaterial();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "J2plasticMaterial"; }
    virtual classType giveClassID() const { return J2MPlasticMaterialClass; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    virtual void computeStressSpaceHardeningVars(FloatArray &answer, GaussPoint *gp,
                                                 const FloatArray &strainSpaceHardeningVariables);
    virtual double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                                       const FloatArray &stressSpaceHardeningVars);
    virtual void computeHardeningReducedModuli(FloatMatrix &answer, GaussPoint *gp,
                                               const FloatArray &strainSpaceHardeningVariables,
                                                 TimeStep *atTime);
    virtual void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                             const FloatArray &stressSpaceHardeningVars);
    virtual void computeStressSpaceHardeningVarsReducedGradient(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp,
                                                                const FloatArray &stressVector,
                                                                const FloatArray &stressSpaceHardeningVars);
    virtual int hasHardening();
    virtual void computeReducedGradientMatrix(FloatMatrix &answer, int isurf,
                                              GaussPoint *gp,
                                              const FloatArray &stressVector,
                                              const FloatArray &stressSpaceHardeningVars);
    virtual void compute3dElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                        TimeStep *atTime);

    // auxiliary function
    double computeJ2InvariantAt(const FloatArray &stressVector);
    int giveSizeOfFullHardeningVarsVector();
    int giveSizeOfReducedHardeningVarsVector(GaussPoint *gp);
    double giveIsotropicHardeningVar(const FloatArray &stressSpaceHardeningVars);
    void giveStressBackVector(FloatArray &answer,
                              const FloatArray &stressSpaceHardeningVars);
};
} // end namespace oofem
#endif // j2mplasticmaterial_h
