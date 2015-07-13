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

#ifndef j2mplasticmaterial_h
#define j2mplasticmaterial_h

#include "mplasticmaterial.h"

///@name Input fields for J2MPlasticMaterial
//@{
#define _IFT_J2MPlasticMaterial_Name "j2mmat"
#define _IFT_J2MPlasticMaterial_ry "ry"
#define _IFT_J2MPlasticMaterial_khm "khm"
#define _IFT_J2MPlasticMaterial_ihm "ihm"
#define _IFT_J2MPlasticMaterial_rma "rma"
//@}

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
    J2MPlasticMaterial(int n, Domain * d);
    virtual ~J2MPlasticMaterial();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveInputRecordName() const { return _IFT_J2MPlasticMaterial_Name; }
    virtual const char *giveClassName() const { return "J2MPlasticMaterial"; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    virtual void computeStressSpaceHardeningVars(FloatArray &answer, GaussPoint *gp,
                                                 const FloatArray &strainSpaceHardeningVariables);
    virtual double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                                       const FloatArray &stressSpaceHardeningVars);
    virtual void computeHardeningReducedModuli(FloatMatrix &answer, GaussPoint *gp,
                                               const FloatArray &strainSpaceHardeningVariables,
                                               TimeStep *tStep);
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
                                        TimeStep *tStep);

    // auxiliary function
    double computeJ2InvariantAt(const FloatArray &stressVector);
    int giveSizeOfFullHardeningVarsVector();
    int giveSizeOfReducedHardeningVarsVector(GaussPoint *gp) const;
    double giveIsotropicHardeningVar(const FloatArray &stressSpaceHardeningVars);
    void giveStressBackVector(FloatArray &answer,
                              const FloatArray &stressSpaceHardeningVars);
};
} // end namespace oofem
#endif // j2mplasticmaterial_h
