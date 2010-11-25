/* $Header: /home/cvs/bp/oofem/sm/src/j2plasticmaterial.h,v 1.5 2003/04/06 14:08:30 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *******************************
//   *** CLASS J2 plastic material
//   *******************************

#ifndef j2plasticmaterial_h
#define j2plasticmaterial_h

#include "plasticmaterial.h"

namespace oofem {
class Domain;

class J2plasticMaterial : public PlasticMaterial
{
    /*
     * This class implements a isotropic  plastic linear material (J2 plasticity condition is used)
     * in a finite element problem. A material
     * is an attribute of a domain. It is usually also attribute of many elements.
     *
     * DESCRIPTION
     * ISOTROPIC PLASTIC Material with J2 plastic condition
     *
     * TASK
     * - Returning standard material stiffness marix for 3d-case.
     * according to current state determined by using data stored
     * in Gausspoint.
     * - Returning a material property (method 'give'). Only for non-standard elements.
     * - Returning real stress state vector(tensor) at gauss point for 3d - case.
     */

protected:
    int kinematicHardeningFlag, isotropicHardeningFlag;
    double kinematicModuli, isotropicModuli;
    //double E, nu; // isotropic material constants
    double k;
public:

    J2plasticMaterial(int n, Domain *d);
    ~J2plasticMaterial();

    IRResultType initializeFrom(InputRecord *ir);
    const char *giveClassName() const { return "J2plasticMaterial"; }
    classType giveClassID()         const { return PerfectlyPlasticMaterialClass; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const;


protected:

    //
    // yield(YC-like functions) and loading(LC-like functions) criteria specific section
    //

    virtual FloatArray *ComputeStressSpaceHardeningVars(GaussPoint *gp,
                                                        FloatArray *strainSpaceHardeningVariables);
    virtual double computeYieldValueAt(GaussPoint *gp, FloatArray *stressVector,
                                       FloatArray *stressSpaceHardeningVars);
    virtual void   computeHardeningReducedModuli(FloatMatrix &answer, GaussPoint *gp,
                                                 FloatArray *strainSpaceHardeningVariables,
                                                 TimeStep *atTime);
    virtual FloatArray *ComputeStressGradient(GaussPoint *gp, FloatArray *stressVector,
                                              FloatArray *stressSpaceHardeningVars);
    virtual FloatArray *ComputeStressSpaceHardeningVarsReducedGradient(GaussPoint *gp,
                                                                       FloatArray *stressVector,
                                                                       FloatArray *stressSpaceHardeningVars);
    virtual int hasHardening();
    virtual void  computeReducedGradientMatrix(FloatMatrix &answer, GaussPoint *gp,
                                               const FloatArray &stressVector,
                                               const FloatArray &stressSpaceHardeningVars);
    virtual void  computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                              const FloatArray &strainIncrement, TimeStep *atTime);
    virtual void  compute3dElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                         TimeStep *atTime);

    // auxiliary function
    double      computeJ2InvariantAt(FloatArray *);
    int         giveSizeOfFullHardeningVarsVector();
    int         giveSizeOfReducedHardeningVarsVector(GaussPoint *);
    double      giveIsotropicHardeningVar(FloatArray *stressSpaceHardeningVars);
    void        giveStressBackVector(FloatArray &answer,
                                     const FloatArray &stressSpaceHardeningVars);
};
} // end namespace oofem
#endif // j2plasticmaterial_h
