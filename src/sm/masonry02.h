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

#ifndef masonry02_h
#define masonry02_h

#include "mplasticmaterial2.h"

namespace oofem {
class Domain;


/**
 * This class implements an interface masonry model based on
 * non associated multisurface plasticity.
 * Model follows the description from
 * Lourenco, P.B., Rots, J.G.: Multisurface Interface Model for Analysis of Masonry Structures
 * as published in Journal of Engineering Mechanics, vol 123, No. 7, 1997.
 */
class Masonry02 : public MPlasticMaterial2
{
protected:
    /// Tensile strength.
    double ft0;
    /// Mode I GF.
    double gfI;
    /// Mode II GF.
    double gfII;
    /// Residual friction angle.
    double tanfir;
    /// Initial friction angle.
    double tanfi0;
    /// Initial cohesion of joint.
    double c0;
    /// Cap mode parameters.
    double Cnn, Css, Cn;
    // double fm;

    /// Elastic properties.
    double kn, ks;

    /// Dilatancy angle.
    double tanpsi;

    /// Cap mode parameters.
    double sic, spc, smc, src;
    double kp, km, kr;
public:

    Masonry02(int n, Domain *d);
    virtual ~Masonry02();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual const char *giveClassName() const { return "Masonry02"; }
    virtual classType giveClassID() const { return Masonry02Class; }

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);
    virtual int giveStressStrainComponentIndOf(MatResponseForm, MaterialMode mmode, int);
    virtual void giveStressStrainMask(IntArray & answer, MatResponseForm, MaterialMode mmode) const;
    virtual int giveSizeOfReducedStressStrainVector(MaterialMode);
    void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &charVector3d);
    void giveFullCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &strainVector);

    virtual int giveSizeOfFullHardeningVarsVector() { return 3; }
    virtual int giveSizeOfReducedHardeningVarsVector(GaussPoint *) { return 3; }

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;


protected:
    virtual int giveMaxNumberOfActiveYieldConds(GaussPoint *gp) { return 2; }
    //
    // yield(YC-like functions) and loading(LC-like functions) criteria specific section
    //

    virtual double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                               const FloatArray &stressSpaceHardeningVars);

    virtual void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                     const FloatArray &stressSpaceHardeningVars);
    virtual void computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp,
                                                     const FloatArray &stress, const FloatArray &dlambda,
                                                     const FloatArray &dplasticStrain, const IntArray &activeConditionMap);
    virtual void computeKGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, FloatArray &fullStressVector,
                                        const FloatArray &strainSpaceHardeningVariables);

    virtual void computeReducedHardeningVarsSigmaGradient(FloatMatrix &answer, GaussPoint *gp, const IntArray &activeConditionMap,
                                                          const FloatArray &fullStressVector,
                                                          const FloatArray &strainSpaceHardeningVars,
                                                          const FloatArray &gamma);
    virtual void computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf,
                                                        const IntArray &activeConditionMap,
                                                        const FloatArray &fullStressVector,
                                                        const FloatArray &strainSpaceHardeningVars,
                                                        const FloatArray &gamma);
    virtual int hasHardening() { return 1; }

    virtual void computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                                const FloatArray &strainSpaceHardeningVariables);
    virtual void computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                                const FloatArray &strainSpaceHardeningVariables);


    void give2dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *atTime);

    virtual void computeReducedElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                             TimeStep *atTime);


    /// Cap mode related functions.
    double computeF3HardeningLaw(double k);
    double computeF3HardeningGradient(double k);
};
} // end namespace oofem
#endif // masonry02_h
