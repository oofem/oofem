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

#include "j2mplasticmaterial.h"
#include "isolinearelasticmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"

#include "structuralcrosssection.h"
#include "mathfem.h"

namespace oofem {
J2MPlasticMaterial :: J2MPlasticMaterial(int n, Domain *d) : MPlasticMaterial(n, d)
{
    //
    // constructor
    //
    kinematicHardeningFlag = isotropicHardeningFlag = 0;
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    this->nsurf = 1;
}

J2MPlasticMaterial :: ~J2MPlasticMaterial()
{
    //
    // destructor
    //
}

IRResultType
J2MPlasticMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    double value;

    MPlasticMaterial :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir);

    IR_GIVE_FIELD(ir, value, IFT_J2MPlasticMaterial_ry, "ry"); // Macro
    k = value / sqrt(3.0);

    //  E = readDouble (initString,"e");
    // nu = readDouble (initString,"nu");
    kinematicModuli = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, kinematicModuli, IFT_J2MPlasticMaterial_khm, "khm"); // Macro

    isotropicModuli = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, isotropicModuli, IFT_J2MPlasticMaterial_ihm, "ihm"); // Macro

    if ( fabs(kinematicModuli) > 1.e-12 ) {
        kinematicHardeningFlag = 1;
    }

    if ( fabs(isotropicModuli) > 1.e-12 ) {
        isotropicHardeningFlag = 1;
    }

    int rma = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, rma, IFT_J2MPlasticMaterial_rma, "rma"); // Macro
    if ( rma == 0 ) {
        this->rmType = mpm_ClosestPoint;
    } else {
        this->rmType = mpm_CuttingPlane;
    }


    return IRRT_OK;
}


MaterialStatus *
J2MPlasticMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates new  material status  corresponding to this class
 */
{
    MPlasticMaterialStatus *status;

    status = new MPlasticMaterialStatus(1, this->giveDomain(), gp);
    return status;
}

void
J2MPlasticMaterial :: computeStressSpaceHardeningVars(FloatArray &answer, GaussPoint *gp,
                                                      const FloatArray &strainSpaceHardeningVariables)
{
    // in full stress strain space

    int i;
    int count = 0, size = this->giveSizeOfFullHardeningVarsVector(), isize, rSize;
    IntArray mask;

    if ( !hasHardening() ) {
        answer.resize(0);
        return;
    }

    answer.resize(size);
    this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode() );
    isize = mask.giveSize();
    rSize = this->giveSizeOfReducedHardeningVarsVector(gp);

    /* kinematic hardening variables are first */
    if ( this->kinematicHardeningFlag ) {
        for ( i = 1; i <= isize; i++ ) {
            // to be consistent with equivalent plastic strain formulation
            // we multiply by (sqrt(2.)*2./3.)
            answer.at( mask.at(i) ) = ( sqrt(2.) * 2. / 3. ) * this->kinematicModuli * strainSpaceHardeningVariables.at(i);
        }

        count = 6;
    }

    if ( this->isotropicHardeningFlag ) {
        answer.at(count + 1) = this->isotropicModuli *
                               strainSpaceHardeningVariables.at(rSize);
    }

    answer.negated();
}



double
J2MPlasticMaterial :: computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                                          const FloatArray &stressSpaceHardeningVars)
{
    double f;
    FloatArray helpVector, backStress;

    if ( this->kinematicHardeningFlag ) {
        if ( stressVector.isNotEmpty() ) {
            this->giveStressBackVector(backStress, stressSpaceHardeningVars);
            helpVector = stressVector;
            helpVector.add(backStress);
        } else {
            return -k;
        }
    } else {
        helpVector = stressVector;
    }

    f = sqrt( this->computeJ2InvariantAt(helpVector) );

    //if (this->kinematicHardeningFlag) delete helpVector;

    return f + sqrt(1. / 3.) * this->giveIsotropicHardeningVar(stressSpaceHardeningVars) - this->k;
}


void
J2MPlasticMaterial :: computeHardeningReducedModuli(FloatMatrix &answer, GaussPoint *gp,
                                                    const FloatArray &strainSpaceHardeningVariables,
                                                    TimeStep *atTime)
{
    /* computes hardening moduli in reduced stress strain space (for kinematic back-stress)*/

    int i;
    int size = this->giveSizeOfReducedHardeningVarsVector(gp);

    if ( !hasHardening() ) {
        answer.resize(0, 0);
        return;
    }

    answer.resize(size, size);
    answer.zero();

    /* kinematic hardening variables are first */
    if ( this->kinematicHardeningFlag ) {
        int ksize = this->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() );
        for ( i = 1; i <= ksize; i++ ) {
            answer.at(i, i) = this->kinematicModuli;
        }
    }

    if ( this->isotropicHardeningFlag ) {
        answer.at(size, size) = this->isotropicModuli;
    }
}


void
J2MPlasticMaterial :: computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                                  const FloatArray &stressSpaceHardeningVars)
{
    /* stress gradient of yield function in full stress - strain space */

    double f, ax, ay, az, sx, sy, sz;
    FloatArray helpVector, backStress;

    answer.resize(6);
    answer.zero();
    if ( this->kinematicHardeningFlag ) {
        if ( stressVector.isNotEmpty() ) {
            this->giveStressBackVector(backStress, stressSpaceHardeningVars);
            helpVector = stressVector;
            helpVector.add(backStress);
        } else {
            return;
        }
    } else {
        helpVector = stressVector;
    }

    f = sqrt( this->computeJ2InvariantAt(helpVector) );
    // check for yield value zero value
    if ( fabs(f) < 1.e-6 ) {
        return;
    }

    ax = helpVector.at(1);
    ay = helpVector.at(2);
    az = helpVector.at(3);

    sx = ( 2. / 3. ) * ax - ( 1. / 3. ) * ay - ( 1. / 3. ) * az;
    sy = ( 2. / 3. ) * ay - ( 1. / 3. ) * ax - ( 1. / 3. ) * az;
    sz = ( 2. / 3. ) * az - ( 1. / 3. ) * ay - ( 1. / 3. ) * ax;

    answer.at(1) = 0.5 * sx / f;
    answer.at(2) = 0.5 * sy / f;
    answer.at(3) = 0.5 * sz / f;
    answer.at(4) = helpVector.at(4) / f;
    answer.at(5) = helpVector.at(5) / f;
    answer.at(6) = helpVector.at(6) / f;
}

void
J2MPlasticMaterial :: computeStressSpaceHardeningVarsReducedGradient(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp,
                                                                     const FloatArray &stressVector,
                                                                     const FloatArray &stressSpaceHardeningVars)
{
    /* computes stress space hardening gradient in reduced stress-strain space */

    int i = 0, kcount = 0, size = this->giveSizeOfReducedHardeningVarsVector(gp);
    //double f,ax,ay,az,sx,sy,sz;
    FloatArray fullKinematicGradient, reducedKinematicGrad;
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );

    if ( !hasHardening() ) {
        answer.resize(0);
        return;
    }

    answer.resize(size);

    /* kinematic hardening variables first */
    if ( this->kinematicHardeningFlag ) {
        this->computeStressGradientVector(fullKinematicGradient, ftype, i, gp, stressVector, stressSpaceHardeningVars);
        crossSection->giveReducedCharacteristicVector(reducedKinematicGrad, gp, fullKinematicGradient);

        kcount = reducedKinematicGrad.giveSize();
    }

    if ( this->kinematicHardeningFlag ) {
        for ( i = 1; i <= kcount; i++ ) {
            answer.at(i) = reducedKinematicGrad.at(i);
        }
    }

    if ( this->isotropicHardeningFlag ) {
        answer.at(size) = sqrt(1. / 3.);
    }
}




int
J2MPlasticMaterial :: hasHardening()
{
    return ( this->kinematicHardeningFlag || this->isotropicHardeningFlag );
}


void
J2MPlasticMaterial :: computeReducedGradientMatrix(FloatMatrix &answer, int isurf,
                                                   GaussPoint *gp,
                                                   const FloatArray &stressVector,
                                                   const FloatArray &stressSpaceHardeningVars)
{
    int i, j, size;
    int imask, jmask;
    FloatArray helpVector, backStress, df(6);
    IntArray mask;
    double f, f32, f12, ax, ay, az;

    this->giveStressStrainMask( mask, FullForm, gp->giveMaterialMode() );
    size = giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) +
           this->giveSizeOfReducedHardeningVarsVector(gp);

    answer.resize(size, size);
    answer.zero();


    if ( stressVector.giveSize() != 0 ) {
        /* kinematic hardening variables first */
        if ( this->kinematicHardeningFlag ) {
            this->giveStressBackVector(backStress, stressSpaceHardeningVars);
            helpVector = stressVector;
            helpVector.add(backStress);
        } else {
            helpVector = stressVector;
        }

        f = this->computeJ2InvariantAt(helpVector);
        f12 = pow(f, 1. / 2.);
        f32 = pow(f, 3. / 2.);

        ax = helpVector.at(1);
        ay = helpVector.at(2);
        az = helpVector.at(3);

        df.at(1) = ( 2. / 3. ) * ax - ( 1. / 3. ) * ay - ( 1. / 3. ) * az;
        df.at(2) = ( 2. / 3. ) * ay - ( 1. / 3. ) * ax - ( 1. / 3. ) * az;
        df.at(3) = ( 2. / 3. ) * az - ( 1. / 3. ) * ay - ( 1. / 3. ) * ax;
        df.at(4) = 2. * helpVector.at(4);
        df.at(5) = 2. * helpVector.at(5);
        df.at(6) = 2. * helpVector.at(6);

        for ( i = 1; i <= 3; i++ ) {
            if ( ( imask = mask.at(i) ) == 0 ) {
                continue;
            }

            for ( j = i; j <= 3; j++ ) {
                if ( ( jmask = mask.at(j) ) == 0 ) {
                    continue;
                }

                if ( i == j ) {
                    answer.at(imask, jmask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j) + 0.5 * ( 1. / f12 ) * ( 4. / 6 );
                } else {
                    answer.at(imask, jmask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j) + 0.5 * ( 1. / f12 ) * ( -2. / 6 );
                    answer.at(jmask, imask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j) + 0.5 * ( 1. / f12 ) * ( -2. / 6 );
                }
            }
        }

        for ( i = 1; i <= 3; i++ ) {
            if ( ( imask = mask.at(i) ) == 0 ) {
                continue;
            }

            for ( j = 4; j <= 6; j++ ) {
                if ( ( jmask = mask.at(j) ) == 0 ) {
                    continue;
                }

                answer.at(imask, jmask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j);
                answer.at(jmask, imask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j);
            }
        }

        for ( i = 4; i <= 6; i++ ) {
            if ( ( imask = mask.at(i) ) == 0 ) {
                continue;
            }

            for ( j = i; j <= 6; j++ ) {
                if ( ( jmask = mask.at(j) ) == 0 ) {
                    continue;
                }

                if ( i == j ) {
                    answer.at(imask, jmask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j) + 0.5 * ( 1. / f12 ) * 2.;
                } else {
                    answer.at(imask, jmask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j);
                    answer.at(jmask, imask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j);
                }
            }
        }
    }

    /* for isotropic hardening: the corresponding part of gradient matrix is zero valued */
}


void
J2MPlasticMaterial :: compute3dElasticModuli(FloatMatrix &answer,
                                             GaussPoint *gp,
                                             TimeStep *atTime)
{
    /* Returns 3d elastic moduli */
    this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, FullForm, ElasticStiffness, gp, atTime);
}


double
J2MPlasticMaterial :: computeJ2InvariantAt(const FloatArray &stressVector)
{
    double answer;
    double v1, v2, v3;

    if ( stressVector.isEmpty() ) {
        return 0.0;
    }

    v1 = ( ( stressVector.at(1) - stressVector.at(2) ) * ( stressVector.at(1) - stressVector.at(2) ) );
    v2 = ( ( stressVector.at(2) - stressVector.at(3) ) * ( stressVector.at(2) - stressVector.at(3) ) );
    v3 = ( ( stressVector.at(3) - stressVector.at(1) ) * ( stressVector.at(3) - stressVector.at(1) ) );

    answer = ( 1. / 6. ) * ( v1 + v2 + v3 ) + stressVector.at(4) * stressVector.at(4) +
             stressVector.at(5) * stressVector.at(5) + stressVector.at(6) * stressVector.at(6);

    return answer;
}


int
J2MPlasticMaterial :: giveSizeOfFullHardeningVarsVector()
{
    /* Returns the size of hardening variables vector */
    int size = 0;

    if ( kinematicHardeningFlag ) {
        size += 6;                      /* size of full stress vector */
    }

    if ( isotropicHardeningFlag ) {
        size += 1;                      /* scalar value */
    }

    return size;
}

int
J2MPlasticMaterial :: giveSizeOfReducedHardeningVarsVector(GaussPoint *gp)
{
    /* Returns the size of hardening variables vector */
    int size = 0;

    if ( kinematicHardeningFlag ) {
        size += this->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() );
    }

    if ( isotropicHardeningFlag ) {
        size += 1;                      /* scalar value */
    }

    return size;
}


void
J2MPlasticMaterial :: giveStressBackVector(FloatArray &answer,
                                           const FloatArray &stressSpaceHardeningVars)
{
    /* returns part of hardening vector corresponding to kinematic hardening */
    if ( this->kinematicHardeningFlag ) {
        answer.resize(6);
        for ( int i = 1; i <= 6; i++ ) {
            answer.at(i) = stressSpaceHardeningVars.at(i);
        }

        return;
    }
    answer.resize(0);
}


double
J2MPlasticMaterial :: giveIsotropicHardeningVar(const FloatArray &stressSpaceHardeningVars)
{
    /* returns value in  hardening vector corresponding to isotropic hardening */
    if ( !isotropicHardeningFlag ) {
        return 0.;
    } else if ( this->kinematicHardeningFlag ) {
        return stressSpaceHardeningVars.at(7);
    } else {
        return stressSpaceHardeningVars.at(1);
    }

    return 0.;
}
} // end namespace oofem
