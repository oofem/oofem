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

#include "j2mplasticmaterial.h"
#include "Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(J2MPlasticMaterial);

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
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    double value;

    result = MPlasticMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;
    result = linearElasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, value, _IFT_J2MPlasticMaterial_ry);
    k = value / sqrt(3.0);

    kinematicModuli = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, kinematicModuli, _IFT_J2MPlasticMaterial_khm);

    isotropicModuli = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, isotropicModuli, _IFT_J2MPlasticMaterial_ihm);

    if ( fabs(kinematicModuli) > 1.e-12 ) {
        kinematicHardeningFlag = 1;
    }

    if ( fabs(isotropicModuli) > 1.e-12 ) {
        isotropicHardeningFlag = 1;
    }

    int rma = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, rma, _IFT_J2MPlasticMaterial_rma);
    if ( rma == 0 ) {
        this->rmType = mpm_ClosestPoint;
    } else {
        this->rmType = mpm_CuttingPlane;
    }

    return IRRT_OK;
}


MaterialStatus *
J2MPlasticMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new MPlasticMaterialStatus(1, this->giveDomain(), gp, this->giveSizeOfReducedHardeningVarsVector(gp));
}

void
J2MPlasticMaterial :: computeStressSpaceHardeningVars(FloatArray &answer, GaussPoint *gp,
                                                      const FloatArray &strainSpaceHardeningVariables)
{
    // in full stress strain space
    int count = 0, size = this->giveSizeOfFullHardeningVarsVector(), isize, rSize;
    IntArray mask;

    if ( !hasHardening() ) {
        answer.clear();
        return;
    }

    answer.resize(size);
    StructuralMaterial :: giveVoigtSymVectorMask( mask, gp->giveMaterialMode() );
    isize = mask.giveSize();
    rSize = this->giveSizeOfReducedHardeningVarsVector(gp);

    /* kinematic hardening variables are first */
    if ( this->kinematicHardeningFlag ) {
        for ( int i = 1; i <= isize; i++ ) {
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
                                                    TimeStep *tStep)
{
    /* computes hardening moduli in reduced stress strain space (for kinematic back-stress)*/

    int size = this->giveSizeOfReducedHardeningVarsVector(gp);

    if ( !hasHardening() ) {
        answer.clear();
        return;
    }

    answer.resize(size, size);
    answer.zero();

    /* kinematic hardening variables are first */
    if ( this->kinematicHardeningFlag ) {
        int ksize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
        for ( int i = 1; i <= ksize; i++ ) {
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

    int kcount = 0, size = this->giveSizeOfReducedHardeningVarsVector(gp);
    //double f,ax,ay,az,sx,sy,sz;
    FloatArray fullKinematicGradient, reducedKinematicGrad;

    if ( !hasHardening() ) {
        answer.clear();
        return;
    }

    answer.resize(size);

    /* kinematic hardening variables first */
    if ( this->kinematicHardeningFlag ) {
        this->computeStressGradientVector(fullKinematicGradient, ftype, isurf, gp, stressVector, stressSpaceHardeningVars);
        StructuralMaterial :: giveReducedSymVectorForm( reducedKinematicGrad, fullKinematicGradient, gp->giveMaterialMode() );

        kcount = reducedKinematicGrad.giveSize();
    }

    if ( this->kinematicHardeningFlag ) {
        for ( int i = 1; i <= kcount; i++ ) {
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
    int size;
    int imask, jmask;
    FloatArray helpVector, backStress, df(6);
    IntArray mask;
    double f, f32, f12, ax, ay, az;

    StructuralMaterial :: giveInvertedVoigtVectorMask( mask, gp->giveMaterialMode() );
    size = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) +
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
        f12 = sqrt(f);
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

        for ( int i = 1; i <= 3; i++ ) {
            if ( ( imask = mask.at(i) ) == 0 ) {
                continue;
            }

            for ( int j = i; j <= 3; j++ ) {
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

        for ( int i = 1; i <= 3; i++ ) {
            if ( ( imask = mask.at(i) ) == 0 ) {
                continue;
            }

            for ( int j = 4; j <= 6; j++ ) {
                if ( ( jmask = mask.at(j) ) == 0 ) {
                    continue;
                }

                answer.at(imask, jmask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j);
                answer.at(jmask, imask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j);
            }
        }

        for ( int i = 4; i <= 6; i++ ) {
            if ( ( imask = mask.at(i) ) == 0 ) {
                continue;
            }

            for ( int j = i; j <= 6; j++ ) {
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
                                             TimeStep *tStep)
{
    /* Returns 3d elastic moduli */
    this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, ElasticStiffness, gp, tStep);
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
J2MPlasticMaterial :: giveSizeOfReducedHardeningVarsVector(GaussPoint *gp) const
{
    /* Returns the size of hardening variables vector */
    int size = 0;

    if ( kinematicHardeningFlag ) {
        size += StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
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

    answer.clear();
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
}
} // end namespace oofem
