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

#include "microplane_m1.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(M1Material);

M1Material :: M1Material(int n, Domain *d) : StructuralMaterial(n, d)
{ E = 0.; nu = 0.; EN = 0.; s0 = 0.; HN = 0.; }

void
M1Material :: giveRealStressVector_PlaneStress(FloatArray &answer,
                                               GaussPoint *gp,
                                               const FloatArray &totalStrain,
                                               TimeStep *tStep)
{
    int i, imp;
    FloatArray sigmaN, deps, sigmaNyield;
    double depsN, epsN;

    answer.resize(3);
    answer.zero();
    sigmaNyield.resize(nmp);
    sigmaNyield.zero();

    M1MaterialStatus *status = static_cast< M1MaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    sigmaN = status->giveNormalMplaneStresses();
    if ( sigmaN.giveSize() < nmp ) {
        sigmaN.resize(nmp);
        sigmaN.zero();
    }
    deps.beDifferenceOf( totalStrain, status->giveStrainVector() );

    for ( imp = 1; imp <= nmp; imp++ ) {
        depsN = N.at(imp, 1) * deps.at(1) + N.at(imp, 2) * deps.at(2) + N.at(imp, 3) * deps.at(3);
        epsN = N.at(imp, 1) * totalStrain.at(1) + N.at(imp, 2) * totalStrain.at(2) + N.at(imp, 3) * totalStrain.at(3);
        sigmaN.at(imp) += EN * depsN;
        double sy = EN * ( s0 + HN * epsN ) / ( EN + HN ); // current microplane yield stress
        if ( sy < 0. ) {
            sy = 0.;
        }
        if ( sigmaN.at(imp) > sy ) {
            sigmaN.at(imp) = sy;
        }
        sigmaNyield.at(imp) = sy;
        for ( i = 1; i <= 3; i++ ) {
            answer.at(i) += N.at(imp, i) * sigmaN.at(imp) * mw.at(imp);
        }
    }

    // update status
    status->letTempStrainVectorBe(totalStrain);
    status->letTempNormalMplaneStressesBe(sigmaN);
    status->letNormalMplaneYieldStressesBe(sigmaNyield);
    status->letTempStressVectorBe(answer);
}

void
M1Material :: giveElasticPlaneStressStiffMtrx(FloatMatrix &answer)
{
    answer.resize(3, 3);
    answer.at(1, 1) = answer.at(2, 2) = E / ( 1. - nu * nu );
    answer.at(1, 2) = answer.at(2, 1) = E * nu / ( 1. - nu * nu );
    answer.at(1, 3) = answer.at(2, 3) = answer.at(3, 1) = answer.at(3, 2) = 0.;
    answer.at(3, 3) = E / ( 2. * ( 1. + nu ) );
}

void
M1Material :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(3, 3);
    if ( rMode == ElasticStiffness ) {
        giveElasticPlaneStressStiffMtrx(answer);
        return;
    }

    M1MaterialStatus *status = static_cast< M1MaterialStatus * >( this->giveStatus(gp) );
    FloatArray sigmaN = status->giveTempNormalMplaneStresses();
    if ( sigmaN.giveSize() != nmp ) {
        sigmaN = status->giveNormalMplaneStresses();
        if ( sigmaN.giveSize() != nmp ) {
            giveElasticPlaneStressStiffMtrx(answer);
            return;
        }
    }
    FloatArray sigmaNyield = status->giveNormalMplaneYieldStresses();

    double D11, D12, D13, D22, D23, aux;
    D11 = D12 = D13 = D22 = D23 = 0.;
    for ( int imp = 1; imp <= nmp; imp++ ) {
        if ( sigmaN.at(imp) < sigmaNyield.at(imp) ) { // otherwise the plane is yielding
            aux = mw.at(imp) * EN;
        } else if ( sigmaNyield.at(imp) > 0. ) {
            aux = mw.at(imp) * EN * HN / ( EN + HN );
        } else {
            aux = 0.;
        }
        D11 += aux * NN.at(imp, 1);
        D12 += aux * NN.at(imp, 3);
        D13 += aux * NN.at(imp, 2);
        D22 += aux * NN.at(imp, 5);
        D23 += aux * NN.at(imp, 4);
    }
    answer.at(1, 1) = D11;
    answer.at(1, 2) = answer.at(2, 1) = answer.at(3, 3) = D12;
    answer.at(1, 3) = answer.at(3, 1) = D13;
    answer.at(2, 2) = D22;
    answer.at(2, 3) = answer.at(3, 2) = D23;
}


IRResultType
M1Material :: initializeFrom(InputRecord *ir)
{
    StructuralMaterial :: initializeFrom(ir);

    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, E, _IFT_M1Material_e);
    EN = 1.5 * E;
    nu = 1. / 3.;
    IR_GIVE_FIELD(ir, s0, _IFT_M1Material_s0);
    HN = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, HN, _IFT_M1Material_hn);

    // specify microplanes
    IR_GIVE_FIELD(ir, nmp, _IFT_M1Material_nmp);
    n.resize(nmp, 2);
    N.resize(nmp, 3);
    NN.resize(nmp, 5);
    mw.resize(nmp);
    for ( int imp = 1; imp <= nmp; imp++ ) {
        double alpha = ( imp - 1 ) * ( M_PI / nmp );
        double c = cos(alpha);
        double s = sin(alpha);
        n.at(imp, 1) = c;
        n.at(imp, 2) = s;
        N.at(imp, 1) = c * c;
        N.at(imp, 2) = s * s;
        N.at(imp, 3) = c * s;
        NN.at(imp, 1) = c * c * c * c;
        NN.at(imp, 2) = c * c * c * s;
        NN.at(imp, 3) = c * c * s * s;
        NN.at(imp, 4) = c * s * s * s;
        NN.at(imp, 5) = s * s * s * s;
        mw.at(imp) = 2. / nmp;
    }

    return IRRT_OK;
}

int
M1Material :: hasMaterialModeCapability(MaterialMode mode)
{
    return mode == _PlaneStress;
}

int
M1Material :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    M1MaterialStatus *status = static_cast< M1MaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        FloatArray plasticStrain(3);
        plasticStrain.zero();
        FloatArray sigmaN = status->giveTempNormalMplaneStresses();
        FloatArray strain = status->giveTempStrainVector();
        double Exx, Eyy, Gamma;
        Exx = Eyy = Gamma = 0.;
        if ( sigmaN.giveSize() == nmp ) {
            for ( int imp = 1; imp <= nmp; imp++ ) {
                double epsN = 0.;
                for ( int i = 1; i <= 3; i++ ) {
                    epsN += strain.at(i) * N.at(imp, i);
                }
                double epsNpl = epsN - sigmaN.at(imp) / EN;
                double aux = epsNpl * mw.at(imp);
                Exx += aux * N.at(imp, 1);
                Eyy += aux * N.at(imp, 2);
                Gamma += aux * N.at(imp, 3);
            }
        }
        plasticStrain.at(1) = 1.5 * Exx - 0.5 * Eyy;
        plasticStrain.at(2) = -0.5 * Exx + 1.5 * Eyy;
        plasticStrain.at(3) = 4. * Gamma;
        StructuralMaterial :: giveFullSymVectorForm( answer, plasticStrain, gp->giveMaterialMode() );
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

////////////////////////////////////////////////////////////////////////////

M1MaterialStatus :: M1MaterialStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), sigN(0)
{}


M1MaterialStatus :: ~M1MaterialStatus()
{ }

void
M1MaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
}

void
M1MaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { sigN ");
    int nm = sigN.giveSize();
    for ( int imp = 1; imp <= nm; imp++ ) {
        fprintf( file, " %f ", sigN.at(imp) );
    }
    fprintf(file, "}\n");
}

void
M1MaterialStatus :: updateYourself(TimeStep *tStep)
{
    sigN = tempSigN;
    StructuralMaterialStatus :: updateYourself(tStep);
}

contextIOResultType
M1MaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    return StructuralMaterialStatus :: saveContext(stream, mode, obj);
}

contextIOResultType
M1MaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    return StructuralMaterialStatus :: restoreContext(stream, mode, obj);
}
} // end namespace oofem
