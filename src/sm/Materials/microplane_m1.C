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

#ifdef microplane_m1_new_implementation
// ========================= new implementation =========================
M1Material :: M1Material(int n, Domain *d) : MicroplaneMaterial(n, d)
{ E = 0.; nu = 0.; EN = 0.; s0 = 0.; HN = 0.; }

IRResultType
M1Material :: initializeFrom(InputRecord *ir)
{
    MicroplaneMaterial :: initializeFrom(ir);

    IRResultType result;                // Required by IR_GIVE_FIELD macro

    if ( nu != 0.25 ) {
        OOFEM_WARNING("Poisson ratio of microplane model M1 must be set to 0.25");
    }
    nu = 0.25; // read by MicroplaneMaterial, but overwritten here
    EN = E / ( 1. - 2. * nu );
    IR_GIVE_FIELD(ir, s0, _IFT_M1Material_s0);
    HN = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, HN, _IFT_M1Material_hn);
    ENtan = EN * HN / ( EN + HN );

    return IRRT_OK;
}

void
M1Material :: giveRealStressVector_3d(FloatArray &answer,
                                      GaussPoint *gp,
                                      const FloatArray &totalStrain,
                                      TimeStep *tStep)
{
    answer.resize(6);
    answer.zero();

    // get the status at the beginning
    M1MaterialStatus *status = static_cast< M1MaterialStatus * >( this->giveStatus(gp) );
    // prepare status at the end
    this->initTempStatus(gp);
    // get the initial values of plastic strains on microplanes (set to zero in the first step)
    FloatArray epspN = status->giveNormalMplanePlasticStrains();
    if ( epspN.giveSize() < numberOfMicroplanes ) {
        epspN.resize(numberOfMicroplanes);
        epspN.zero();
    }

    // loop over microplanes
    FloatArray sigN(numberOfMicroplanes);
    IntArray plState(numberOfMicroplanes);
    for ( int imp = 1; imp <= numberOfMicroplanes; imp++ ) {
        Microplane *mPlane = this->giveMicroplane(imp - 1, gp);
        //IntegrationPointStatus *mPlaneStatus =  this->giveMicroplaneStatus(mPlane);
        double epsN = computeNormalStrainComponent(mPlane, totalStrain);
        // evaluate trial stress on the microplane
        double sigTrial = EN * ( epsN - epspN.at(imp) );
        // evaluate the yield stress (from total microplane strain, not from its plastic part)
        double sigYield = EN * ( s0 + HN * epsN ) / ( EN + HN );
        if ( sigYield < 0. ) {
            sigYield = 0.;
        }
        // check whether the yield stress is exceeded and set the microplane stress
        if ( sigTrial > sigYield ) { //plastic yielding
            sigN.at(imp) = sigYield;
            epspN.at(imp) = epsN - sigYield / EN;
            plState.at(imp) = 1;
        } else {
            sigN.at(imp) = sigTrial;
            plState.at(imp) = 0;
        }
        // add the contribution of the microplane to macroscopic stresses
        for ( int i = 1; i <= 6; i++ ) {
            answer.at(i) += N [ imp - 1 ] [ i - 1 ] * sigN.at(imp) * microplaneWeights [ imp - 1 ];
        }
    }
    // multiply the integral over unit hemisphere by 6
    answer.times(6);

    // update status
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->letTempNormalMplaneStressesBe(sigN);
    status->letTempNormalMplanePlasticStrainsBe(epspN);
    status->letPlasticStateIndicatorsBe(plState);
}

void
M1Material :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
{
    answer.resize(6, 6);
    answer.zero();
    // elastic stiffness matrix
    if ( mode == ElasticStiffness ) {
        MicroplaneMaterial :: give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
        return;
    }

    M1MaterialStatus *status = static_cast< M1MaterialStatus * >( this->giveStatus(gp) );
    IntArray plasticState = status->givePlasticStateIndicators();
    if ( plasticState.giveSize() != numberOfMicroplanes ) {
        MicroplaneMaterial :: give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
        return;
    }
    // tangent stiffness matrix
    double aux, D11 = 0., D12 = 0., D13 = 0., D14 = 0., D15 = 0., D16 = 0., D22 = 0., D23 = 0., D24 = 0., D25 = 0., D26 = 0., D33 = 0., D34 = 0., D35 = 0., D36 = 0.;
    // loop over microplanes
    for ( int im = 0; im < numberOfMicroplanes; im++ ) {
        if ( plasticState.at(im + 1) ) {
            aux = ENtan * microplaneWeights [ im ];
        } else {
            aux = EN * microplaneWeights [ im ];
        }
        D11 += aux * N [ im ] [ 0 ] * N [ im ] [ 0 ];
        D12 += aux * N [ im ] [ 0 ] * N [ im ] [ 1 ];
        D13 += aux * N [ im ] [ 0 ] * N [ im ] [ 2 ];
        D14 += aux * N [ im ] [ 0 ] * N [ im ] [ 3 ];
        D15 += aux * N [ im ] [ 0 ] * N [ im ] [ 4 ];
        D16 += aux * N [ im ] [ 0 ] * N [ im ] [ 5 ];

        D22 += aux * N [ im ] [ 1 ] * N [ im ] [ 1 ];
        D23 += aux * N [ im ] [ 1 ] * N [ im ] [ 2 ];
        D24 += aux * N [ im ] [ 1 ] * N [ im ] [ 3 ];
        D25 += aux * N [ im ] [ 1 ] * N [ im ] [ 4 ];
        D26 += aux * N [ im ] [ 1 ] * N [ im ] [ 5 ];

        D33 += aux * N [ im ] [ 2 ] * N [ im ] [ 2 ];
        D34 += aux * N [ im ] [ 2 ] * N [ im ] [ 3 ];
        D35 += aux * N [ im ] [ 2 ] * N [ im ] [ 4 ];
        D36 += aux * N [ im ] [ 2 ] * N [ im ] [ 5 ];
    }
    answer.at(1, 1) = D11;
    answer.at(1, 2) = answer.at(2, 1) = answer.at(6, 6) = D12;
    answer.at(1, 3) = answer.at(3, 1) = answer.at(5, 5) = D13;
    answer.at(1, 4) = answer.at(4, 1) = answer.at(5, 6) = answer.at(6, 5) = D14;
    answer.at(1, 5) = answer.at(5, 1) = D15;
    answer.at(1, 6) = answer.at(6, 1) = D16;

    answer.at(2, 2) = D22;
    answer.at(2, 3) = answer.at(3, 2) = answer.at(4, 4) = D23;
    answer.at(2, 4) = answer.at(4, 2) = D24;
    answer.at(2, 5) = answer.at(5, 2) = answer.at(4, 6) = answer.at(6, 4) = D25;
    answer.at(2, 6) = answer.at(6, 2) = D26;

    answer.at(3, 3) = answer.at(3, 3) = D33;
    answer.at(3, 4) = answer.at(4, 3) = D34;
    answer.at(3, 5) = answer.at(5, 3) = D35;
    answer.at(3, 6) = answer.at(6, 3) = answer.at(4, 5)  = answer.at(5, 4) = D36;

    answer.times(6.);
}

int
M1Material :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    M1MaterialStatus *status = static_cast< M1MaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        // plastic strain is computed as total strain minus elastic strain
        // (note that integration of microplane plastic strains would give a different result)
        answer = status->giveStrainVector();
        FloatArray sig = status->giveStressVector();
        double aux = nu * ( sig.at(1) + sig.at(2) + sig.at(3) );
        double G = E / ( 2. * ( 1. + nu ) );
        answer.at(1) -= ( ( 1. + nu ) * sig.at(1) - aux ) / E;
        answer.at(2) -= ( ( 1. + nu ) * sig.at(2) - aux ) / E;
        answer.at(3) -= ( ( 1. + nu ) * sig.at(3) - aux ) / E;
        answer.at(4) -= sig.at(4) / G;
        answer.at(5) -= sig.at(5) / G;
        answer.at(6) -= sig.at(6) / G;
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}
////////////////////////////////////////////////////////////////////////////

M1MaterialStatus :: M1MaterialStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), sigN(0), epspN(0), plasticState(0)
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
        fprintf( file, " %g ", sigN.at(imp) );
    }
    fprintf(file, " epspN ");
    for ( int imp = 1; imp <= nm; imp++ ) {
        fprintf( file, " %g ", epspN.at(imp) );
    }
    fprintf(file, " plast ");
    for ( int imp = 1; imp <= nm; imp++ ) {
        fprintf( file, " %d ", plasticState.at(imp) );
    }
    fprintf(file, "}\n");
}

void
M1MaterialStatus :: updateYourself(TimeStep *tStep)
{
    sigN = tempSigN;
    epspN = tempEpspN;
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

#else
// ========================= old implementation =========================
M1Material :: M1Material(int n, Domain *d) : StructuralMaterial(n, d)
{ E = 0.; nu = 0.; EN = 0.; s0 = 0.; HN = 0.; }

void
M1Material :: giveRealStressVector_PlaneStress(FloatArray &answer,
                                               GaussPoint *gp,
                                               const FloatArray &totalStrain,
                                               TimeStep *tStep)
{
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

    for ( int imp = 1; imp <= nmp; imp++ ) {
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
        for ( int i = 1; i <= 3; i++ ) {
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
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

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

#endif // end of old implementation
} // end namespace oofem
