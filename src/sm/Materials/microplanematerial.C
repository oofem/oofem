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

#include "microplanematerial.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"

namespace oofem {

double
MicroplaneMaterial :: giveMicroplaneIntegrationWeight(int mnumber) const
{
    return microplaneWeights [ mnumber - 1 ];
}

double
MicroplaneMaterial :: computeNormalStrainComponent(int mnumber,
                                                   const FloatArray &macroStrain) const
{
    double en = 0.0;
    for ( int i = 0; i < 6; i++ ) {
        en += N [ mnumber - 1 ] [ i ] * macroStrain[ i ];
    }

    return en;
}

double
MicroplaneMaterial :: computeNormalVolumetricStrainComponent(const FloatArray &macroStrain) const
{
    return ( macroStrain.at(1) + macroStrain.at(2) + macroStrain.at(3) ) / 3.0;
}

double
MicroplaneMaterial :: computeNormalDeviatoricStrainComponent(int mnumber,
                                                             const FloatArray &macroStrain) const
{
    return this->computeNormalStrainComponent(mnumber, macroStrain) -
           this->computeNormalVolumetricStrainComponent(macroStrain);
}

double
MicroplaneMaterial :: computeShearMStrainComponent(int mnumber,
                                                   const FloatArray &macroStrain) const
{
    double em = 0.0;
    for ( int i = 0; i < 6; i++ ) {
        em += M [ mnumber - 1 ] [ i ] * macroStrain[ i ];
    }

    return em;
}

double
MicroplaneMaterial :: computeShearLStrainComponent(int mnumber,
                                                   const FloatArray &macroStrain) const
{
    double el = 0.0;

    for ( int i = 0; i < 6; i++ ) {
        el += L [ mnumber - 1 ] [ i ] * macroStrain[ i ];
    }

    return el;
}

MicroplaneState
MicroplaneMaterial :: computeStrainVectorComponents(int mnumber,
                                                    const FloatArray &macroStrain) const
{
    MicroplaneState e;

    for ( int i = 0; i < 6; i++ ) {
        e.n += N [ mnumber - 1 ] [ i ] * macroStrain [ i ];
        e.m += M [ mnumber - 1 ] [ i ] * macroStrain [ i ];
        e.l += L [ mnumber - 1 ] [ i ] * macroStrain [ i ];
    }

    e.v = ( macroStrain.at(1) + macroStrain.at(2) + macroStrain.at(3) ) / 3.0;
    return e;
}


FloatMatrixF<6,6>
MicroplaneMaterial :: give3dMaterialStiffnessMatrix(MatResponseMode mode,
                                                    GaussPoint *gp,
                                                    TimeStep *tStep) const
{
    FloatMatrixF<6,6> answer;
    // elastic stiffness matrix
    answer.at(4, 4) = answer.at(5, 5) = answer.at(6, 6) = E / ( 2. + 2. * nu );
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = E * ( 1. - nu ) / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
    answer.at(1, 2) = answer.at(2, 1) = answer.at(1, 3) = answer.at(3, 1) =
            answer.at(2, 3) = answer.at(3, 2) = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
    return answer;
}

void
MicroplaneMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    // elastic constants
    IR_GIVE_FIELD(ir, E, _IFT_MicroplaneMaterial_e);
    IR_GIVE_FIELD(ir, nu, _IFT_MicroplaneMaterial_n);
    // number of microplanes
    IR_GIVE_FIELD(ir, numberOfMicroplanes, _IFT_MicroplaneMaterial_nmp);
    this->initializeData(numberOfMicroplanes);
}


void
MicroplaneMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
    input.setField(this->numberOfMicroplanes, _IFT_MicroplaneMaterial_nmp);
}


void
MicroplaneMaterial :: initializeData(int numberOfMicroplanes)
{
    int ij [ 6 ] [ 2 ] = { { 1, 1 }, { 2, 2 }, { 3, 3 }, { 2, 3 }, { 3, 1 }, { 1, 2 } };

    // kronecker delta
    Kronecker [ 0 ] = 1.;
    Kronecker [ 1 ] = 1.;
    Kronecker [ 2 ] = 1.;
    Kronecker [ 3 ] = 0.;
    Kronecker [ 4 ] = 0.;
    Kronecker [ 5 ] = 0.;

    // definition of intergation points and weights
    microplaneNormals.resize(numberOfMicroplanes);
    microplaneWeights.resize(numberOfMicroplanes);
    M.resize(numberOfMicroplanes);
    N.resize(numberOfMicroplanes);
    L.resize(numberOfMicroplanes);

    if ( numberOfMicroplanes == 28 ) {
        for ( int i = 0; i < 4; i++ ) {
            microplaneWeights [ i ]    = 0.0160714276E0;
            microplaneWeights [ i + 4 ]  = 0.0204744730E0;
            microplaneWeights [ i + 8 ]  = 0.0204744730E0;
            microplaneWeights [ i + 12 ] = 0.0204744730E0;
            microplaneWeights [ i + 16 ] = 0.0158350505E0;
            microplaneWeights [ i + 20 ] = 0.0158350505E0;
            microplaneWeights [ i + 24 ] = 0.0158350505E0;
        }

        double help [ 7 ] [ 3 ] = {
            { .577350259E0, .577350259E0, .577350259E0 }, { .935113132E0, .250562787E0, .250562787E0 },
            { .250562787E0, .935113132E0, .250562787E0 }, { .250562787E0, .250562787E0, .935113132E0 },
            { .186156720E0, .694746614E0, .694746614E0 }, { .694746614E0, .186156720E0, .694746614E0 },
            { .694746614E0, .694746614E0, .186156720E0 }
        };
        for ( int i = 0; i < 7; i++ ) {
            int i4 = i * 4;
            microplaneNormals [ i4  ] [ 0 ] = help [ i ] [ 0 ];
            microplaneNormals [ i4  ] [ 1 ] = help [ i ] [ 1 ];
            microplaneNormals [ i4  ] [ 2 ] = help [ i ] [ 2 ];
            microplaneNormals [ i4 + 1 ] [ 0 ] = help [ i ] [ 0 ];
            microplaneNormals [ i4 + 1 ] [ 1 ] = help [ i ] [ 1 ];
            microplaneNormals [ i4 + 1 ] [ 2 ] = -help [ i ] [ 2 ];
            microplaneNormals [ i4 + 2 ] [ 0 ] = help [ i ] [ 0 ];
            microplaneNormals [ i4 + 2 ] [ 1 ] = -help [ i ] [ 1 ];
            microplaneNormals [ i4 + 2 ] [ 2 ] = help [ i ] [ 2 ];
            microplaneNormals [ i4 + 3 ] [ 0 ] = help [ i ] [ 0 ];
            microplaneNormals [ i4 + 3 ] [ 1 ] = -help [ i ] [ 1 ];
            microplaneNormals [ i4 + 3 ] [ 2 ] = -help [ i ] [ 2 ];
        }

        /*
         * // compute projection tensors for each microplane
         * int ii,jj;
         * FloatArray n(3), m(3), l(3);
         * double aux;
         *
         * for (mPlane=0; mPlane < numberOfMicroplanes; mPlane++) {
         * //prohozemo poradi os normal, aby smer 3 byl na prvnim miste
         * n.at(1) = microplaneNormals[mPlane][2];
         * n.at(2) = microplaneNormals[mPlane][1];
         * n.at(3) = microplaneNormals[mPlane][0];
         *
         * if ((mPlane+1)%3 == 0) {
         * aux = sqrt (n.at(1)*n.at(1) + n.at(2)*n.at(2));
         * m.at(1) = n.at(2)/aux;
         * m.at(2) = -n.at(1)/aux;
         * m.at(3) = 0.0;
         * } else if ((mPlane+1)%3 == 1) {
         * aux = sqrt (n.at(2)*n.at(2) + n.at(3)*n.at(3));
         * m.at(1) = 0.0;
         * m.at(2) = n.at(3)/aux;
         * m.at(3) = -n.at(2)/aux;
         * } else  {
         * aux = sqrt (n.at(1)*n.at(1) + n.at(3)*n.at(3));
         * m.at(1) = n.at(3)/aux;
         * m.at(2) = 0.0;
         * m.at(3) = -n.at(1)/aux;
         * }
         *
         *
         * l.beVectorProductOf (m,n);
         *
         * for (int i=0; i<6; i++) {
         * ii = ij[i][0];
         * jj = ij[i][1];
         * N[mPlane][i] = n.at(ii)*n.at(jj);
         * M[mPlane][i] = 0.5*(m.at(ii)*n.at(jj) + m.at(jj)*n.at(ii)) ;
         * L[mPlane][i] = 0.5*(l.at(ii)*n.at(jj) + l.at(jj)*n.at(ii)) ;
         * }
         *
         *
         * } */
    } else if ( numberOfMicroplanes == 21 ) {
        for ( int i = 0; i < 3; i++ ) {
            microplaneWeights [ i ]    = 0.02652141274;
        }

        for ( int i = 3; i < 9; i++ ) {
            microplaneWeights [ i ]    = 0.01993014153;
        }

        for ( int i = 9; i < 21; i++ ) {
            microplaneWeights [ i ]    = 0.02507124272;
        }


        microplaneNormals [ 0 ] [ 0 ] =  microplaneNormals [ 1 ] [ 1 ] =  microplaneNormals [ 2 ] [ 2 ] = 1.;
        microplaneNormals [ 0 ] [ 1 ] =  microplaneNormals [ 0 ] [ 2 ] =  microplaneNormals [ 1 ] [ 0 ] =
                                                                             microplaneNormals [ 1 ] [ 2 ] = microplaneNormals [ 2 ] [ 0 ] = microplaneNormals [ 2 ] [ 1 ] = 0.;
        microplaneNormals [ 3 ] [ 0 ] =  microplaneNormals [ 3 ] [ 1 ] =  microplaneNormals [ 4 ] [ 0 ] =
                                                                             microplaneNormals [ 5 ] [ 0 ] = microplaneNormals [ 5 ] [ 2 ] = microplaneNormals [ 6 ] [ 0 ] =
                                                                                                                                                 microplaneNormals [ 7 ] [ 1 ] = microplaneNormals [ 7 ] [ 2 ] = microplaneNormals [ 8 ] [ 1 ] =
                                                                                                                                                                                                                     0.7071067812;
        microplaneNormals [ 4 ] [ 1 ] =  microplaneNormals [ 6 ] [ 2 ] =  microplaneNormals [ 8 ] [ 2 ] =
                                                                             -0.7071067812;
        microplaneNormals [ 3 ] [ 2 ] =  microplaneNormals [ 4 ] [ 2 ] =  microplaneNormals [ 5 ] [ 1 ] =
                                                                             microplaneNormals [ 6 ] [ 1 ] =  microplaneNormals [ 7 ] [ 0 ] =  microplaneNormals [ 8 ] [ 0 ] = 0.;


        double help [ 3 ] [ 3 ] = { { 0.3879072746, 0.3879072746, 0.8360956240 },
                                    { 0.3879072746, 0.8360956240, 0.3879072746 },
                                    { 0.8360956240, 0.3879072746, 0.3879072746 } };

        for ( int i = 0; i < 3; i++ ) {
            int i4 = i * 4;
            microplaneNormals [ 9 + i4 ] [ 0 ] = help [ i ] [ 0 ];
            microplaneNormals [ 9 + i4 ] [ 1 ] = help [ i ] [ 1 ];
            microplaneNormals [ 9 + i4 ] [ 2 ] = help [ i ] [ 2 ];

            microplaneNormals [ 10 + i4 ] [ 0 ] = help [ i ] [ 0 ];
            microplaneNormals [ 10 + i4 ] [ 1 ] = help [ i ] [ 1 ];
            microplaneNormals [ 10 + i4 ] [ 2 ] = -help [ i ] [ 2 ];

            microplaneNormals [ 11 + i4 ] [ 0 ] = help [ i ] [ 0 ];
            microplaneNormals [ 11 + i4 ] [ 1 ] = -help [ i ] [ 1 ];
            microplaneNormals [ 11 + i4 ] [ 2 ] = help [ i ] [ 2 ];

            microplaneNormals [ 12 + i4 ] [ 0 ] = help [ i ] [ 0 ];
            microplaneNormals [ 12 + i4 ] [ 1 ] = -help [ i ] [ 1 ];
            microplaneNormals [ 12 + i4 ] [ 2 ] = -help [ i ] [ 2 ];
        }


        /*  // compute projection tensors for each microplane
         * int ii,jj;
         * FloatArray n(3), m(3), l(3);
         * double aux;
         *
         * for (mPlane=0; mPlane < numberOfMicroplanes; mPlane++) {
         * n.at(1) = microplaneNormals[mPlane][0];
         * n.at(2) = microplaneNormals[mPlane][1];
         * n.at(3) = microplaneNormals[mPlane][2];
         *
         * if ((mPlane+1)%3 == 0) {
         *  aux = sqrt (n.at(1)*n.at(1) + n.at(2)*n.at(2));
         *  if(fabs(aux) > 1.0e-8){
         *  m.at(1) = n.at(2)/aux;
         *  m.at(2) = -n.at(1)/aux;
         *  m.at(3) = 0.0;
         *  }
         *  else {
         *    m.at(1) = 1.0;
         *    m.at(2) = 0.0;
         *    m.at(3) = 0.0;
         *  }
         * } else if ((mPlane+1)%3 == 1) {
         *  aux = sqrt (n.at(2)*n.at(2) + n.at(3)*n.at(3));
         *  if(fabs(aux) > 1.0e-8){
         *  m.at(1) = 0.0;
         *  m.at(2) = n.at(3)/aux;
         *  m.at(3) = -n.at(2)/aux;
         *  }
         *  else {
         *    m.at(1) = 0.0;
         *    m.at(2) = 1.0;
         *    m.at(3) = 0.0;
         *  }
         * } else {
         *  aux = sqrt (n.at(1)*n.at(1) + n.at(3)*n.at(3));
         *  if(fabs(aux) > 1.0e-8){
         *    m.at(1) = n.at(3)/aux;
         *    m.at(2) = 0.0;
         *    m.at(3) = -n.at(1)/aux;
         *  }
         *  else {
         *    m.at(1) = 0.0;
         *    m.at(2) = 0.0;
         *    m.at(3) = 1.0;
         *  }
         * }
         *
         *
         * l.beVectorProductOf (m,n);
         *
         * for (i=0; i<6; i++) {
         *  ii = ij[i][0];
         *  jj = ij[i][1];
         *
         *  N[mPlane][i] = n.at(ii)*n.at(jj);
         *  M[mPlane][i] = 0.5*(m.at(ii)*n.at(jj) + m.at(jj)*n.at(ii)) ;
         *  L[mPlane][i] = 0.5*(l.at(ii)*n.at(jj) + l.at(jj)*n.at(ii)) ;
         *
         *  // printf("\n i= %d,N= %lf,L= %lf,M= %lf \n",i,N[mPlane][i],L[mPlane][i],M[mPlane][i]);
         * }
         *
         * }
         *
         */
    } else if ( numberOfMicroplanes == 61 ) {
        double help [ 61 ] [ 4 ] = {
            { 1.000000000000, 0.000000000000, 0.000000000000, 0.00795844204678 },
            { 0.745355992500, 0.0, 0.666666666667, 0.00795844204678 },
            { 0.745355992500, -0.577350269190, -0.333333333333, 0.00795844204678 },
            { 0.745355992500, 0.577350269190, -0.333333333333, 0.00795844204678 },
            { 0.333333333333, 0.577350269190, 0.745355992500, 0.00795844204678 },
            { 0.333333333333, -0.577350269190, 0.745355992500, 0.00795844204678 },
            { 0.333333333333, -0.934172358963, 0.127322003750, 0.00795844204678 },
            { 0.333333333333, -0.356822089773, -0.872677996250, 0.00795844204678 },
            { 0.333333333333, 0.356822089773, -0.872677996250, 0.00795844204678 },
            { 0.333333333333, 0.934172358963, 0.127322003750, 0.00795844204678 },
            { 0.794654472292, -0.525731112119, 0.303530999103, 0.0105155242892 },
            { 0.794654472292, 0.0, -0.607061998207, 0.0105155242892 },
            { 0.794654472292, 0.525731112119, 0.303530999103, 0.0105155242892 },
            { 0.187592474085, 0.0, 0.982246946377, 0.0105155242892 },
            { 0.187592474085, -0.850650808352, -0.491123473188, 0.0105155242892 },
            { 0.187592474085, 0.850650808352, -0.491123473188, 0.0105155242892 },
            { 0.934172358963, 0.0, 0.356822089773, 0.0100119364272 },
            { 0.934172358963, -0.309016994375, -0.178411044887, 0.0100119364272 },
            { 0.934172358963, 0.309016994375, -0.178411044887, 0.0100119364272 },
            { 0.577350269190, 0.309016994375, 0.755761314076, 0.0100119364272 },
            { 0.577350269190, -0.309016994375, 0.755761314076, 0.0100119364272 },
            { 0.577350269190, -0.809016994375, -0.110264089708, 0.0100119364272 },
            { 0.577350269190, -0.5, -0.645497224368, 0.0100119364272 },
            { 0.577350269190, 0.5, -0.645497224368, 0.0100119364272 },
            { 0.577350269190, 0.809016994375, -0.110264089708, 0.0100119364272 },
            { 0.356822089773, -0.809016994375, 0.467086179481, 0.0100119364272 },
            { 0.356822089773, 0.0, -0.934172358963, 0.0100119364272 },
            { 0.356822089773, 0.809016994375, 0.467086179481, 0.0100119364272 },
            { 0.0, 0.5, 0.866025403784, 0.0100119364272 },
            { 0.0, -1.0, 0.0, 0.0100119364272 },
            { 0.0, 0.5, -0.866025403784, 0.0100119364272 },
            { 0.947273580412, -0.277496978165, 0.160212955043, 0.00690477957966 },
            { 0.812864676392, -0.277496978165, 0.512100034157, 0.00690477957966 },
            { 0.595386501297, -0.582240127941, 0.553634669695, 0.00690477957966 },
            { 0.595386501297, -0.770581752342, 0.227417407053, 0.00690477957966 },
            { 0.812864676392, -0.582240127941, -0.015730584514, 0.00690477957966 },
            { 0.492438766306, -0.753742692223, -0.435173546254, 0.00690477957966 },
            { 0.274960591212, -0.942084316623, -0.192025554687, 0.00690477957966 },
            { -0.076926487903, -0.942084316623, -0.326434458707, 0.00690477957966 },
            { -0.076926487903, -0.753742692223, -0.652651721349, 0.00690477957966 },
            { 0.274960591212, -0.637341166847, -0.719856173359, 0.00690477957966 },
            { 0.947273580412, 0.0, -0.320425910085, 0.00690477957966 },
            { 0.812864676392, -0.304743149777, -0.496369449643, 0.00690477957966 },
            { 0.595386501297, -0.188341624401, -0.781052076747, 0.00690477957966 },
            { 0.595386501297, 0.188341624401, -0.781052076747, 0.00690477957966 },
            { 0.812864676392, 0.304743149777, -0.496369449643, 0.00690477957966 },
            { 0.492438766306, 0.753742692223, -0.435173546254, 0.00690477957966 },
            { 0.274960591212, 0.637341166847, -0.719856173359, 0.00690477957966 },
            { -0.076926487903, 0.753742692223, -0.652651721349, 0.00690477957966 },
            { -0.076926487903, 0.942084316623, -0.326434458707, 0.00690477957966 },
            { 0.274960591212, 0.942084316623, -0.192025554687, 0.00690477957966 },
            { 0.947273580412, 0.277496978165, 0.160212955043, 0.00690477957966 },
            { 0.812864676392, 0.582240127941, -0.015730584514, 0.00690477957966 },
            { 0.595386501297, 0.770581752342, 0.227417407053, 0.00690477957966 },
            { 0.595386501297, 0.582240127941, 0.553634669695, 0.00690477957966 },
            { 0.812864676392, 0.277496978165, 0.512100034157, 0.00690477957966 },
            { 0.492438766306, 0.0, 0.870347092509, 0.00690477957966 },
            { 0.274960591212, 0.304743149777, 0.911881728046, 0.00690477957966 },
            { -0.076926487903, 0.188341624401, 0.979086180056, 0.00690477957966 },
            { -0.076926487903, -0.188341624401, 0.979086180056, 0.00690477957966 },
            { 0.274960591212, -0.304743149777, 0.911881728046, 0.00690477957966 }
        };

        for ( int i = 0; i < numberOfMicroplanes; i++ ) {
            microplaneNormals [ i ] [ 0 ] = help [ i ] [ 0 ];
            microplaneNormals [ i ] [ 1 ] = help [ i ] [ 1 ];
            microplaneNormals [ i ] [ 2 ] = help [ i ] [ 2 ];
            microplaneWeights [ i ]   = help [ i ] [ 3 ];
        }

        /*
         * // compute projection tensors for each microplane
         * int ii,jj;
         * FloatArray n(3), m(3), l(3);
         * double aux;
         *
         * for (int mPlane=0; mPlane < numberOfMicroplanes; mPlane++) {
         * n.at(1) = microplaneNormals[mPlane][0];
         * n.at(2) = microplaneNormals[mPlane][1];
         * n.at(3) = microplaneNormals[mPlane][2];
         *
         * //if ((mPlane+1)%3 == 0) {
         *  aux = sqrt (n.at(1)*n.at(1) + n.at(2)*n.at(2));
         *  if(fabs(aux) > 1.0e-8){
         * m.at(1) = n.at(2)/aux;
         * m.at(2) = -n.at(1)/aux;
         * m.at(3) = 0.0;
         *  }
         *  else {
         * m.at(1) = 1.0;
         * m.at(2) = 0.0;
         * m.at(3) = 0.0;
         *  }
         *
         * l.beVectorProductOf (m,n);
         *
         *
         *
         * for (i=0; i<6; i++) {
         *  ii = ij[i][0];
         *  jj = ij[i][1];
         *
         *  N[mPlane][i] = n.at(ii)*n.at(jj);
         *  M[mPlane][i] = 0.5*(m.at(ii)*n.at(jj) + m.at(jj)*n.at(ii)) ;
         *  L[mPlane][i] = 0.5*(l.at(ii)*n.at(jj) + l.at(jj)*n.at(ii)) ;
         *
         *
         * }
         *
         *
         *
         * }
         */
    } else {
        OOFEM_ERROR("Unsupported number of microplanes");
    }

    // compute projection tensors for each microplane
    FloatArray m(3), l(3);

    for ( int mPlane = 0; mPlane < numberOfMicroplanes; mPlane++ ) {
        auto n = microplaneNormals [ mPlane ];

        if ( ( mPlane + 1 ) % 3 == 0 ) {
            double aux = sqrt( n.at(1) * n.at(1) + n.at(2) * n.at(2) );
            if ( fabs(aux) > 1.0e-12 ) {
                m.at(1) = n.at(2) / aux;
                m.at(2) = -n.at(1) / aux;
                m.at(3) = 0.0;
            } else {
                m.at(1) = 1.0;
                m.at(2) = 0.0;
                m.at(3) = 0.0;
            }
        } else if ( ( mPlane + 1 ) % 3 == 1 ) {
            double aux = sqrt( n.at(2) * n.at(2) + n.at(3) * n.at(3) );
            if ( fabs(aux) > 1.0e-12 ) {
                m.at(1) = 0.0;
                m.at(2) = n.at(3) / aux;
                m.at(3) = -n.at(2) / aux;
            } else {
                m.at(1) = 0.0;
                m.at(2) = 1.0;
                m.at(3) = 0.0;
            }
        } else {
            double aux = sqrt( n.at(1) * n.at(1) + n.at(3) * n.at(3) );
            if ( fabs(aux) > 1.0e-12 ) {
                m.at(1) = n.at(3) / aux;
                m.at(2) = 0.0;
                m.at(3) = -n.at(1) / aux;
            } else {
                m.at(1) = 0.0;
                m.at(2) = 0.0;
                m.at(3) = 1.0;
            }
        }


        l.beVectorProductOf(m, n);

        for ( int i = 0; i < 6; i++ ) {
            int ii = ij [ i ] [ 0 ];
            int jj = ij [ i ] [ 1 ];

            N [ mPlane ] [ i ] = n.at(ii) * n.at(jj);
            M [ mPlane ] [ i ] = 0.5 * ( m.at(ii) * n.at(jj) + m.at(jj) * n.at(ii) );
            L [ mPlane ] [ i ] = 0.5 * ( l.at(ii) * n.at(jj) + l.at(jj) * n.at(ii) );
        }
    }
}
} // end namespace oofem
