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

#include "microplanematerial.h"
#include "microplane.h"
#include "gausspnt.h"
#include "flotarry.h"
#include "contextioerr.h"
#include "mathfem.h"

namespace oofem {
Microplane *
MicroplaneMaterial :: giveMicroplane(int i, GaussPoint *masterGp)
{
    //
    // return the i-th slave integration  point (microplane) of master gp
    // if slaves of  gp don't exists - create them

    //
    GaussPoint *slave = masterGp->giveSlaveGaussPoint(i);
    if ( slave == NULL ) {
        // check for proper dimensions - slave can be NULL if index too high or if not
        // slaves previously defined
        if ( ( i < 0 ) || ( i > ( this->numberOfMicroplanes - 1 ) ) ) {
            _error("giveMicroplane: no such microplane defined");
        }

        // create new slave record in masterGp
        // (requires that this is friend of gp)

        MaterialMode slaveMode, masterMode = masterGp->giveMaterialMode();
        slaveMode = this->giveCorrespondingSlaveMaterialMode(masterMode);

        masterGp->numberOfGp = this->numberOfMicroplanes;
        masterGp->gaussPointArray = new GaussPoint * [ numberOfMicroplanes ];
        for ( int j = 0; j < numberOfMicroplanes; j++ ) {
            masterGp->gaussPointArray [ j ] = new Microplane(masterGp->giveIntegrationRule(), j + 1, slaveMode);
        }

        slave = masterGp->gaussPointArray [ i ];
    }

    return ( Microplane * ) slave;
}


MaterialMode
MicroplaneMaterial :: giveCorrespondingSlaveMaterialMode(MaterialMode masterMode)
{
    if ( masterMode == _3dMat ) {
        return _3dMicroplane;
    }

    return _Unknown;
}

MaterialStatus *
MicroplaneMaterial :: giveMicroplaneStatus(GaussPoint *gp)
/*
 * returns material status in gp corresponding to specific material class
 */
{
    MaterialStatus *status;

    status = gp->giveMaterialStatus();
    if ( status == NULL ) {
        // create a new one
        status = this->CreateMicroplaneStatus(gp);
        // if newly created status is null
        // dont include it. specific instance
        // does not have status.
        if ( status != NULL ) {
            gp->setMaterialStatus(status);
        }
    }

    return status;
}

void
MicroplaneMaterial :: initTempStatus(GaussPoint *gp)
{
    int mPlaneIndex;
    Microplane *mPlane;
    MaterialStatus *status;

    // init master
    this->giveStatus(gp)->initTempStatus();

    // init master microplanes
    for ( mPlaneIndex = 0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++ ) {
        mPlane = this->giveMicroplane(mPlaneIndex, gp);
        status = this->giveMicroplaneStatus(mPlane);
        if ( status ) {
            status->initTempStatus();
        }
    }
}

contextIOResultType
MicroplaneMaterial :: saveIPContext(DataStream *stream, ContextMode mode, GaussPoint* gp)
{
    contextIOResultType iores;
    int mPlaneIndex;
    Microplane *mPlane;
    MaterialStatus *status;

    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( gp == NULL ) {
        THROW_CIOERR(CIO_BADOBJ);
    }

    // save master
    if ( ( iores = StructuralMaterial :: saveIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // save microplanes
    for ( mPlaneIndex = 0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++ ) {
        mPlane = this->giveMicroplane(mPlaneIndex, gp);
        status =  this->giveMicroplaneStatus(mPlane);
        if ( status ) {
            if ( ( iores = status->saveContext(stream, mode, gp) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        }
    }

    return CIO_OK;
}

contextIOResultType
MicroplaneMaterial :: restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp)
{
    contextIOResultType iores;
    int mPlaneIndex;
    Microplane *mPlane;
    MaterialStatus *status;

    // corresponding gp is passed in obj
    if ( gp == NULL ) {
        THROW_CIOERR(CIO_BADOBJ);
    }

    // save master
    if ( ( iores = StructuralMaterial :: restoreIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // save microplanes
    for ( mPlaneIndex = 0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++ ) {
        mPlane = this->giveMicroplane(mPlaneIndex, gp);
        status =  this->giveMicroplaneStatus(mPlane);
        if ( status ) {
            if ( ( iores = status->restoreContext(stream, mode, gp) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        }
    }

    return CIO_OK;
}




void
MicroplaneMaterial ::  giveMicroplaneNormal(FloatArray &answer, Microplane *mplane)
{
    int mnumber = mplane->giveNumber();
    answer.resize(3);

    answer.at(1) = microplaneNormals [ mnumber - 1 ] [ 0 ];
    answer.at(2) = microplaneNormals [ mnumber - 1 ] [ 1 ];
    answer.at(3) = microplaneNormals [ mnumber - 1 ] [ 2 ];
}

double
MicroplaneMaterial :: giveMicroplaneIntegrationWeight(Microplane *mplane)
{
    int mnumber = mplane->giveNumber();

    return microplaneWeights [ mnumber - 1 ];
}

double
MicroplaneMaterial :: computeNormalStrainComponent(Microplane *mplane,
                                                   const FloatArray &macroStrain)
{
    int mnumber = mplane->giveNumber();
    int i;
    double en = 0.0;

    for ( i = 0; i < 6; i++ ) {
        en += N [ mnumber - 1 ] [ i ] * macroStrain.at(i + 1);
    }

    return en;
}

double
MicroplaneMaterial :: computeNormalVolumetricStrainComponent(Microplane *mplane,
                                                             const FloatArray &macroStrain)
{
    return ( macroStrain.at(1) + macroStrain.at(2) + macroStrain.at(3) ) / 3.0;
}

double
MicroplaneMaterial :: computeNormalDeviatoricStrainComponent(Microplane *mplane,
                                                             const FloatArray &macroStrain)
{
    return this->computeNormalStrainComponent(mplane, macroStrain) -
           this->computeNormalVolumetricStrainComponent(mplane, macroStrain);
}

double
MicroplaneMaterial :: computeShearMStrainComponent(Microplane *mplane,
                                                   const FloatArray &macroStrain)
{
    int mnumber = mplane->giveNumber();
    int i;
    double em = 0.0;

    for ( i = 0; i < 6; i++ ) {
        em += M [ mnumber - 1 ] [ i ] * macroStrain.at(i + 1);
    }

    return em;
}

double
MicroplaneMaterial :: computeShearLStrainComponent(Microplane *mplane,
                                                   const FloatArray &macroStrain)
{
    int mnumber = mplane->giveNumber();
    int i;
    double el = 0.0;

    for ( i = 0; i < 6; i++ ) {
        el += L [ mnumber - 1 ] [ i ] * macroStrain.at(i + 1);
    }

    return el;
}

void
MicroplaneMaterial :: computeStrainVectorComponents(FloatArray &answer, Microplane *mplane,
                                                    const FloatArray &macroStrain)
{
    int mnumber = mplane->giveNumber();
    double en = 0., ev = 0., em = 0., el = 0.;
    int i;

    for ( i = 0; i < 6; i++ ) {
        en += N [ mnumber - 1 ] [ i ] * macroStrain.at(i + 1);
        em += M [ mnumber - 1 ] [ i ] * macroStrain.at(i + 1);
        el += L [ mnumber - 1 ] [ i ] * macroStrain.at(i + 1);
    }

    ev = ( macroStrain.at(1) + macroStrain.at(2) + macroStrain.at(3) ) / 3.0;

    answer.resize(4);
    answer.zero();
    answer.at(1) = ev;
    answer.at(2) = en;
    answer.at(3) = el;
    answer.at(4) = em;
}

IRResultType
MicroplaneMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StructuralMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, numberOfMicroplanes, IFT_MicroplaneMaterial_nmp, "nmp"); // Macro
    this->initializeData(numberOfMicroplanes);
    return IRRT_OK;
}


int
MicroplaneMaterial :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    StructuralMaterial :: giveInputRecordString(str, keyword);
    sprintf(buff, " nmp %d", this->numberOfMicroplanes);
    str += buff;

    return 1;
}


void
MicroplaneMaterial :: initializeData(int numberOfMicroplanes)
{
    int i, i4, mPlane;
    int ij [ 6 ] [ 2 ] = { { 1, 1 }, { 2, 2 }, { 3, 3 }, { 2, 3 }, { 3, 1 }, { 1, 2 } };

    // kronecker delta
    Kronecker [ 0 ] = 1.;
    Kronecker [ 1 ] = 1.;
    Kronecker [ 2 ] = 1.;
    Kronecker [ 3 ] = 0.;
    Kronecker [ 4 ] = 0.;
    Kronecker [ 5 ] = 0.;

    // definition of intergation points and weights

    if ( numberOfMicroplanes == 28 ) {
        for ( i = 0; i < 4; i++ ) {
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
        for ( i = 0; i < 7; i++ ) {
            i4 = i * 4;
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
         * for (i=0; i<6; i++) {
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
        for ( i = 0; i < 3; i++ ) {
            microplaneWeights [ i ]    = 0.02652141274;
        }

        for ( i = 3; i < 9; i++ ) {
            microplaneWeights [ i ]    = 0.01993014153;
        }

        for ( i = 9; i < 21; i++ ) {
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

        for ( i = 0; i < 3; i++ ) {
            i4 = i * 4;
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
         * } else  {
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

        for ( i = 0; i < numberOfMicroplanes; i++ ) {
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
         * for (mPlane=0; mPlane < numberOfMicroplanes; mPlane++) {
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
        _error("initializeData: Unsupported number of microplanes");
    }

    // compute projection tensors for each microplane
    int ii, jj;
    FloatArray n(3), m(3), l(3);
    double aux;

    for ( mPlane = 0; mPlane < numberOfMicroplanes; mPlane++ ) {
        n.at(1) = microplaneNormals [ mPlane ] [ 0 ];
        n.at(2) = microplaneNormals [ mPlane ] [ 1 ];
        n.at(3) = microplaneNormals [ mPlane ] [ 2 ];

        if ( ( mPlane + 1 ) % 3 == 0 ) {
            aux = sqrt( n.at(1) * n.at(1) + n.at(2) * n.at(2) );
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
            aux = sqrt( n.at(2) * n.at(2) + n.at(3) * n.at(3) );
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
            aux = sqrt( n.at(1) * n.at(1) + n.at(3) * n.at(3) );
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

        for ( i = 0; i < 6; i++ ) {
            ii = ij [ i ] [ 0 ];
            jj = ij [ i ] [ 1 ];

            N [ mPlane ] [ i ] = n.at(ii) * n.at(jj);
            M [ mPlane ] [ i ] = 0.5 * ( m.at(ii) * n.at(jj) + m.at(jj) * n.at(ii) );
            L [ mPlane ] [ i ] = 0.5 * ( l.at(ii) * n.at(jj) + l.at(jj) * n.at(ii) );
        }
    }
}
} // end namespace oofem
