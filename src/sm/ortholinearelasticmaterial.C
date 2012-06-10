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

#include "linearelasticmaterial.h"
#include "ortholinearelasticmaterial.h"
#include "structuralelement.h"
#include "material.h"
#include "structuralms.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#include "mathfem.h"

namespace oofem {
#define ZERO_LENGTH 1.e-6

IRResultType
OrthotropicLinearElasticMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    double value;
    int j, size;
    FloatArray triplets;


    this->LinearElasticMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_ex, "ex"); // Macro
    propertyDictionary->add(Ex, value);

    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_ey, "ey"); // Macro
    propertyDictionary->add(Ey, value);

    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_ez, "ez"); // Macro
    propertyDictionary->add(Ez, value);


    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_nyyz, "nyyz"); // Macro
    propertyDictionary->add(NYyz, value);

    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_nyxz, "nyxz"); // Macro
    propertyDictionary->add(NYxz, value);

    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_nyxy, "nyxy"); // Macro
    propertyDictionary->add(NYxy, value);


    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_gyz, "gyz"); // Macro
    propertyDictionary->add(Gyz, value);

    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_gxz, "gxz"); // Macro
    propertyDictionary->add(Gxz, value);

    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_gxy, "gxy"); // Macro
    propertyDictionary->add(Gxy, value);



    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_talphax, "talphax"); // Macro
    propertyDictionary->add(tAlphax, value);

    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_talphay, "talphay"); // Macro
    propertyDictionary->add(tAlphay, value);

    IR_GIVE_FIELD(ir, value, IFT_OrthotropicLinearElasticMaterial_talphaz, "talphaz"); // Macro
    propertyDictionary->add(tAlphaz, value);

    // check for suspicious parameters
    // ask for dependent parameters (symmetry conditions) and check if reasonable
    /*
     * nyzx = this->give(NYzx);
     * nyzy = this->give(NYzy);
     * nyyx = this->give(NYyx);
     * if ( ( nyzx < 0. ) || ( nyzx > 0.5 ) || ( nyzy < 0. ) || ( nyzy > 0.5 ) || ( nyyx < 0. ) || ( nyyx > 0.5 ) ) {
     *  _warning2("instanciateFrom: suspicious parameters", 1);
     * }
     */

    // Read local coordinate system of principal axes of ortotrophy
    // in localCoordinateSystem the unity vectors are stored
    // COLUMNWISE (this is exception, but allows faster numerical
    // implementation)
    // if you wish to align local material orientation with element, use "lcs" keyword as an element parameter

    // try to read lcs section
    triplets.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, triplets, IFT_OrthotropicLinearElasticMaterial_lcs, "lcs"); // Macro

    size = triplets.giveSize();
    if ( !( ( size == 0 ) || ( size == 6 ) ) ) {
        _warning2( "instanciateFrom: Warning: lcs in material %d is not properly defined, will be assumed as global",
                  this->giveNumber() );
    }

    if ( size == 6 ) {
        cs_type = localCS;
        double n1 = 0.0, n2 = 0.0;

        localCoordinateSystem = new FloatMatrix(3, 3);
        for ( j = 1; j <= 3; j++ ) {
            localCoordinateSystem->at(j, 1) = triplets.at(j);
            n1 += triplets.at(j) * triplets.at(j);
            localCoordinateSystem->at(j, 2) = triplets.at(j + 3);
            n2 += triplets.at(j + 3) * triplets.at(j + 3);
        }

        n1 = sqrt(n1);
        n2 = sqrt(n2);
        for ( j = 1; j <= 3; j++ ) { // normalize e1' e2'
            localCoordinateSystem->at(j, 1) /= n1;
            localCoordinateSystem->at(j, 2) /= n2;
        }

        // vector e3' computed from vector product of e1', e2'
        localCoordinateSystem->at(1, 3) =
            ( localCoordinateSystem->at(2, 1) * localCoordinateSystem->at(3, 2) -
             localCoordinateSystem->at(3, 1) * localCoordinateSystem->at(2, 2) );
        localCoordinateSystem->at(2, 3) =
            ( localCoordinateSystem->at(3, 1) * localCoordinateSystem->at(1, 2) -
             localCoordinateSystem->at(1, 1) * localCoordinateSystem->at(3, 2) );
        localCoordinateSystem->at(3, 3) =
            ( localCoordinateSystem->at(1, 1) * localCoordinateSystem->at(2, 2) -
             localCoordinateSystem->at(2, 1) * localCoordinateSystem->at(1, 2) );
    }

    // try to read ElementCS section
    if ( cs_type == unknownCS ) {
        triplets.resize(0);
        IR_GIVE_OPTIONAL_FIELD(ir, triplets, IFT_OrthotropicLinearElasticMaterial_scs, "scs"); // cs for shells.
        // first three numbers are direction of normal n - see orthoelasticmaterial.h for description
        // shellCS  - coordinate system of principal axes is specified in shell  coordinate system
        //            this is defined as follows: principal z-axis is perpendicular to mid-section
        //            x-axis is perpendicular to z-axis and normal to user specified vector n.
        //            (so x-axis is parallel to plane, with n beeing normal to this plane).
        //            y-axis is then perpendicular both to x and z axes.
        //            WARNING: this definition of cs is valid only for plates and shells
        //            when vector n is paralel to z-axis an error occurs and program is terminated.
        //
        size = triplets.giveSize();
        if ( !( ( size == 0 ) || ( size == 3 ) ) ) {
            _warning2( "instanciateFrom: scs in material %d is not properly defined, will be assumed as global",
                      this->giveNumber() );
        }

        if ( size == 3 ) {
            cs_type = shellCS;
            triplets.normalize();
            helpPlaneNormal = new FloatArray(3);
            for ( j = 1; j < 4; j++ ) {
                helpPlaneNormal->at(j) = triplets.at(j);
            }

            //
            // store normal defining help plane into row matrix
            // localCoordinateSystemmust be computed on demand from specific element
            //
        }
    } //

    if ( cs_type == unknownCS ) {
        //
        // if no cs defined assume global one
        //
        cs_type = localCS;
        localCoordinateSystem = new FloatMatrix(3, 3);
        for ( j = 1; j < 4; j++ ) {
            localCoordinateSystem->at(j, j) = 1.0;
        }
    }

    return IRRT_OK;
}


double
OrthotropicLinearElasticMaterial :: give(int aProperty, GaussPoint *gp)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
{
    if ( aProperty == NYzx ) {
        return this->give(NYxz, gp) * this->give(Ez, gp) / this->give(Ex, gp);
    }

    if ( aProperty == NYzy ) {
        return this->give(NYyz, gp) * this->give(Ez, gp) / this->give(Ey, gp);
    }

    if ( aProperty == NYyx ) {
        return this->give(NYxy, gp) * this->give(Ey, gp) / this->give(Ex, gp);
    }

    return this->Material :: give(aProperty, gp);
}


void
OrthotropicLinearElasticMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                                  MatResponseForm form,
                                                                  MatResponseMode mode,
                                                                  GaussPoint *gp,
                                                                  TimeStep *atTime)
//
// forceElasticResponse ignored - always elastic
//
{
    FloatMatrix rotationMatrix;

    this->give3dLocalMaterialStiffnessMatrix(answer, form, mode, gp, atTime);

    this->giveRotationMatrix(rotationMatrix, gp);
    answer.rotatedWith(rotationMatrix);
}


void
OrthotropicLinearElasticMaterial :: give3dLocalMaterialStiffnessMatrix(FloatMatrix &answer,
                                                                       MatResponseForm form,
                                                                       MatResponseMode mode,
                                                                       GaussPoint *gp,
                                                                       TimeStep *atTime)
{
    double eksi, nxz, nyz, nxy, nzx, nzy, nyx;
    int i, j;

    nxz = this->give(NYxz, gp);
    nyz = this->give(NYyz, gp);
    nxy = this->give(NYxy, gp);
    nzx = this->give(NYzx, gp);
    nzy = this->give(NYzy, gp);
    nyx = this->give(NYyx, gp);

    eksi = 1. - ( nxy * nyx + nyz * nzy + nzx * nxz ) - ( nxy * nyz * nzx + nyx * nzy * nxz );

    answer.resize(6, 6);
    answer.zero();
    // switched letters from original oofem -> now produces same material stiffness matrix as Abaqus method
    answer.at(1, 1) =  this->give(Ex, gp) * ( 1. - nyz * nzy ) / eksi;
    answer.at(1, 2) =  this->give(Ey, gp) * ( nxy + nxz * nzy ) / eksi;
    answer.at(1, 3) =  this->give(Ez, gp) * ( nxz + nyz * nxy ) / eksi;
    answer.at(2, 2) =  this->give(Ey, gp) * ( 1. - nxz * nzx ) / eksi;
    answer.at(2, 3) =  this->give(Ez, gp) * ( nyz + nyx * nxz ) / eksi;
    answer.at(3, 3) =  this->give(Ez, gp) * ( 1. - nyx * nxy ) / eksi;

    // define the lower triangle
    for ( i = 1; i < 4; i++ ) {
        for ( j = 1; j < i; j++ ) {
            answer.at(i, j) = answer.at(j, i);
        }
    }

    answer.at(4, 4) =  this->give(Gyz, gp);
    answer.at(5, 5) =  this->give(Gxz, gp);
    answer.at(6, 6) =  this->give(Gxy, gp);
}


void
OrthotropicLinearElasticMaterial :: giveTensorRotationMatrix(FloatMatrix &answer, GaussPoint *gp)
//
// returns [3,3] rotation matrix from local principal axes of material
// to local axes used at gp (element) level
//
{
    int elementCsFlag;
    FloatMatrix elementCs;
    StructuralElement *element = ( StructuralElement * ) gp->giveElement();

    if ( gp->giveMaterialMode() == _1dMat ) { //do not rotate 1D materials on trusses and beams
        answer.resize(3, 3);
        answer.beUnitMatrix();
        return;
    }

    elementCsFlag = element->giveLocalCoordinateSystem(elementCs);
    //
    // in localCoordinateSystem the directional cosines are stored columwise (exception)
    // in elementCs rowwise.
    //
    if ( this->cs_type == localCS ) {
        //
        // in localCoordinateSystem are stored directional cosines
        //
        if ( elementCsFlag ) {
            answer.beProductOf(elementCs, *this->localCoordinateSystem);
        } else {
            answer = *this->localCoordinateSystem;
        }
    } else if ( this->cs_type == shellCS ) {
        FloatArray elementNormal, helpx, helpy;
        localCoordinateSystem = new FloatMatrix(3, 3);

        element->computeMidPlaneNormal(elementNormal, gp);
        helpx.beVectorProductOf(*(this->helpPlaneNormal), elementNormal);
        // test if localCoordinateSystem is uniquely
        // defined by elementNormal and helpPlaneNormal
        if ( helpx.computeNorm() < ZERO_LENGTH ) {
            _error("GiveTensorRotationMatrix: element normal parallel to plane normal encountered");
        }

        helpy.beVectorProductOf(elementNormal, helpx);
        for (int i = 1; i < 4; i++ ) {
            localCoordinateSystem->at(i, 1) = helpx.at(i);
            localCoordinateSystem->at(i, 2) = helpy.at(i);
            localCoordinateSystem->at(i, 3) = elementNormal.at(i);
        }

        //
        // possible rotation about local z-axis should be considered in future
        //
        /*
         * //
         * // GiveZRotationMtrx assembles rotMtrx from cs rotated from curent about rotAngle
         * // to current cs
         * //
         * zRotMtrx = GiveZRotationMtrx (rotAngle); // rotAngle supplied by user
         * rotatedLocalCoordinateSystem = localCoordinateSystem->Times (zRotMtrx);
         * delete localCoordinateSystem;
         * localCoordinateSystem = rotatedLocalCoordinateSystem;
         */
        if ( elementCsFlag ) {
            answer.beProductOf(elementCs, *this->localCoordinateSystem);
        } else {
            answer = *this->localCoordinateSystem;
        }

        delete localCoordinateSystem;
        localCoordinateSystem = NULL;
    } else {
        _error("GiveTensorRotationMatrix - internal error no cs defined");
    }
    // t at (i,j) contains cosine of angle between elementAxis(i) and localMaterialAxis(j).
}


void
OrthotropicLinearElasticMaterial :: giveRotationMatrix(FloatMatrix &answer, GaussPoint *gp)
//
// returns [6,6] rotation matrix from local principal axes of material
// to local axes used at the gp (element) level for beams and trusses
// at element level is implemented next transformation to global cs.
//
//
{
    FloatMatrix t;
    this->giveTensorRotationMatrix(t, gp);
    this->giveStrainVectorTranformationMtrx(answer, t);
}


void
OrthotropicLinearElasticMaterial :: giveThermalDilatationVector(FloatArray &answer,
                                                                GaussPoint *gp, TimeStep *tStep)
//
// returns a FloatArray(3) of coefficients of thermal dilatation in direction
// of each (local) axis given by element lcs.
//
{
    FloatMatrix transf;
    FloatArray help(6);
    help.at(1) = this->give(tAlphax, gp);
    help.at(2) = this->give(tAlphay, gp);
    help.at(3) = this->give(tAlphaz, gp);

    this->giveRotationMatrix(transf, gp);
    answer.beProductOf(transf, help);
}


MaterialStatus *
OrthotropicLinearElasticMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates new  material status  corresponding to this class
 */
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}
} // end namespace oofem
