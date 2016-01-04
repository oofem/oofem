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

#include "Materials/linearelasticmaterial.h"
#include "Materials/ortholinearelasticmaterial.h"
#include "../sm/Elements/structuralelement.h"
#include "material.h"
#include "../sm/Materials/structuralms.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
#define ZERO_LENGTH 1.e-6

REGISTER_Material(OrthotropicLinearElasticMaterial);

IRResultType
OrthotropicLinearElasticMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    double value;
    int size;
    FloatArray triplets;


    result = LinearElasticMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_ex);
    propertyDictionary.add(Ex, value);

    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_ey);
    propertyDictionary.add(Ey, value);

    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_ez);
    propertyDictionary.add(Ez, value);


    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_nyyz);
    propertyDictionary.add(NYyz, value);

    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_nyxz);
    propertyDictionary.add(NYxz, value);

    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_nyxy);
    propertyDictionary.add(NYxy, value);


    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_gyz);
    propertyDictionary.add(Gyz, value);

    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_gxz);
    propertyDictionary.add(Gxz, value);

    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_gxy);
    propertyDictionary.add(Gxy, value);



    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_talphax);
    propertyDictionary.add(tAlphax, value);

    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_talphay);
    propertyDictionary.add(tAlphay, value);

    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_talphaz);
    propertyDictionary.add(tAlphaz, value);

    // check for suspicious parameters
    // ask for dependent parameters (symmetry conditions) and check if reasonable
    /*
     * nyzx = this->give(NYzx);
     * nyzy = this->give(NYzy);
     * nyyx = this->give(NYyx);
     * if ( ( nyzx < 0. ) || ( nyzx > 0.5 ) || ( nyzy < 0. ) || ( nyzy > 0.5 ) || ( nyyx < 0. ) || ( nyyx > 0.5 ) ) {
     *  OOFEM_WARNING("suspicious parameters", 1);
     * }
     */

    // Read local coordinate system of principal axes of ortotrophy
    // in localCoordinateSystem the unity vectors are stored
    // COLUMNWISE (this is exception, but allows faster numerical
    // implementation)
    // if you wish to align local material orientation with element, use "lcs" keyword as an element parameter

    // try to read lcs section
    triplets.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, triplets, _IFT_OrthotropicLinearElasticMaterial_lcs);

    size = triplets.giveSize();
    if ( !( ( size == 0 ) || ( size == 6 ) ) ) {
        OOFEM_WARNING("Warning: lcs in material %d is not properly defined, will be assumed as global",
                  this->giveNumber() );
    }

    if ( size == 6 ) {
        cs_type = localCS;
        double n1 = 0.0, n2 = 0.0;

        localCoordinateSystem = new FloatMatrix(3, 3);
        for ( int j = 1; j <= 3; j++ ) {
            localCoordinateSystem->at(j, 1) = triplets.at(j);
            n1 += triplets.at(j) * triplets.at(j);
            localCoordinateSystem->at(j, 2) = triplets.at(j + 3);
            n2 += triplets.at(j + 3) * triplets.at(j + 3);
        }

        n1 = sqrt(n1);
        n2 = sqrt(n2);
        for ( int j = 1; j <= 3; j++ ) { // normalize e1' e2'
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
        triplets.clear();
        IR_GIVE_OPTIONAL_FIELD(ir, triplets, _IFT_OrthotropicLinearElasticMaterial_scs); // cs for shells.
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
            OOFEM_WARNING("scs in material %d is not properly defined, will be assumed as global",
                      this->giveNumber() );
        }

        if ( size == 3 ) {
            cs_type = shellCS;
            triplets.normalize();
            helpPlaneNormal = new FloatArray(triplets);

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
        localCoordinateSystem->beUnitMatrix();
    }

    return IRRT_OK;
}


void
OrthotropicLinearElasticMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    Material :: giveInputRecord(input);


    input.setField(propertyDictionary.at(Ex), _IFT_OrthotropicLinearElasticMaterial_ex);
    input.setField(propertyDictionary.at(Ey), _IFT_OrthotropicLinearElasticMaterial_ey);
    input.setField(propertyDictionary.at(Ez), _IFT_OrthotropicLinearElasticMaterial_ez);

    input.setField(propertyDictionary.at(NYyz), _IFT_OrthotropicLinearElasticMaterial_nyyz);
    input.setField(propertyDictionary.at(NYxz), _IFT_OrthotropicLinearElasticMaterial_nyxz);
    input.setField(propertyDictionary.at(NYxy), _IFT_OrthotropicLinearElasticMaterial_nyxy);

    input.setField(propertyDictionary.at(Gyz), _IFT_OrthotropicLinearElasticMaterial_gyz);
    input.setField(propertyDictionary.at(Gxz), _IFT_OrthotropicLinearElasticMaterial_gxz);
    input.setField(propertyDictionary.at(Gxy), _IFT_OrthotropicLinearElasticMaterial_gxy);

    input.setField(propertyDictionary.at(tAlphax), _IFT_OrthotropicLinearElasticMaterial_talphax);
    input.setField(propertyDictionary.at(tAlphay), _IFT_OrthotropicLinearElasticMaterial_talphay);
    input.setField(propertyDictionary.at(tAlphaz), _IFT_OrthotropicLinearElasticMaterial_talphaz);


    ///@todo Should add optional arguments:
    // _IFT_OrthotropicLinearElasticMaterial_lcs
    // _IFT_OrthotropicLinearElasticMaterial_scs
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
                                                                  MatResponseMode mode,
                                                                  GaussPoint *gp,
                                                                  TimeStep *tStep)
//
// forceElasticResponse ignored - always elastic
//
{
    FloatMatrix rotationMatrix;

    this->give3dLocalMaterialStiffnessMatrix(answer, mode, gp, tStep);

    this->giveRotationMatrix(rotationMatrix, gp);
    answer.rotatedWith(rotationMatrix);
}


void
OrthotropicLinearElasticMaterial :: give3dLocalMaterialStiffnessMatrix(FloatMatrix &answer,
                                                                       MatResponseMode mode,
                                                                       GaussPoint *gp,
                                                                       TimeStep *tStep)
{
    double eksi, nxz, nyz, nxy, nzx, nzy, nyx;

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
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < i; j++ ) {
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
    StructuralElement *element = static_cast< StructuralElement * >( gp->giveElement() );

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
            answer.beProductOf(elementCs, * this->localCoordinateSystem);
        } else {
            answer = * this->localCoordinateSystem;
        }
    } else if ( this->cs_type == shellCS ) {
        FloatArray elementNormal, helpx, helpy;
        localCoordinateSystem = new FloatMatrix(3, 3);

        element->computeMidPlaneNormal(elementNormal, gp);
        helpx.beVectorProductOf(* ( this->helpPlaneNormal ), elementNormal);
        // test if localCoordinateSystem is uniquely
        // defined by elementNormal and helpPlaneNormal
        if ( helpx.computeNorm() < ZERO_LENGTH ) {
            OOFEM_ERROR("element normal parallel to plane normal encountered");
        }

        helpy.beVectorProductOf(elementNormal, helpx);
        for ( int i = 1; i < 4; i++ ) {
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
            answer.beProductOf(elementCs, * this->localCoordinateSystem);
        } else {
            answer = * this->localCoordinateSystem;
        }

        delete localCoordinateSystem;
        localCoordinateSystem = NULL;
    } else {
        OOFEM_ERROR("internal error no cs defined");
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

} // end namespace oofem
