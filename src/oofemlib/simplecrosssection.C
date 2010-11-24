/* $Header: /home/cvs/bp/oofem/oofemlib/src/simplecrosssection.C,v 1.11 2003/04/14 16:00:47 bp Exp $ */
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

//   file SIMPLECROSSSECTION.CC

#include "simplecrosssection.h"
#include "gausspnt.h"
#include "material.h"
#include "structuralmaterial.h"
#include "flotarry.h"

namespace oofem {

void
SimpleCrossSection ::  giveRealStresses(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                        const FloatArray &totalStrain, TimeStep *tStep)
//
// this function returns a real stresses corresponding to
// given totalStrain according to stressStrain mode stored
// in each gp.
// IMPORTANT:
//
{
    //FloatArray *reducedStressIncrement;
    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );

    if ( mat->hasMaterialModeCapability(mode) ) {
        StructuralCrossSection :: giveRealStresses(answer, form, gp, totalStrain, tStep);
        return;
    } else {
        _error("giveRealStresses : unsupported mode");
    }

    return;
}


void
SimpleCrossSection :: giveCharMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                        MatResponseForm form,
                                                        MatResponseMode rMode,
                                                        GaussPoint *gp,
                                                        StructuralMaterial *mat,
                                                        TimeStep *tStep)
//
// only interface to material class, forcing returned matrix to be in reduced form.
//
{
    this->giveMaterialStiffnessMatrixOf(answer, form, rMode, gp, mat, tStep);
    return;
}



void
SimpleCrossSection :: giveMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                    MatResponseForm form,
                                                    MatResponseMode rMode,
                                                    GaussPoint *gp,
                                                    StructuralMaterial *mat,
                                                    TimeStep *tStep)
//
// may be only simple interface to material class, forcing returned matrix to be in reduced form.
// otherwise special methods called to obtain required stiffness from 3d case.
//
{
    // Material *mat = gp->giveElement()->giveMaterial();
    if ( mat->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
        mat->giveCharacteristicMatrix(answer, form, rMode, gp, tStep);
        return;
    } else {
        _error("GiveMaterialStiffnessMatrixOf: unsupported StressStrainMode");
    }
}

void
SimpleCrossSection :: giveFullCharacteristicVector(FloatArray &answer,
                                                   GaussPoint *gp,
                                                   const FloatArray &strainVector)
//
// returns full 3d strain vector from strainVector in reducedMode
// based on StressStrainMode in gp
// strainVector {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// enhaced method in order to support cases with integral bending (2dplate, 3dshell..)
// in such cases full strain vector has the form:
// strainVectorShell {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy, kappa_x,kappa_y,kappa_y,kappa_yz,kappa_xz,kappa_xy}
//
// enhanced method in order to support 3dbeam elements
// in such cases full strain vector has the form:
// strainVector {eps_x, gamma_xz, gamma_xy, \der{phi_x}{x}, kappa_y, kappa_z}
//
// enhance support also for PlaneStressRot case with full strain vector of form
// {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy,(omega_xy-(dv/dx-du/dy)*0.5)}
//
//
{
    // FloatArray *answer;
    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( gp->giveMaterial() );
    IntArray indx;
    int i, j, answerSize = 0;


    //if (mode ==  _3dShell) {answer =  strainVector; return ;}
    //if (mode ==  _3dBeam)  {answer =  strainVector; return ;}
    if ( ( mode == _3dShell ) || ( mode ==  _3dBeam ) || ( mode == _2dPlate ) || ( mode == _2dBeam ) || ( mode == _PlaneStressRot ) ) {
        if ( ( mode == _3dShell ) || ( mode ==  _3dBeam ) || ( mode == _2dPlate ) || ( mode == _2dBeam ) ) {
            answerSize = 12;
        }

        if ( mode == _PlaneStressRot ) {
            answerSize = 7;
        }

        //answer = new FloatArray (answerSize);
        answer.resize(answerSize);
        answer.zero();

        mat->giveStressStrainMask( indx, ReducedForm, gp->giveMaterialMode() );
        for ( i = 1; i <= indx.giveSize(); i++ ) {
            if ( ( j = indx.at(i) ) ) {
                answer.at(j) = strainVector.at(i);
            }
        }

        //delete indx;
        return;
    } else {
        StructuralCrossSection :: giveFullCharacteristicVector(answer, gp, strainVector);
        return;
    }
}


void
SimpleCrossSection :: giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                      const FloatArray &charVector3d)
//
// returns reduced stressVector or strainVector from full 3d vector reduced
// to vector required by gp->giveStressStrainMode()
//
{
    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );
    IntArray indx;
    int size = charVector3d.giveSize();
    int i, j;
    //FloatArray* answer;

    if ( ( mode == _3dShell ) || ( mode == _3dBeam ) || ( mode == _2dPlate ) || ( mode == _2dBeam ) ) {
        if ( size != 12 ) {
            _error("giveReducedCharacteristicVector - charVector3d size mismatch");
            exit(0);
        }

        mat->giveStressStrainMask( indx, ReducedForm, gp->giveMaterialMode() );
        //answer = new FloatArray (indx->giveSize());
        answer.resize( indx.giveSize() );
        answer.zero();

        for ( i = 1; i <= indx.giveSize(); i++ ) {
            if ( ( j = indx.at(i) ) ) {
                answer.at(i) = charVector3d.at(j);
            }
        }

        //delete indx;
        return;
    } else if ( mode == _3dBeam ) {
        if ( size != 6 ) {
            _error("giveReducedCharacteristicVector - charVector3d size mismatch");
            exit(0);
        }

        answer =  charVector3d;
        return;
    } else if ( mode == _PlaneStressRot ) {
        if ( size != 7 ) {
            _error("giveReducedCharacteristicVector - charVector3d size mismatch");
            exit(0);
        }

        mat->giveStressStrainMask( indx, ReducedForm, gp->giveMaterialMode() );
        //answer = new FloatArray (indx->giveSize());
        answer.resize( indx.giveSize() );
        answer.zero();

        for ( i = 1; i <= indx.giveSize(); i++ ) {
            if ( ( j = indx.at(i) ) ) {
                answer.at(i) = charVector3d.at(j);
            }
        }

        //delete indx;
        return;
    } else {
        StructuralCrossSection :: giveReducedCharacteristicVector(answer, gp, charVector3d);
        return;
    }
}



IRResultType
SimpleCrossSection :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    double value;

    this->CrossSection :: initializeFrom(ir);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_SimpleCrossSection_thick, "thick"); // Macro
    propertyDictionary->add(CS_Thickness, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_SimpleCrossSection_width, "width"); // Macro
    propertyDictionary->add(CS_Width, value);

    if ( ir->hasField(IFT_SimpleCrossSection_area, "area") ) {
        IR_GIVE_FIELD(ir, value, IFT_SimpleCrossSection_area, "area"); // Macro
        propertyDictionary->add(CS_Area, value);
    }

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_SimpleCrossSection_iy, "iy"); // Macro
    propertyDictionary->add(CS_InertiaMomentY, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_SimpleCrossSection_iz, "iz"); // Macro
    propertyDictionary->add(CS_InertiaMomentZ, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_SimpleCrossSection_ik, "ik"); // Macro
    propertyDictionary->add(CS_TorsionMomentX, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_SimpleCrossSection_shearcoeff, "beamshearcoeff"); // Macro
    propertyDictionary->add(CS_BeamShearCoeff, value);

    return IRRT_OK;
}


int
SimpleCrossSection :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    CrossSection :: giveInputRecordString(str, keyword);
    sprintf( buff, " thick %e width %e area %e iy %e iz %e ik %e beamshearcoeff %e",
            this->give(CS_Thickness), this->give(CS_Width), this->give(CS_Area),
            this->give(CS_InertiaMomentY), this->give(CS_InertiaMomentZ), this->give(CS_TorsionMomentX),
            this->give(CS_BeamShearCoeff) );
    str += buff;

    return 1;
}



double
SimpleCrossSection :: give(CrossSectionProperty aProperty)
{
    printf("property = %d",aProperty);
    double value = 0.0;

     if ( aProperty == CS_Area ) {
        if ( propertyDictionary->includes(CS_Area) ) {
            return propertyDictionary->at(CS_Area);
        } else {
            return this->give(CS_Thickness) * this->give(CS_Width);
        }
    }
    else if ( propertyDictionary->includes(aProperty) ) {
        value = propertyDictionary->at(aProperty);
    } else {
        OOFEM_ERROR("SimpleCrossSection :: give: property not defined");
    }

    return value;
}

void
SimpleCrossSection :: computeStressIndependentStrainVector(FloatArray &answer,
                                                           GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
//
// returns initial strain vector induced by stress independent effects
// like temperatue or shrinkage.
// takes into account form of load vector assumed by engngModel (Incremental or Total Load form).
//
{
    StructuralMaterial *mat = ( StructuralMaterial * ) gp->giveElement()->giveMaterial();
    /*
     * MaterialMode matmode = gp-> giveMaterialMode ();
     * FloatArray et, e0, fullAnswer;
     * double thick, width;
     *
     * if ((matmode == _2dBeam) || (matmode == _3dBeam) || (matmode == _3dShell) || (matmode == _2dPlate)) {
     *
     * StructuralElement *elem = (StructuralElement*)gp->giveElement();
     * elem -> computeResultingIPTemperatureAt (et, stepN, gp, mode);
     * FloatArray redAnswer;
     *
     * if (et.giveSize() == 0) {answer.resize(0); return ;}
     * if (et.giveSize() < 1) {
     * _error ("computeStressIndependentStrainVector - Bad format of TemperatureLoad");
     * exit (1);
     * }
     * mat->giveThermalDilatationVector (e0, gp,stepN);
     *
     * if (matmode == _2dBeam) {
     * answer.resize (3);
     * answer.zero();
     * answer.at(1) = e0.at(1) * (et.at(1)- mat->giveReferenceTemperature());
     * if (et.giveSize() > 1) {
     * thick = this->give(THICKNESS);
     * answer.at(2) = e0.at(1) * et.at(2)/ thick;   // kappa_x
     * }
     * } else if (matmode == _3dBeam) {
     * answer.resize (6);
     * answer.zero();
     *
     * answer.at(1) = e0.at(1) * (et.at(1)- mat->giveReferenceTemperature());
     * if (et.giveSize() > 1) {
     * thick = this->give(THICKNESS);
     * width = this->give(WIDTH);
     * answer.at(5) = e0.at(1) * et.at(2)/ thick;   // kappa_y
     * if (et.giveSize() > 2)
     * answer.at(6) = e0.at(1) * et.at(3)/ width;   // kappa_z
     * }
     * } else if (matmode == _2dPlate) {
     *
     * if (et.giveSize() > 1) {
     * answer.resize (5);
     * answer.zero();
     *
     * thick = this->give(THICKNESS);
     * if (et.giveSize() > 1) {
     * answer.at(1) = e0.at(1) * et.at(2)/ thick;   // kappa_x
     * answer.at(2) = e0.at(2) * et.at(2)/ thick;   // kappa_y
     * }
     * }
     * } else if (matmode == _3dShell) {
     * answer.resize (8);
     * answer.zero();
     *
     * answer.at(1) = e0.at(1) * (et.at(1)- mat->giveReferenceTemperature());
     * answer.at(2) = e0.at(2) * (et.at(1)- mat->giveReferenceTemperature());
     * if (et.giveSize() > 1) {
     * thick = this->give(THICKNESS);
     * answer.at(4) = e0.at(1) * et.at(2)/ thick;   // kappa_x
     * answer.at(5) = e0.at(2) * et.at(2)/ thick;   // kappa_y
     * }
     * } else _error ("Unsupported material mode");
     * } else {
     * mat->computeStressIndependentStrainVector (answer, gp, stepN, mode);
     * }
     */
    mat->computeStressIndependentStrainVector(answer, gp, stepN, mode);
}

} // end namespace oofem
