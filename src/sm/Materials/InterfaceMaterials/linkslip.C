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

#include "linkslip.h"
//#include "linearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
//#include "CrossSections/structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "Elements/structuralelement.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"
// #ifdef __TM_MODULE
//  #include "latticetransportelement.h"
//  #include "latticetransmat.h"
//  #include "latticelammat.h"
//  #include "pore.h"
//  #include "discretetransportproblem.h"
// #endif

namespace oofem {
REGISTER_Material(LinkSlip);

/// constructor which creates a dummy material without a status and without random extension interface
LinkSlip :: LinkSlip(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{
}


void
LinkSlip :: initializeFrom(InputRecord &ir)
{
    StructuralInterfaceMaterial :: initializeFrom(ir);;
    
    //axial stiffness
    IR_GIVE_FIELD(ir, kNormal, _IFT_LinkSlip_kn); // Macro

    //Ratio of lateral to axial stiffness
    IR_GIVE_OPTIONAL_FIELD(ir, kLateral, _IFT_LinkSlip_kl); // Macro

    //Two models available
    //0 default. Linear elastic-perfect plastic
    //1 bond slip according first part of CEB model and then constant
    //2 bond slip according to CEB model
    IR_GIVE_OPTIONAL_FIELD(ir, type, _IFT_LinkSlip_type); // Macro
    
    IR_GIVE_FIELD(ir, tauMax, _IFT_LinkSlip_t0); // Macro

    if ( type == 1 || type == 2 ) {
        IR_GIVE_FIELD(ir, s1, _IFT_LinkSlip_s1);
        IR_GIVE_FIELD(ir, alpha, _IFT_LinkSlip_alpha);
    }
    
    if ( type == 2 ) {
        IR_GIVE_FIELD(ir, s2, _IFT_LinkSlip_s2);
        IR_GIVE_FIELD(ir, s3, _IFT_LinkSlip_s3);
        IR_GIVE_FIELD(ir, tauFinal, _IFT_LinkSlip_tf);
    }

    if ( type == 1 || type == 2 ) {
        if ( (type == 1 || type == 2) && this->kNormal < this->tauMax/this->s1 ) {
            this->kNormal = this->tauMax/this->s1;
            OOFEM_WARNING("Parameter kN adjusted");
        }
    }
}


MaterialStatus *
LinkSlip :: CreateStatus(GaussPoint *gp) const
{
    return new LinkSlipStatus(gp);
}
  
double
LinkSlip :: evaluateBondStress(const double kappa) const
{
    if ( this->type == 0 ) { //elastic-perfectly plastic
        return this->tauMax;
    } else if ( this->type == 1 ) { //modified bond model      
        if ( kappa <= 0. ){
            return 0.;
        }
        if ( kappa <= s1 ) {
            return tauMax * pow(kappa/s1, alpha);
        }
        return tauMax;
    } else if ( this->type == 2 ) {
        if ( kappa <= 0. ) {
            return 0.;
        }
        if ( kappa <= s1 ) {
            return tauMax * pow(kappa/s1, alpha);
        }
        if ( kappa <= s2 ) {
            return tauMax;
        }
        if ( kappa <= s3 ) {
            return tauMax - (tauMax-tauFinal) * (kappa-s2) / (s3-s2);
        }
        return tauFinal;
    } else { //unknown type
        OOFEM_ERROR("Unknown bond model type. Type should be 0, 1 or 2.");
    }
    
    // return 0.; //Should not be here. 
}


FloatArrayF<3>
LinkSlip :: giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< LinkSlipStatus * >( this->giveStatus(gp) );

    //For axial (first) component, strain has the meanig of slip. Stress is traction.
    
    const auto &oldTraction = status->giveTraction();
    const auto &oldJump = status->giveJump();

    //evaluate tempKappa (no elastic strain in axial direction)
    double tempKappa = status->giveKappa() + fabs(jump.at(1)-oldJump.at(1));

    FloatArrayF<3> traction;

    //trial stress in axial direction
    traction.at(1) = oldTraction.at(1) + (jump.at(1)-oldJump.at(1))*this->kNormal;
    
    
    double f = fabs(traction.at(1)) - evaluateBondStress(tempKappa);

    if ( f > 0 ) { //plastic response.
        //Reduced stress by increasing plastic strain.
        traction.at(1) = evaluateBondStress(tempKappa);;
    }

    //Compute the lateral stress components
    for ( int i = 2; i <= 3; i++ ) { // only diagonal terms matter
        traction.at(i) = this->kLateral * jump.at(i);
    }

    //Set temp values in status needed for dissipation
    status->letTempKappaBe(tempKappa);
    status->letTempJumpBe(jump);
    status->letTempTractionBe(traction);

    return traction;
}

Interface *
LinkSlip :: giveInterface(InterfaceType type)
{
    return nullptr;
}

FloatMatrixF<3,3>
LinkSlip :: give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    return diag<3>({this->kNormal, this->kLateral, this->kLateral});
}

LinkSlipStatus :: LinkSlipStatus(GaussPoint *g) :  StructuralInterfaceMaterialStatus(g)
{
}
  

void
LinkSlipStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    StructuralInterfaceMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
}
  

void
LinkSlipStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    StructuralInterfaceMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "  jump ");
    for ( auto &val : this->jump ) {
        fprintf(file, " %.4e", val );
    }

    fprintf(file, "\n              traction ");
    for ( auto &val : this->traction ) {
        fprintf(file, " %.4e", val );
    }
    fprintf(file, "\n");
    
    fprintf(file, "kappa %.8e\n", this->kappa);
    return;
}
  
void
LinkSlipStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    // save parent class status
  
    StructuralInterfaceMaterialStatus :: saveContext(stream, mode);
  
    // write a raw data
    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    
}


void
LinkSlipStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralInterfaceMaterialStatus :: restoreContext(stream, mode);
  
    // read raw data
    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
 
void
LinkSlipStatus :: updateYourself(TimeStep *atTime)
{
    StructuralInterfaceMaterialStatus :: updateYourself(atTime);
    this->kappa = this->tempKappa;
}

 
int
LinkSlip :: giveIPValue(FloatArray &answer,
                        GaussPoint *gp,
                        InternalStateType type,
                        TimeStep *atTime)
{
    auto status = static_cast< LinkSlipStatus * >( this->giveStatus(gp) );
    if ( type == IST_InterfaceJump ) {
        answer.resize(3);
        answer.zero();
        answer = status->giveJump();
        return 1;
    } else if ( type == IST_InterfaceTraction ) {
        answer.resize(3);
        answer.zero();
        answer = status->giveTraction();
        return 1;
    } else {
        return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, atTime);
    }
//     return Material :: giveIPValue(answer, gp, type, atTime);
}
}
