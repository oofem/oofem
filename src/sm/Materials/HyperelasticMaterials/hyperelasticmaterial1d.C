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
 *               Copyright (C) 1993 - 2020   Borek Patzak
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

#include "hyperelasticmaterial1d.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"



namespace oofem {
REGISTER_Material(HyperelasticMaterial1d);

HyperelasticMaterial1d::HyperelasticMaterial1d(int n, Domain *d) : StructuralMaterial(n, d)
{ }

FloatArrayF< 1 >
HyperelasticMaterial1d::giveFirstPKStressVector_1d(const FloatArrayF< 1 > &vF, GaussPoint *gp, TimeStep *tStep) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    FloatArrayF<1> vP;
    if(this->hyperelasticMaterialType == HEMT_Biot) {
      vP.at(1) = this->E * (vF.at(1) - 1.);
    } else if(this->hyperelasticMaterialType == HEMT_StVenantKirchhoff) {
      vP.at(1) = vF.at(1) * this->E * 0.5 * (vF.at(1) * vF.at(1) - 1.);
    } else {
      OOFEM_ERROR("Unknow material type");      
    }
    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(vP);
    //
    return vP;
}



FloatMatrixF< 1, 1 >
HyperelasticMaterial1d::give1dStressStiffnessMatrix_dPdF(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double e = this->E;
    if(this->hyperelasticMaterialType == HEMT_Biot) {
      return {e};
    } else if(this->hyperelasticMaterialType == HEMT_StVenantKirchhoff) {
      FloatArray vF(status->giveTempFVector() );
      FloatMatrix d(1,1);
      d.at(1,1)= e * 0.5 * (vF.at(1) * vF.at(1) - 1.) +  e * vF.at(1)*vF.at(1);
      return d;
    } else {
      OOFEM_ERROR("Unknow material type");      
    }
    

}




  

MaterialStatus *
HyperelasticMaterial1d::CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(gp);
}


void
HyperelasticMaterial1d::initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    int hyperTypeRecord = 0; // default 
    IR_GIVE_OPTIONAL_FIELD(ir, hyperTypeRecord, _IFT_HyperelasticMaterial1d_type);
    // specify the type of hyperelasticmaterial
    if ( hyperTypeRecord == 0 ) {
        this->hyperelasticMaterialType = HEMT_Biot;
    } else if ( hyperTypeRecord == 1 ) {
        this->hyperelasticMaterialType = HEMT_StVenantKirchhoff;
    } else {
      throw ValueInputException(ir, _IFT_HyperelasticMaterial1d_type, "Unknown hyperelastic material type");
    }

    
    IR_GIVE_FIELD(ir, this->E, _IFT_HyperelasticMaterial1d_E);
    
}
} // end namespace oofem
