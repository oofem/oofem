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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#include "beam3dsubsoil.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "fei3dlinelin.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Beam3dSubsoil);

FEI3dLineLin Beam3dSubsoil :: interp;

Beam3dSubsoilElementInterface::Beam3dSubsoilElementInterface()
{
}
  
Beam3dSubsoil :: Beam3dSubsoil(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{
    
    numberOfGaussPoints = 3;
    numberOfDofMans = 2;
    mode = _Winkler;
    local_formulation = true;

    
}

FEInterpolation *
Beam3dSubsoil:: giveInterpolation() const
{
  return &interp;
}
  
void
Beam3dSubsoil :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 5) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

// Returns the strain-displacement matrix {B} of the receiver,
// evaluated at gp.
void
Beam3dSubsoil :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
  Beam3dSubsoilElementInterface *ei = (Beam3dSubsoilElementInterface*) this->domain->giveElement(this->master)->giveInterface(Beam3dSubsoilElementInterfaceType);
  FloatMatrix T;

  if (ei) {
    


    if (mode == _Winkler) {
      FloatMatrix Ng;
      answer.clear();
      ei->B3SSI_getNMatrix (answer, gp); // interpolation matrix (in local c.s.)
      if (!local_formulation) { // formulated (stifness given) in global c.s.
	ei->B3SSI_getNodalGtoLRotationMatrix(T);
        // transform N to global form
	FloatMatrix h(answer);
        answer.beTProductOf(T, h); // N = T^T N
      }
    } else {
      OOFEM_ERROR ("Unsupported mode");
    }
  } else {
    OOFEM_ERROR ("master does not implement Beam3dSubsoilElementInterfaceType");
  }
}

bool
Beam3dSubsoil::computeGtoLRotationMatrix(FloatMatrix &answer)
{
  if (local_formulation) { // formulated (stifness given) in local c.s.
    Beam3dSubsoilElementInterface* ei = (Beam3dSubsoilElementInterface*) (this->domain->giveElement(this->master)->giveInterface(Beam3dSubsoilElementInterfaceType));
    ei->B3SSI_getGtoLRotationMatrix(answer);
    return true;
  } else {
    // global formulation
    answer.clear();
    return false;
  }
}


void
Beam3dSubsoil :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
  answer.clear();
}


void
Beam3dSubsoil :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  answer.beDiagonal(this->springConstants);
}


IRResultType
Beam3dSubsoil :: initializeFrom(InputRecord *ir)
{
  IRResultType result;
  this->numberOfGaussPoints = 4;
  result = StructuralElement :: initializeFrom(ir);
  IR_GIVE_FIELD(ir, master, _IFT_Beam3dSubsoil_master);
  IR_GIVE_FIELD(ir, springConstants, _IFT_Beam3dSubsoil_springConstants);
  if (springConstants.giveSize() != 6) {
    OOFEM_ERROR ("Size of springConstants parameter should be 6");
    return IRRT_BAD_FORMAT;
  }

    
  int val = 0;
  IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_Beam3dSubsoil_localFormulation);
  local_formulation = val; // default is false = global formulation

  return IRRT_OK;
}


void
Beam3dSubsoil :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {D_u, D_v, D_w, R_u, R_v, R_w};
}


double
Beam3dSubsoil :: computeVolumeAround(GaussPoint *gp)
{
  return ((Beam3dSubsoilElementInterface*) (this->domain->giveElement(this->master)->giveInterface(Beam3dSubsoilElementInterfaceType)))->B3SSI_computeVolumeAround(gp);
}

void
Beam3dSubsoil :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
  OOFEM_ERROR("Mass matrix not provided");
}


void
Beam3dSubsoil::giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // This function can be quite costly to do inside the loops when one has many slave dofs.
  FloatMatrix K;
  FloatArray u;
  this->computeVectorOf(VM_Total, tStep, u);
  this->computeStiffnessMatrix(K, TangentStiffness, tStep);
  answer.beProductOf(K,u);
}


} // end namespace oofem
