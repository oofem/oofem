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

#include "../sm/Elements/linedistributedspring.h"
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
#include "load.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(LineDistributedSpring);

FEI3dLineLin LineDistributedSpring :: interp_lin;

LineDistributedSpring :: LineDistributedSpring(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this),
    SPRNodalRecoveryModelInterface()
{
    numberOfGaussPoints = 2;
    numberOfDofMans = 2;
}


FEInterpolation *
LineDistributedSpring :: giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}


FEInterpolation *
LineDistributedSpring :: giveInterpolation() const { return & interp_lin; }


void
LineDistributedSpring :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 5) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
LineDistributedSpring :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
{
  OOFEM_ERROR("Body load not supported, use surface load instead");
}


void
LineDistributedSpring :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [3x3] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FloatArray n;
    int ndofs = this->dofs.giveSize();

    this->interp_lin.evalN( n, gp->giveNaturalCoordinates(),  FEIElementGeometryWrapper(this) );

    answer.resize(ndofs, ndofs*2);
    answer.zero();

    for (int idof=1; idof<=ndofs; idof++) {
      answer.at(idof, idof) = n.at(1); 
      answer.at(idof, ndofs+idof) = n.at(2); 
    }
}


void
LineDistributedSpring :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
  int ndofs = this->dofs.giveSize();
  answer.resize(ndofs);

  for (int idof=1; idof<=ndofs; idof++) {
    answer.at(idof) = strain.at(idof)*this->springStiffnesses.at(idof);
  }
  //this->giveStructuralCrossSection()->giveGeneralizedStress_DistributedSpring(answer, gp, strain, tStep);
}


void
LineDistributedSpring :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  //this->giveStructuralCrossSection()->give2dPlateSubSoilStiffMtrx(answer, rMode, gp, tStep);
  int ndofs = this->dofs.giveSize();
  answer.resize(ndofs, ndofs);
  answer.beDiagonal(this->springStiffnesses);
}


void
LineDistributedSpring::giveInternalForcesVector(FloatArray &answer,
                                                TimeStep *tStep, int useUpdatedGpRecord)
{
  FloatArray u;
  FloatMatrix k;

  this->computeVectorOf(VM_Total, tStep, u);
  this->computeStiffnessMatrix(k, TangentStiffness, tStep);
  answer.beProductOf(k, u);

}


  
IRResultType
LineDistributedSpring :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD (ir, dofs, _IFT_LineDistributedSpring_Dofs);
    IR_GIVE_FIELD (ir, springStiffnesses, _IFT_LineDistributedSpring_Stifnesses);

    if (dofs.giveSize() != springStiffnesses.giveSize()) {
      OOFEM_ERROR ("dofs and k params size mismatch");
    }
    // from element
    return StructuralElement::initializeFrom(ir);
}

int
LineDistributedSpring::checkConsistency()
{
  // skip StructuralElement consistency as there is checjk for structurak cross section capability
  return Element::checkConsistency();
}

void LineDistributedSpring :: printOutputAt(FILE *File, TimeStep *tStep)
{
  FloatArray stress, strain;

  fprintf(File, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );
    
  for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
    fprintf(File, "\tGP 1.%d :\n", gp->giveNumber());
    this->computeStrainVector(strain, gp, tStep);
    this->computeStressVector(stress, strain, gp, tStep);
    fprintf(File, "\t\tStrain");
    for (int i=1; i<=strain.giveSize(); i++) fprintf (File, " %e", strain.at(i));
    fprintf(File, "\n\t\tStress");
    for (int i=1; i<=stress.giveSize(); i++) fprintf (File, " %e", stress.at(i));
    fprintf(File, "\n");
  }

}


  
void
LineDistributedSpring :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
   answer = this->dofs;
}


double
LineDistributedSpring :: computeVolumeAround(GaussPoint *gp)
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs( this->interp_lin.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ * weight;
}


void
LineDistributedSpring :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
  OOFEM_ERROR("Mass matrix not provided");
}


int
LineDistributedSpring :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
  return StructuralElement :: giveIPValue(answer, gp, type, tStep);
}

Interface *
LineDistributedSpring :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    }

    return NULL;
}

void
LineDistributedSpring :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(2);
    for ( int i = 1; i < 3; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
LineDistributedSpring :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i < 3; i++ ) {
        if ( pap == this->giveNode(i)->giveNumber() ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}


} // end namespace oofem
