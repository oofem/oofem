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
 *               Copyright (C) 1993 - 2023   Borek Patzak
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


#include "Contact/thermalcontactbc.h"
#include "set.h"
#include "domain.h"
#include "node.h"
#include "floatmatrix.h"
#include "function.h"

namespace oofem {


void ThermalContactBoundaryCondition :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    ContactBoundaryCondition::initializeFrom(ir);


    IR_GIVE_FIELD(ir, this->gapConductanceFunctionNumber, _IFT_ThermalContactBoundaryCondition_gapConductanceFunctionNumber);
   
    /*IR_GIVE_FIELD(ir, this->gapConductance, _IFT_ThermalContactBoundaryCondition_gapConductance);
    IR_GIVE_FIELD(ir, this->conductanceDistance, _IFT_ThermalContactBoundaryCondition_conductanceDistance);
    if(this->gapConductance.giveSize() != 2 && this->conductanceDistance.giveSize() !=2) {
      OOFEM_ERROR("Wrong size of gap conductance input parameter");
    }
    */
}


  
double
ThermalContactBoundaryCondition :: evaluateGapConductance(double gap)
{
  /*
  if( gap >this->conductanceDistance.at(2)) {
    return 0;
  } else {
    return (this->gapConductance.at(1) + ((this->gapConductance.at(2)-this->gapConductance.at(1))/(this->conductanceDistance.at(2)-this->conductanceDistance.at(1)))*gap );
  }
  */
  return domain->giveFunction(this->gapConductanceFunctionNumber)->evaluateAtTime(gap);
   
}


void
ThermalContactBoundaryCondition :: computeTangentFromContact(FloatMatrix &answer, ContactPair *contactPair, TimeStep *tStep)
{

    auto gap = contactPair->giveNormalGap();
    if(gap < 0) {
      gap = 0;
    }
    auto normal = contactPair->giveNormalVector();
    auto A = normal.computeNorm();
   
    //
    FloatMatrix N;
    contactPair->computeNmatrix(N);
    //
    answer.beTProductOf(N,N);
    //
    answer.times(this->evaluateGapConductance(gap) * A);
}


void
ThermalContactBoundaryCondition :: computeInternalForcesFromContact(FloatArray &answer, ContactPair *contactPair, TimeStep *tStep)
{
  auto gap = contactPair->giveNormalGap();
  if(gap < 0) {
    gap = 0;
  }
  auto normal = contactPair->giveNormalVector();
  auto A = normal.computeNorm();

  //
  FloatMatrix N;
  contactPair->computeNmatrix(N);
  //
  FloatArray T, Ts;
  contactPair->giveMasterContactPoint()->giveUnknownVector(T, {T_f}, VM_Total, tStep);
  contactPair->giveSlaveContactPoint()->giveUnknownVector(Ts, {T_f}, VM_Total, tStep);
  //  
  T.copySubVector(Ts,T.giveSize()+1);
  Ts.beProductOf(N,T);
  //
  answer.beRowOf(N,1);
  answer.times( Ts.at(1) * this->evaluateGapConductance(gap) * A);
}



      
  
} // namespace oofem
