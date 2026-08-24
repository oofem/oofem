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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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

#include "latticelinearelastic.h"
#include "latticematstatus.h"
#include "latticestructuralmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "CrossSections/structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "Elements/LatticeElements/lattice3d.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dofmanager.h"
#include "dof.h"
#include "dofiditem.h"

#ifdef __TM_MODULE
 #include "tm/Elements/LatticeElements/latticetransportelement.h"
#endif

namespace oofem {
REGISTER_Material(LatticeLinearElastic);

// constructor which creates a dummy material without a status and without random extension interface
LatticeLinearElastic :: LatticeLinearElastic(int n, Domain *d, double e, double a1, double a2, double a3) :
    LatticeStructuralMaterial(n, d),
    eNormalMean(e),
    alphaOne(a1),
    alphaTwo(a2),
      alphaThree(a3)
{}


bool
LatticeLinearElastic :: hasMaterialModeCapability(MaterialMode mode) const
{
    return ( mode == _3dLattice );
}


void
LatticeLinearElastic :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    LatticeStructuralMaterial :: initializeFrom(ir);
    RandomMaterialExtensionInterface :: initializeFrom(ir);

    // eNormalMean = spring Em, emMacro = macroscopic Em. At least one of `e`, `em` must be set.
    const bool emGiven = ir->hasField(_IFT_LatticeLinearElastic_em);
    const bool eGiven  = ir->hasField(_IFT_LatticeLinearElastic_e);
    if ( emGiven ) {
        IR_GIVE_FIELD(ir, emMacro, _IFT_LatticeLinearElastic_em);
    }
    if ( eGiven ) {
        IR_GIVE_FIELD(ir, eNormalMean, _IFT_LatticeLinearElastic_e);
    } else if ( emGiven ) {
        eNormalMean = emMacro;
    } else {
        OOFEM_ERROR("LatticeLinearElastic: either 'e' or 'em' (or both) must be supplied");
    }
    if ( !emGiven ) {
        emMacro = eNormalMean;
    }

    // Spring ratios a1 (shear), a2 (bending), a3 (torsion); track explicit input.
    const bool a1Given = ir->hasField(_IFT_LatticeLinearElastic_a1);
    const bool a2Given = ir->hasField(_IFT_LatticeLinearElastic_a2);
    const bool a3Given = ir->hasField(_IFT_LatticeLinearElastic_a3);

    alphaOne = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaOne, _IFT_LatticeLinearElastic_a1);

    alphaTwo = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaTwo, _IFT_LatticeLinearElastic_a2);

    alphaThree = alphaTwo;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaThree, _IFT_LatticeLinearElastic_a3);

    // Nu shortcut: spring ratios from Griffiths & Mustoe (2001); explicit a1/a2/a3 win.
    nu = 0.;
    if ( ir->hasField(_IFT_LatticeLinearElastic_nu) ) {
        IR_GIVE_FIELD(ir, nu, _IFT_LatticeLinearElastic_nu);
        nuWasGiven = true;
        if ( nu >= 1. ) {
            OOFEM_ERROR("LatticeLinearElastic: nu must be < 1");
        }
        const bool applyContinuumScaling = ( emGiven && !eGiven );
        if ( applyContinuumScaling ) {
            eNormalMean = emMacro / ( 1. - nu );
        }
        if ( !a1Given ) alphaOne   = ( 1. - 3. * nu ) / ( 1. + nu );
        if ( !a2Given ) alphaTwo   = applyContinuumScaling ? 1. / ( 1. + nu ) : 1.;
        if ( !a3Given ) alphaThree = 1. / ( 2. * ( 1. + nu ) );
    }

    localRandomType = 0; //Default: No local random field
    IR_GIVE_OPTIONAL_FIELD(ir, localRandomType, _IFT_LatticeLinearElastic_localrandomtype); // Macro
    if ( localRandomType == 1 ) { //Gaussian random generator
        coefficientOfVariation = 0.;
        IR_GIVE_FIELD(ir, coefficientOfVariation, _IFT_LatticeLinearElastic_cov); // Macro
    }



    this->cAlpha = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, cAlpha, _IFT_LatticeLinearElastic_calpha);

    this->biotCoefficient = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, biotCoefficient, _IFT_LatticeLinearElastic_bio);
}

std::unique_ptr<MaterialStatus> 
LatticeLinearElastic :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<LatticeMaterialStatus>(gp);
}

MaterialStatus *
LatticeLinearElastic :: giveStatus(GaussPoint *gp) const
{
    if ( gp->hasMaterialStatus() ) {
        return static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    }
        // create a new one
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->setMaterialStatus(this->CreateStatus(gp)) );
        this->_generateStatusVariables(gp);
    return status;
    }


double
LatticeLinearElastic :: giveCouplingPressure(GaussPoint *gp, TimeStep *tStep)
{
    double waterPressure = 0.;

#ifdef __TM_MODULE
    // Coupled (staggered) run: read the pore pressure from the dual transport
    // element in the transport slave problem (couplingNumbers.at(1)). Average
    // its current-step nodal pressures (the converged transport solution of
    // this step), the same pressure state the boundary LatticeNeumannCoupling
    // uses, so both coupling paths stay consistent.
    EngngModel *master = this->domain->giveEngngModel()->giveMasterEngngModel();
    if ( master ) {
        auto *sp = dynamic_cast< StaggeredProblem * >( master );
        if ( sp ) {
            IntArray coupledModels;
            sp->giveCoupledModels(coupledModels);
            auto *se = static_cast< LatticeStructuralElement * >( gp->giveElement() );
            if ( se->giveCouplingFlag() == 1 && coupledModels.giveSize() >= 2 &&
                 coupledModels.at(2) != 0 ) {
                IntArray couplingNumbers;
                se->giveCouplingNumbers(couplingNumbers);
                if ( couplingNumbers.giveSize() >= 1 && couplingNumbers.at(1) != 0 ) {
                    Element *te = sp->giveSlaveProblem( coupledModels.at(2) )->giveDomain(1)->giveElement( couplingNumbers.at(1) );
                    int nnodes = te->giveNumberOfDofManagers();
                    for ( int i = 1; i <= nnodes; i++ ) {
                        waterPressure += te->giveDofManager(i)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
                    }
                    if ( nnodes > 0 ) {
                        waterPressure /= nnodes;
                    }
                }
            }
            return waterPressure;
        }
    }
#endif

    // Standalone run: static pore-pressure field carried on the element.
    FloatArray pressures;
    static_cast< LatticeStructuralElement * >( gp->giveElement() )->givePressures(pressures);
    for ( int i = 1; i <= pressures.giveSize(); i++ ) {
        waterPressure += pressures.at(i) / pressures.giveSize();
    }
    return waterPressure;
}


FloatArrayF< 6 >
LatticeLinearElastic :: giveLatticeStress3d(const FloatArrayF< 6 > &strain,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
{
    auto status = static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    // subtract stress independent part
    auto reducedStrain = strain;
    FloatArray indepStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
    if ( indepStrain.giveSize() > 0 ) {
        reducedStrain -= FloatArrayF< 6 >(indepStrain);
    }

    auto stiffnessMatrix = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);
    //stress are sectional forces
    auto stress = dot(stiffnessMatrix, reducedStrain);

    // Poromechanical coupling, concrete-mechanics convention: fluid pressure is
    // a stress (compression negative, suction positive), so the total normal
    // stress = effective + biot * (fluid pressure as stress). The pore pressure
    // is read live from the coupled transport problem in a staggered run, or
    // from the element's static field when run standalone.
    stress.at(1) += this->biotCoefficient * this->giveCouplingPressure(gp, tStep);

    //Set all temp values
    status->letTempLatticeStrainBe(strain);
    status->letTempLatticeStressBe(stress);
    // Expose the axial (normal) stress so a staggered Dirichlet coupling can read
    // it via Lattice2d::giveTempNormalStress (cf. latticedamage).
    status->setTempNormalLatticeStress(stress.at(1) );

    return stress;
}


void LatticeLinearElastic :: giveRandomParameters(FloatArray &param)
{
    param.resize(3);
    param.zero();
    param.at(1) = localRandomType;

    if ( localRandomType == 1 ) { //Gaussian
        param.at(2) = coefficientOfVariation;
    } else {
        OOFEM_ERROR("Error: Unknown local random type:\n randomtype 1 = Gaussian\n");
    }
}


Interface *
LatticeLinearElastic :: giveInterface(InterfaceType type)
{
    return nullptr;
}


FloatMatrixF< 6, 6 >
LatticeLinearElastic :: give3dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime) const
{
    //Needed to make sure that status exists before random values are requested for elastic stiffness. Problem is that gp->giveMaterialStatus does not check if status exist already
  static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

  //Reduce Young's modulus based on temperature
  double reductionFactor =1.;
  if(this->tCrit !=0.){
    reductionFactor = computeTemperatureReductionFactor(gp,atTime,VM_Total);
  }
  
  // Shell mode: a1y (comp 2, out-of-plane shear) = G = Em/(2(1+nu)); the (1-nu)
  // cancels the eNormalMean 1/(1-nu). a1z (comp 3) keeps the Griffiths-Mustoe ratio.
  double a3eff = this->alphaThree;
  double a1y = this->alphaOne;
  const double a1z = this->alphaOne;
  if ( nuWasGiven ) {
      Lattice3d *elem = dynamic_cast< Lattice3d * >( gp->giveElement() );
      if ( elem != nullptr && elem->isShellElement() ) {
          a3eff = 1. / ( 1. + nu );
          a1y = ( 1. - nu ) / ( 2. * ( 1. + nu ) );
      }
  }

    FloatArrayF< 6 >d = {
        1.,
      a1y,
      a1z,
      a3eff,
      this->alphaTwo,
      this->alphaTwo
    };

    return diag(d * this->give(eNormal_ID, gp) * this->eNormalMean * reductionFactor);
}

FloatMatrixF< 3, 3 >
LatticeLinearElastic :: give2dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime) const
{
  //Needed to make sure that status exists before random values are requested for elastic stiffness. Problem is that gp->giveMaterialStatus does not check if status exist already
  static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

  FloatArrayF< 3 >d = {
        1.,
        this->alphaOne, // shear
        this->alphaTwo, // bending
    };

    return diag(d * this->give(eNormal_ID, gp) * this->eNormalMean);
}


FloatArrayF< 6 >
LatticeLinearElastic :: giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const
//
// returns a FloatArray(6) of initial strain vector
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    double alpha = this->give(tAlpha, gp);

    //Option to add a eigendisplacement instead of strain
    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    alpha += this->cAlpha / length;

    return {
               alpha, 0., 0., 0., 0., 0.
    };
}


double
LatticeLinearElastic :: give(int aProperty, GaussPoint *gp) const
{
   this->giveStatus(gp);
  
    double answer;
    if ( RandomMaterialExtensionInterface :: give(aProperty, gp, answer) ) {
        if ( answer < 0.1 ) { //Introduce cut off to avoid numerical problems
            answer = 0.1;
        } else if ( answer > 10 ) {
            answer = 10;
        }
        return answer;
    } else if ( aProperty == eNormal_ID ) {
        return 1.;
    } else if ( aProperty == 'E' ) {
        return this->eNormalMean;
    } else {
        return LatticeStructuralMaterial :: give(aProperty, gp);
    }
}
}
