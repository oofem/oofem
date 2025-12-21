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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include <vector>
#include <cstdio>
#include <iostream>

#include "errorcheckingexportmodule.h"
#include "engngm.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "timestep.h"
#include "classfactory.h"
#include "dof.h"
#include "oofemtxtinputrecord.h"
#ifdef _USE_XML
	#include "xmlinputrecord.h"
	#include "xmldatareader.h"
#endif
#include "mathfem.h"
#ifdef __SM_MODULE
 #include "sm/EngineeringModels/structengngmodel.h"
 #include "sm/Elements/Beams/beam2d.h"
 #include "sm/Elements/Beams/beam3d.h"
#endif

namespace oofem {
REGISTER_ExportModule(ErrorCheckingExportModule)


bool
ErrorCheckingRule :: checkValue(double computedValue)
{
    return computedValue >= value - tolerance && computedValue <= value + tolerance;
}


NodeErrorCheckingRule :: NodeErrorCheckingRule(const std :: string &line, double tol) :
  ErrorCheckingRule(tol)
{
/*
  char unknown;
  std :: string  kwd;
  OOFEMTXTInputRecord rec (0, line), *ptr = &rec;
  rec.giveRecordKeywordField(kwd);
  if (kwd == "#NODE") {
    IR_GIVE_FIELD (ptr,tstep, "tStep");
    IR_GIVE_OPTIONAL_FIELD (ptr,tsubstep, "tStepVer");
    IR_GIVE_FIELD (ptr,number, "number");
    IR_GIVE_FIELD (ptr,dofid, "dofid");
    IR_GIVE_FIELD (ptr,unknown, "unknown");
    IR_GIVE_FIELD (ptr,value, "value");
    IR_GIVE_OPTIONAL_FIELD (ptr,tolerance, "tolerance");
  } else {
    OOFEM_ERROR("Something wrong in the error checking rule: %s\n", line.c_str());
  }
*/
    char unknown;
    int ret = std :: sscanf(line.c_str(), "#NODE tStep %d number %d dof %d unknown %c value %le tolerance %le",
                  &tstep, & number, & dofid, & unknown, & value, & tolerance);
    if ( ret < 2 ) {
      ret = std :: sscanf(line.c_str(), "#NODE tStep %d tStepVer %d number %d dof %d unknown %c value %le tolerance %le",
                          &tstep, &tsubstep, & number, & dofid, & unknown, & value, & tolerance);
    }
    if ( ret < 5 ) {
        OOFEM_ERROR("Something wrong in the error checking rule: %s\n", line.c_str());
    }

    if ( unknown == 'd' ) {
        mode = VM_Total;
    } else if ( unknown == 'v' ) {
        mode = VM_Velocity;
    } else if ( unknown == 'a' ) {
        mode = VM_Acceleration;
    } else {
        OOFEM_ERROR("Can't recognize unknown '%c'", unknown);
    }
}

NodeErrorCheckingRule :: NodeErrorCheckingRule(InputRecord& ir, double tol): ErrorCheckingRule(tol){
    std::string unknown;
    ir.giveField(tstep,"tStep");
    ir.giveOptionalField(tsubstep,"tStepVer");
    ir.giveField(number,"number");
    ir.giveField(dofid,"dof");
    ir.giveField(unknown,"unknown");
    ir.giveField(value,"value");
    ir.giveOptionalField(tolerance,"tolerance");
    if ( unknown == "d" ) {
        mode = VM_Total;
    } else if ( unknown == "v" ) {
        mode = VM_Velocity;
    } else if ( unknown == "a" ) {
        mode = VM_Acceleration;
    } else {
        OOFEM_ERROR("Can't recognize unknown '%s' (must be one of: d, v, a)", unknown.c_str());
    }
}

bool 
NodeErrorCheckingRule :: getValue(double &answer, Domain *domain, TimeStep *tStep)
{
    DofManager *dman = domain->giveGlobalDofManager(number);
    if ( !dman ) {
        if ( domain->giveEngngModel()->isParallel() ) {
            return false;
        } else {
            OOFEM_WARNING("Dof manager %d not found.", number);
            return false;
        }
    }

    if ( dman->giveParallelMode() == DofManager_remote || dman->giveParallelMode() == DofManager_null ) {
        return false;
    }

    Dof *dof = dman->giveDofWithID(dofid);
    answer = dof->giveUnknown(mode, tStep);
    return true;
}


bool
NodeErrorCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    // Rule doesn't apply yet.
    if ( tStep->giveNumber() != tstep || tStep->giveVersion() != tsubstep )  {
        return true;
    }

    DofManager *dman = domain->giveGlobalDofManager(number);
    if ( !dman ) {
        if ( domain->giveEngngModel()->isParallel() ) {
            return true;
        } else {
            OOFEM_WARNING("Dof manager %d not found.", number);
            return false;
        }
    }

    if ( dman->giveParallelMode() == DofManager_remote || dman->giveParallelMode() == DofManager_null ) {
        return true;
    }

    Dof *dof = dman->giveDofWithID(dofid);

    double dmanValue = dof->giveUnknown(mode, tStep);
    bool check = checkValue(dmanValue);
    if ( !check ) {
        OOFEM_WARNING("Check failed in %s: tstep %d, node %d, dof %d, mode %d:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      domain->giveEngngModel()->giveOutputBaseFileName().c_str(), tstep, number, dofid, mode,
                      dmanValue, value, fabs(dmanValue-value), tolerance );
    }
    return check;
}


InternalElementDofManErrorCheckingRule::InternalElementDofManErrorCheckingRule (const std :: string &line, double tol) :
    ErrorCheckingRule(tol)
{
    char unknown;
    int ret = std :: sscanf(line.c_str(), "#ELEMENTNODE tStep %d number %d dofman %d dof %d unknown %c value %le tolerance %le",
                  &tstep, & number, &idofman, & dofid, & unknown, & value, & tolerance);
    if ( ret < 5 ) {
        OOFEM_ERROR("Something wrong in the error checking rule: %s\n", line.c_str());
    }

    if ( unknown == 'd' ) {
        mode = VM_Total;
    } else if ( unknown == 'v' ) {
        mode = VM_Velocity;
    } else if ( unknown == 'a' ) {
        mode = VM_Acceleration;
    } else {
        OOFEM_ERROR("Can't recognize unknown '%c'", unknown);
    }
}

InternalElementDofManErrorCheckingRule::InternalElementDofManErrorCheckingRule(InputRecord& ir, double tol): ErrorCheckingRule(tol) {
    std::string unknown;
    ir.giveField(tstep,"tStep");
    ir.giveField(number,"number");
    ir.giveField(idofman,"dofman");
    ir.giveField(dofid,"dof");
    ir.giveField(unknown,"unknown");
    ir.giveField(value,"value");
    if(unknown=="d") mode=VM_Total;
    else if(unknown=="v") mode=VM_Velocity;
    else if(unknown=="a") mode=VM_Acceleration;
    else OOFEM_ERROR("Can't recognize unknown '%s' (must be one of: d, v, a)",unknown.c_str());
}


bool
InternalElementDofManErrorCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    // Rule doesn't apply yet.
    if ( tStep->giveNumber() != tstep || tStep->giveVersion() != tsubstep )  {
        return true;
    }

    Element *element = domain->giveGlobalElement(number);
    if ( !element ) {
        if ( domain->giveEngngModel()->isParallel() ) {
            return true;
        } else {
            OOFEM_WARNING("Element %d not found.", number);
            return false;
        }
    }
    if ( element->giveParallelMode() != Element_local ) {
        return true;
    }
    DofManager *dman = element->giveInternalDofManager(idofman);
    if ( !dman ) {
        OOFEM_WARNING("Internal DofManager %d on element %d not found.", idofman, number);
        return false;
    }

    Dof *dof = dman->giveDofWithID(dofid);

    double dmanValue = dof->giveUnknown(mode, tStep);
    bool check = checkValue(dmanValue);
    if ( !check ) {
        OOFEM_WARNING("Check failed in %s: tstep %d, node %d, dof %d, mode %d:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      domain->giveEngngModel()->giveOutputBaseFileName().c_str(), tstep, number, dofid, mode,
                      dmanValue, value, fabs(dmanValue-value), tolerance );
    }
    return check;
}

bool
InternalElementDofManErrorCheckingRule :: getValue(double &answer, Domain *domain, TimeStep *tStep)
{
    Element *element = domain->giveGlobalElement(number);
    if ( !element ) {
        if ( domain->giveEngngModel()->isParallel() ) {
            return false;
        } else {
            OOFEM_WARNING("Element %d not found.", number);
            return false;
        }
    }
    if ( element->giveParallelMode() != Element_local ) {
        return false;
    }
    DofManager *dman = element->giveInternalDofManager(idofman);
    if ( !dman ) {
        OOFEM_WARNING("Internal DofManager %d on element %d not found.", idofman, number);
        return false;
    }

    Dof *dof = dman->giveDofWithID(dofid);

    answer = dof->giveUnknown(mode, tStep);
    return true;
}

ElementErrorCheckingRule :: ElementErrorCheckingRule(const std :: string &line, double tol) :
    ErrorCheckingRule(tol), irule(0)
{
    int istnum;
    int ret = std :: sscanf(line.c_str(), "#ELEMENT tStep %d number %d irule %d gp %d keyword %d component %d value %le tolerance %le",
                  &tstep, & number, & irule, & gpnum, & istnum, & component, & value, & tolerance);
    if ( ret < 2 ) {
      ret = std :: sscanf(line.c_str(), "#ELEMENT tStep %d tStepVer %d number %d irule %d gp %d keyword %d component %d value %le tolerance %le",
                          &tstep, &tsubstep, & number, & irule, & gpnum, & istnum, & component, & value, & tolerance);
      if (ret < 4) {
        ret = std :: sscanf(line.c_str(), "#ELEMENT tStep %d tStepVer %d  number %d gp %d keyword %d component %d value %le tolerance %le",
                              &tstep, &tsubstep, & number, & gpnum, & istnum, & component, & value, & tolerance);
        }
    } else if ( ret < 3 ) {
        ret = std :: sscanf(line.c_str(), "#ELEMENT tStep %d number %d gp %d keyword %d component %d value %le tolerance %le",
                    &tstep, & number, & gpnum, & istnum, & component, & value, & tolerance);
    }
    if ( ret < 6 ) {
        OOFEM_ERROR("Something wrong in the error checking rule: %s\n", line.c_str());
    }
    ist = (InternalStateType)istnum;
}

ElementErrorCheckingRule :: ElementErrorCheckingRule(InputRecord& ir, double tol): ErrorCheckingRule(tol) {
    ir.giveField(tstep,"tStep");
    ir.giveOptionalField(tsubstep,"tStepVer");
    ir.giveField(number,"number");
    ir.giveOptionalField(irule,"irule");
    ir.giveField(gpnum,"gp");
    int istnum;
    ir.giveField(istnum,"keyword");
    ist = (InternalStateType)istnum;
    ir.giveField(component,"component");
    ir.giveField(value,"value");
    ir.giveOptionalField(tolerance,"tolerance");
}

bool
ElementErrorCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    // Rule doesn't apply yet.
     if ( tStep->giveNumber() != tstep || tStep->giveVersion() != tsubstep )  {
        return true;
    }

    FloatArray ipval;
    Element *element = domain->giveGlobalElement(number);
    // std::cerr<<"Element "<<number<<" @ "<<(void*)element<<std::endl;
    if ( !element ) {
        if ( domain->giveEngngModel()->isParallel() ) {
            return true;
        } else {
            OOFEM_WARNING("Element %d not found.", number);
            return false;
        }
    }
        if ( element->giveParallelMode() != Element_local ) {
            return true;
        }

    // note! GPs are numbered from 0 internally, but written with 1-index, inconsistent!
    GaussPoint *gp = element->giveIntegrationRule(irule)->getIntegrationPoint(gpnum-1);
    element->giveIPValue(ipval, gp, ist, tStep);

    if ( component > ipval.giveSize() || component < 1 ) {
        OOFEM_WARNING("Check failed in %s: element %d, gpnum %d, ist %d, component %d:\n"
                      "Component not found!", 
                      domain->giveEngngModel()->giveOutputBaseFileName().c_str(), number, gpnum, ist, component);
        ipval.printYourself();
        return false;
    }

    double elementValue = ipval.at(component);
    bool check = checkValue(elementValue);
    if ( !check ) {
        OOFEM_WARNING("Check failed in %s: tstep %d, element %d, gpnum %d, ist %d, component %d:\n"                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      domain->giveEngngModel()->giveOutputBaseFileName().c_str(), tstep, number, gpnum, ist, component,
                      elementValue, value, fabs(elementValue-value), tolerance );
        ipval.printYourself();
    }
    return check;
}


bool
ElementErrorCheckingRule :: getValue(double & answer, Domain *domain, TimeStep *tStep)
{
    FloatArray ipval;
    Element *element = domain->giveGlobalElement(number);
    if ( !element ) {
        if ( domain->giveEngngModel()->isParallel() ) {
            return false;
        } else {
            OOFEM_WARNING("Element %d not found.", number);
            return false;
        }
    }
        if ( element->giveParallelMode() != Element_local ) {
            return false;
        }

    // note! GPs are numbered from 0 internally, but written with 1-index, inconsistent!
    GaussPoint *gp = element->giveIntegrationRule(irule)->getIntegrationPoint(gpnum-1);
    element->giveIPValue(ipval, gp, ist, tStep);

    if ( component > ipval.giveSize() || component < 1 ) {
        OOFEM_WARNING("Check failed in %s: element %d, gpnum %d, ist %d, component %d:\n"
                      "Component not found!", 
                      domain->giveEngngModel()->giveOutputBaseFileName().c_str(), number, gpnum, ist, component);
        ipval.printYourself();
        return false;
    }

    answer = ipval.at(component);
    return true;
}


BeamElementErrorCheckingRule :: BeamElementErrorCheckingRule(const std :: string &line, double tol) :
    ErrorCheckingRule(tol)
{
    int istnum;
    int ret = std :: sscanf(line.c_str(), "#BEAM_ELEMENT tStep %d number %d keyword %d component %d value %le tolerance %le",
                  &tstep, & number, & istnum, & component, & value, & tolerance);
    if ( ret < 2 ) {
      ret = std :: sscanf(line.c_str(), "#BEAM_ELEMENT tStep %d tStepVer %d number %d keyword %d component %d value %le tolerance %le",
                  &tstep, &tsubstep, & number, & istnum, & component, & value, & tolerance);
    }
    if ( ret < 5 ) {
        OOFEM_ERROR("Something wrong in the error checking rule: %s\n", line.c_str());
    }
    ist = (BeamElementErrorCheckingRule::BeamElementValueType)istnum;
}

BeamElementErrorCheckingRule :: BeamElementErrorCheckingRule(InputRecord& ir, double tol): ErrorCheckingRule(tol) {
    ir.giveField(tstep,"tStep");
    ir.giveOptionalField(tsubstep,"tStepVer");
    ir.giveField(number,"number");
    int istnum;
    ir.giveField(istnum,"keyword");
    ist = (BeamElementErrorCheckingRule::BeamElementValueType)istnum;
    ir.giveField(component,"component");
    ir.giveField(value,"value");
    ir.giveOptionalField(tolerance,"tolerance");
}

bool
BeamElementErrorCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    // Rule doesn't apply yet.
    if ( tStep->giveNumber() != tstep || tStep->giveVersion() != tsubstep )  {
        return true;
    }

    FloatArray val;
    Element *element = domain->giveGlobalElement(number);
    if ( !element ) {
        if ( domain->giveEngngModel()->isParallel() ) {
            return true;
        } else {
            OOFEM_WARNING("Element %d not found.", number);
            return false;
        }
    }
    if ( element->giveParallelMode() != Element_local ) {
        return true;
    }

    if (ist == BET_localEndDisplacement) {
        element->computeVectorOf(VM_Total, tStep, val);
    } else if (ist ==  BET_localEndForces) {
#ifdef __SM_MODULE 
        if(Beam2d* b = dynamic_cast<Beam2d*>(element)) b->giveEndForcesVector(val, tStep);
        else if(Beam3d* b = dynamic_cast<Beam3d*>(element)) b->giveEndForcesVector(val, tStep);
        else {
            OOFEM_WARNING("Element %d has no beam interface.", number);
            return false;
        }
#else
        OOFEM_WARNING("Element %d has no beam interface.", number);
        return false;
#endif
    }

    if ( component > val.giveSize() || component < 1 ) {
        OOFEM_WARNING("Check failed in %s: beam_element %d, ist %d, component %d:\n"
                      "Component not found!",
                      domain->giveEngngModel()->giveOutputBaseFileName().c_str(), number, ist, component);
        val.printYourself();
        return false;
    }

    double elementValue = val.at(component);
    bool check = checkValue(elementValue);
    if ( !check ) {
        OOFEM_WARNING("Check failed in %s: tstep %d, beam_element %d, ist %d, component %d:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      domain->giveEngngModel()->giveOutputBaseFileName().c_str(), tstep, number, ist, component,
                      elementValue, value, fabs(elementValue-value), tolerance );
        val.printYourself();
    }
    return check;
}


bool
BeamElementErrorCheckingRule :: getValue(double& answer, Domain *domain, TimeStep *tStep)
{
    FloatArray val;
    Element *element = domain->giveGlobalElement(number);
    if ( !element ) {
        if ( domain->giveEngngModel()->isParallel() ) {
            return false;
        } else {
            OOFEM_WARNING("Element %d not found.", number);
            return false;
        }
    }
    if ( element->giveParallelMode() != Element_local ) {
        return false;
    }

    if (ist == BET_localEndDisplacement) {
        element->computeVectorOf(VM_Total, tStep, val);
    } else if (ist ==  BET_localEndForces) {
#ifdef __SM_MODULE 
        if(Beam2d* b = dynamic_cast<Beam2d*>(element)) b->giveEndForcesVector(val, tStep);
        else if(Beam3d* b = dynamic_cast<Beam3d*>(element)) b->giveEndForcesVector(val, tStep);
        else {
            OOFEM_WARNING("Element %d has no beam interface.", number);
            return false;
        }
#else
        OOFEM_WARNING("Element %d has no beam interface.", number);
        return false;
#endif
    }

    if ( component > val.giveSize() || component < 1 ) {
        OOFEM_WARNING("Check failed in %s: beam_element %d, ist %d, component %d:\n"
                      "Component not found!",
                      domain->giveEngngModel()->giveOutputBaseFileName().c_str(), number, ist, component);
        val.printYourself();
        return false;
    }

    answer = val.at(component);
    return true;
}


ReactionErrorCheckingRule :: ReactionErrorCheckingRule(const std :: string &line, double tol) :
    ErrorCheckingRule(tol)
{
    int ret = std :: sscanf(line.c_str(), "#REACTION tStep %d number %d dof %d value %le tolerance %le", 
                  &tstep, & number, & dofid, & value, & tolerance);
    if ( ret < 2 ) {
      ret = std :: sscanf(line.c_str(), "#REACTION tStep %d tStepVer %d number %d dof %d value %le tolerance %le", 
                  &tstep, &tsubstep, & number, & dofid, & value, & tolerance);
    }
    if ( ret < 4 ) {
        OOFEM_ERROR("Something wrong in the error checking rule: %s\n", line.c_str());
    }
}

ReactionErrorCheckingRule :: ReactionErrorCheckingRule(InputRecord& ir, double tol): ErrorCheckingRule(tol) {
    ir.giveField(tstep,"tStep");
    ir.giveOptionalField(tsubstep,"tStepVer");
    ir.giveField(number,"number");
    ir.giveField(dofid,"dof");
    ir.giveField(value,"value");
    ir.giveOptionalField(tolerance,"tolerance");
}

bool
ReactionErrorCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    // Rule doesn't apply yet.
    if ( tStep->giveNumber() != tstep || tStep->giveVersion() != tsubstep )  {
        return true;
    }

#ifdef __SM_MODULE
    //EngngModel *emodel = domain->giveEngngModel();
    StructuralEngngModel *emodel = static_cast< StructuralEngngModel* >( domain->giveEngngModel() );
    ///@todo These would be easier to use if they just returned some more useable structure. Perhaps embed inside a PrimaryField, or just some std::pair/tuple.
    ///Now we instead have two calls and have to decipher the information from that.
    //std :: map< int, std :: map< int, double > > reactionForces; // reactionForces[nodeNumber][dofid] like this
    //emodel->computeReaction(std :: reactionForces, tStep, domain->giveNumber());
    FloatArray reactionForces;
    IntArray restrDofMans, restrDofs, eqn;
    emodel->buildReactionTable(restrDofMans, restrDofs, eqn, tStep, domain->giveNumber());
    emodel->computeReaction(reactionForces, tStep, domain->giveNumber());

    bool found = false;
    int index;
    for ( index = 1; index <= restrDofs.giveSize(); ++index ) {
        if ( restrDofs.at(index) == dofid && domain->giveNode(restrDofMans.at(index))->giveLabel() == number ) {
            found = true;
            break;
        }
    }
    if ( !found ) {
        if ( domain->giveEngngModel()->isParallel() ) {
            return true;
        } else {
            OOFEM_WARNING("Reaction force node: %d dof: %d not found.", number, dofid);
            return false;
        }
    }

    double reactionForce = reactionForces.at(index);
    bool check = checkValue(reactionForce);
    if ( !check ) {
        OOFEM_WARNING("Check failed in %s: tstep %d, reaction forces number %d, dof %d:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      domain->giveEngngModel()->giveOutputBaseFileName().c_str(), tstep, number, dofid,
                      reactionForce, value, fabs(reactionForce-value), tolerance );
    }
    return check;
#else
    OOFEM_WARNING("Reaction forces only supported for structural problems yet");
    return false;
#endif
}

bool
ReactionErrorCheckingRule :: getValue(double &answer, Domain *domain, TimeStep *tStep)
{

#ifdef __SM_MODULE
    //EngngModel *emodel = domain->giveEngngModel();
    StructuralEngngModel *emodel = static_cast< StructuralEngngModel* >( domain->giveEngngModel() );
    ///@todo These would be easier to use if they just returned some more useable structure. Perhaps embed inside a PrimaryField, or just some std::pair/tuple.
    ///Now we instead have two calls and have to decipher the information from that.
    //std :: map< int, std :: map< int, double > > reactionForces; // reactionForces[nodeNumber][dofid] like this
    //emodel->computeReaction(std :: reactionForces, tStep, domain->giveNumber());
    FloatArray reactionForces;
    IntArray restrDofMans, restrDofs, eqn;
    emodel->buildReactionTable(restrDofMans, restrDofs, eqn, tStep, domain->giveNumber());
    emodel->computeReaction(reactionForces, tStep, domain->giveNumber());

    bool found = false;
    int index;
    for ( index = 1; index <= restrDofs.giveSize(); ++index ) {
        if ( restrDofs.at(index) == dofid && domain->giveNode(restrDofMans.at(index))->giveLabel() == number ) {
            found = true;
            break;
        }
    }
    if ( !found ) {
        if ( domain->giveEngngModel()->isParallel() ) {
            return false;
        } else {
            OOFEM_WARNING("Reaction force node: %d dof: %d not found.", number, dofid);
            return false;
        }
    }

    answer = reactionForces.at(index);
    return true;
#else
    OOFEM_WARNING("Reaction forces only supported for structural problems yet");
    return false;
#endif
}

LoadLevelErrorCheckingRule :: LoadLevelErrorCheckingRule(const std :: string &line, double tol) :
    ErrorCheckingRule(tol)
{
    int ret = std :: sscanf(line.c_str(), "#LOADLEVEL tStep %d value %le tolerance %le", 
                  &tstep, & value, & tolerance);
    if ( ret < 2 ) {
        OOFEM_ERROR("Something wrong in the error checking rule: %s\n", line.c_str());
    }
}

LoadLevelErrorCheckingRule :: LoadLevelErrorCheckingRule(InputRecord& ir, double tol): ErrorCheckingRule(tol) {
    ir.giveField(tstep,"tStep");
    ir.giveField(value,"value");
    ir.giveOptionalField(tolerance,"tolerance");
}


bool
LoadLevelErrorCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    // Rule doesn't apply yet.
    if ( tStep->giveNumber() != tstep ) {
        return true;
    }

    double loadLevel = domain->giveEngngModel()->giveLoadLevel();
    bool check = checkValue(loadLevel);
    if ( !check ) {
        OOFEM_WARNING("Check failed in %s: tstep %d, load level:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      domain->giveEngngModel()->giveOutputBaseFileName().c_str(), tstep,
                      loadLevel, value, fabs(loadLevel-value), tolerance );
    }
    return check;
}

bool
LoadLevelErrorCheckingRule :: getValue(double& answer, Domain *domain, TimeStep *tStep)
{

    answer = domain->giveEngngModel()->giveLoadLevel();
    return true;
}

EigenValueErrorCheckingRule :: EigenValueErrorCheckingRule(const std :: string &line, double tol) :
    ErrorCheckingRule(tol)
{
    int ret = std :: sscanf(line.c_str(), "#EIGVAL tStep %d EigNum %d value %le tolerance %le", 
                  &tstep, & number, & value, & tolerance);
    if ( ret < 3 ) {
        OOFEM_ERROR("Something wrong in the error checking rule: %s\n", line.c_str());
    }
}

EigenValueErrorCheckingRule :: EigenValueErrorCheckingRule(InputRecord& ir, double tol): ErrorCheckingRule(tol) {
    ir.giveField(tstep,"tStep");
    ir.giveField(number,"EigNum");
    ir.giveField(value,"value");
    ir.giveOptionalField(tolerance,"tolerance");
}

bool
EigenValueErrorCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    // Rule doesn't apply yet.
    if ( tStep->giveNumber() != tstep ) {
        return true;
    }

    double eig = domain->giveEngngModel()->giveEigenValue(number);
    bool check = checkValue(eig);
    if ( !check ) {
        OOFEM_WARNING("Check failed in %s: tstep %d, eigen value %d:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      domain->giveEngngModel()->giveOutputBaseFileName().c_str(), tstep, number,
                      eig, value, fabs(eig-value), tolerance );
    }
    return check;
}
bool
EigenValueErrorCheckingRule :: getValue(double&answer, Domain *domain, TimeStep *tStep)
{
    answer = domain->giveEngngModel()->giveEigenValue(number);
    return true;
}

TimeCheckingRule :: TimeCheckingRule(const std :: string &line, double tol) :
    ErrorCheckingRule(tol)
{
}

TimeCheckingRule::TimeCheckingRule(InputRecord& ir, double tol): ErrorCheckingRule(tol){ }

bool
TimeCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    return true;
}
bool
TimeCheckingRule :: getValue(double&answer, Domain *domain, TimeStep *tStep)
{
    answer = tStep->giveTargetTime();
    return true;
}



////////////////////////////////////////////////////////////////////////////////////////////////////


ErrorCheckingExportModule :: ErrorCheckingExportModule(int n, EngngModel *e) : ExportModule(n, e)
{
}

void
ErrorCheckingExportModule :: initializeFrom(InputRecord &ir)
{
    ExportModule :: initializeFrom(ir);

    allPassed = true;
    this->errorCheckingRules.clear();

    filename = std::string("");

    if ( ir.hasField(_IFT_ErrorCheckingExportModule_filename) ) {
        IR_GIVE_FIELD(ir, this->filename, _IFT_ErrorCheckingExportModule_filename);
    }
    else {
        filename = emodel->giveReferenceFileName();
    }
    #ifdef _USE_XML
        /* we need to cast to XMLInputRecord just to get the reader object */
        XMLInputRecord* xmlrec=dynamic_cast<XMLInputRecord*>(&ir);
        if(xmlrec) this->readRulesFromRecords(*(xmlrec->giveReader()),ir);
        else
    #endif
    this->readRulesFromTextFile(ir);

    this->writeIST.clear();
    writeChecks = ir.hasField(_IFT_ErrorCheckingExportModule_writeIST);
    if ( writeChecks ) {
        IR_GIVE_FIELD(ir, this->writeIST, _IFT_ErrorCheckingExportModule_writeIST);
    }

    if ( errorCheckingRules.size() == 0 && !writeChecks ) {
        OOFEM_WARNING("No rules found (possibly wrong file or syntax).");
    }

    IR_GIVE_OPTIONAL_FIELD(ir, this->extractorMode, _IFT_ErrorCheckingExportModule_extractormode);
}

void ErrorCheckingExportModule::readRulesFromTextFile(InputRecord& ir){
    // Reads all the rules;
    std :: ifstream inputStream(this->filename);
    if ( !inputStream ) {
        throw ValueInputException(ir, _IFT_ErrorCheckingExportModule_filename, "Couldn't open file");
    }
    double tol = 0.;
    if ( this->scanToErrorChecks(inputStream,  tol) ) {
        for (;;) {
            std :: unique_ptr< ErrorCheckingRule > rule = this->giveErrorCheck(inputStream, tol);
            if ( !rule ) {
                break;
            }
            errorCheckingRules.push_back(std :: move(rule));
        }
    }
}

void ErrorCheckingExportModule::readRulesFromRecords(DataReader& dr, InputRecord& ir){
    double tol=1e-6;
    ir.giveOptionalField(tol,"tolerance");
    DataReader::GroupRecords ruleRecs=dr.giveGroupRecords("",/*whatever*/DataReader::IR_elemRec,-1);
    for(auto& rir: ruleRecs){
        std::string n;
        rir.giveRecordKeywordField(n);
        // std::cerr<<"Check rule of type "<<n<<std::endl;
        std::unique_ptr<ErrorCheckingRule> rule;
        if (n=="NODE") { rule=std::make_unique<NodeErrorCheckingRule>(rir,tol); }
        if (n=="ELEMENT") { rule=std::make_unique<ElementErrorCheckingRule>(rir,tol); }
        if (n=="REACTION") { rule=std::make_unique<ReactionErrorCheckingRule>(rir,tol); }
        if (n=="BEAM_ELEMENT") { rule=std::make_unique<BeamElementErrorCheckingRule>(rir,tol); }
        if (n=="EIGVAL") { rule=std::make_unique<EigenValueErrorCheckingRule>(rir,tol); }
        if (n=="LOADLEVEL") { rule=std::make_unique<LoadLevelErrorCheckingRule>(rir,tol); }
        if (n=="TIME") { rule=std::make_unique<TimeCheckingRule>(rir,tol); }
        if (n=="ELEMENTNODE") { rule=std::make_unique<InternalElementDofManErrorCheckingRule>(rir,tol); }
        if(rule) errorCheckingRules.push_back(std::move(rule));
        else { std::cerr<<"No rule for "<<n<<" created (not yet implemented for XML?)."<<std::endl; }
    }
}


void
ErrorCheckingExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
#if 0
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }
#endif

    if (!this->extractorMode) {
        // Error checking rules are hardcoded to domain 1 always.
        Domain *domain = emodel->giveDomain(1);

        OOFEM_LOG_INFO("Checking rules...\n");
        for ( auto &rule: this->errorCheckingRules ) {
            this->allPassed &= rule->check(domain, tStep);
        }

        if ( !tStep->isNotTheLastStep() ) {
            if ( !this->allPassed ) {
                OOFEM_ERROR("Rule not passed, exiting with error");
            }
        }

        if ( this->writeChecks ) {
            this->writeCheck(domain, tStep);
        }
    } else {

        // Error checking rules are hardcoded to domain 1 always.
        Domain *domain = emodel->giveDomain(1);

        double value;
        for ( auto &rule: this->errorCheckingRules ) {
            bool result = rule->getValue(value, domain, tStep);
            if (result) {
                fprintf (outputFile, "%+12.8e ", value);
            } else {
                fprintf (outputFile, "%12s","-");
            }
        }
        fprintf(outputFile, "\n");
    }
}

void
ErrorCheckingExportModule::initialize()
{
    if (this->extractorMode) {
        char filename [100];    
        sprintf( filename, "%s.m%d", this->emodel->giveOutputBaseFileName().c_str(), this->number);
        this->outputFile = fopen(filename, "w");
    }
}

void ErrorCheckingExportModule::terminate()
{
    if (this->extractorMode) {
        fclose (this->outputFile);   
    }
}


void
ErrorCheckingExportModule :: writeCheck(Domain *domain, TimeStep *tStep)
{
    if ( tStep->isTheFirstStep() ) {
        std :: cout << "#%BEGIN_CHECK% tolerance 1.e-3\n";
    }

    for ( auto &dman : domain->giveDofManagers() ) {
        for ( Dof *dof: *dman ) {
            if ( dof->giveEqn() < 0 ) {
                continue;
            }
            std :: cout << "#NODE tStep " << tStep->giveNumber();
            std :: cout << " number " << dman->giveNumber();
            std :: cout << " dof " << dof->giveDofID(); 
            std :: cout << " unknown " << 'd';
            std :: cout << " value " << dof->giveUnknown(VM_Total, tStep);
            std :: cout << std :: endl;
        }
    }

    for ( auto &element : domain->giveElements() ) {
        IntegrationRule *iRule = element->giveDefaultIntegrationRulePtr();
        FloatArray ipval;
        for ( int ist: this->writeIST ) {
            for ( int gpnum = 0; gpnum < iRule->giveNumberOfIntegrationPoints(); ++gpnum ){
                GaussPoint *gp = iRule->getIntegrationPoint(gpnum);
                element->giveIPValue(ipval, gp, (InternalStateType)ist, tStep);
                for ( int component = 1; component <= ipval.giveSize(); ++component ) {
                    std :: cout << "#ELEMENT tStep " << tStep->giveNumber();
                    std :: cout << " number " << element->giveNumber();
                    std :: cout << " gp " << gpnum+1;
                    std :: cout << " keyword " << ist; ///@note This writes IST number and not "stresses"
                    std :: cout << " component " << component;
                    std :: cout << " value " << ipval.at(component);
                    std :: cout << std :: endl;
                }
            }
        }
    }

    if ( !tStep->isNotTheLastStep() ) {
        std :: cout << "#%END_CHECK%" << std :: endl;
    }
}


bool
ErrorCheckingExportModule :: scanToErrorChecks(std :: ifstream &stream, double &errorTolerance)
{
    while ( !stream.eof() ) {
        std :: string line;
        std :: getline(stream, line);
        if ( line.compare(0, 14, "#%BEGIN_CHECK%") == 0 ) {
            errorTolerance = 1e-6;
            if ( line.size() >= 25 ) {
                errorTolerance = std :: stod( line.substr(25) );
            }
            return true;
        }
    }
    return false;
}

std::unique_ptr<ErrorCheckingRule>
ErrorCheckingExportModule :: giveErrorCheck(std :: ifstream &stream, double errorTolerance)
{
    std :: string line;
    while ( !stream.eof() ) {
        std :: getline(stream, line);
        if ( line.compare(0, 12, "#%END_CHECK%") == 0 ) {
            return NULL;
        }
        if (line.size() > 2 && line[0] == '#' && line[1] != '#') {
            break;
        }
        
    }

    if ( line.compare(0, 5, "#NODE") == 0 ) {
        return std::make_unique<NodeErrorCheckingRule>(line, errorTolerance);
    } else if ( line.compare(0, 12, "#ELEMENTNODE") == 0 ) {
        return std::make_unique<InternalElementDofManErrorCheckingRule>(line, errorTolerance);
    } else if ( line.compare(0, 8, "#ELEMENT") == 0 ) {
        return std::make_unique<ElementErrorCheckingRule>(line, errorTolerance);
    } else if ( line.compare(0, 13, "#BEAM_ELEMENT") == 0 ) {
        return std::make_unique<BeamElementErrorCheckingRule>(line, errorTolerance);
    } else if ( line.compare(0, 9, "#REACTION") == 0 ) {
        return std::make_unique<ReactionErrorCheckingRule>(line, errorTolerance);
    } else if ( line.compare(0, 10, "#LOADLEVEL") == 0 ) {
        return std::make_unique<LoadLevelErrorCheckingRule>(line, errorTolerance);
    } else if ( line.compare(0, 7, "#EIGVAL") == 0 ) {
        return std::make_unique<EigenValueErrorCheckingRule>(line, errorTolerance);
    } else if ( line.compare(0, 5, "#TIME") == 0 ) {
        return std::make_unique<TimeCheckingRule>(line, errorTolerance);
    } else {
        OOFEM_ERROR("Unsupported rule '%s'", line.c_str());
    }
}

} // end namespace oofem
