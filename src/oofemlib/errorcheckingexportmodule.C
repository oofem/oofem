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
#include "irresulttype.h"
#ifdef __SM_MODULE
 #include "../sm/EngineeringModels/structengngmodel.h"
 #include "../sm/Elements/Beams/beam2d.h"
 #include "../sm/Elements/Beams/beam3d.h"
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
  IRResultType result;
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

bool
NodeErrorCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    // Rule doesn't apply yet.
    if ( (tStep->giveNumber() != tstep ) || (tStep->giveVersion() != tsubstep) )  {
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
        OOFEM_WARNING("Check failed in: tstep %d, node %d, dof %d, mode %d:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      tstep, number, dofid, mode,
                      dmanValue, value, fabs(dmanValue-value), tolerance );
    }
    return check;
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

bool
ElementErrorCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    // Rule doesn't apply yet.
     if ( (tStep->giveNumber() != tstep ) || (tStep->giveVersion() != tsubstep) )  {
        return true;
    }

    FloatArray ipval;
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

    // note! GPs are numbered from 0 internally, but written with 1-index, inconsistent!
    GaussPoint *gp = element->giveIntegrationRule(irule)->getIntegrationPoint(gpnum-1);
    element->giveIPValue(ipval, gp, ist, tStep);

    if ( component > ipval.giveSize() || component < 1 ) {
        OOFEM_WARNING("Check failed in: element %d, gpnum %d, ist %d, component %d:\n"
                      "Component not found!",
                      number, gpnum, ist, component);
        ipval.printYourself();
        return false;
    }

    double elementValue = ipval.at(component);
    bool check = checkValue(elementValue);
    if ( !check ) {
        OOFEM_WARNING("Check failed in: tstep %d, element %d, gpnum %d, ist %d, component %d:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      tstep, number, gpnum, ist, component,
                      elementValue, value, fabs(elementValue-value), tolerance );
        ipval.printYourself();
    }
    return check;
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

bool
BeamElementErrorCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    // Rule doesn't apply yet.
    if ( (tStep->giveNumber() != tstep ) || (tStep->giveVersion() != tsubstep) )  {
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
      if(Beam2d* b = dynamic_cast<Beam2d*>(element)) b->giveEndForcesVector(val, tStep);
      else if(Beam3d* b = dynamic_cast<Beam3d*>(element)) b->giveEndForcesVector(val, tStep);
      else {
        OOFEM_WARNING("Element %d has no beam interface.", number);
        return false;
      }
    }

    if ( component > val.giveSize() || component < 1 ) {
        OOFEM_WARNING("Check failed in: beam_element %d, ist %d, component %d:\n"
                      "Component not found!",
                      number, ist, component);
        val.printYourself();
        return false;
    }

    double elementValue = val.at(component);
    bool check = checkValue(elementValue);
    if ( !check ) {
        OOFEM_WARNING("Check failed in: tstep %d, beam_element %d, ist %d, component %d:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      tstep, number, ist, component,
                      elementValue, value, fabs(elementValue-value), tolerance );
        val.printYourself();
    }
    return check;
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

bool
ReactionErrorCheckingRule :: check(Domain *domain, TimeStep *tStep)
{
    // Rule doesn't apply yet.
    if ( (tStep->giveNumber() != tstep ) || (tStep->giveVersion() != tsubstep) )  {
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
        OOFEM_WARNING("Check failed in: tstep %d, reaction forces number %d, dof %d:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      tstep, number, dofid,
                      reactionForce, value, fabs(reactionForce-value), tolerance );
    }
    return check;
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
        OOFEM_WARNING("Check failed in: tstep %d, load level:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      tstep,
                      loadLevel, value, fabs(loadLevel-value), tolerance );
    }
    return check;
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
        OOFEM_WARNING("Check failed in: tstep %d, eigen value %d:\n"
                      "value is %.8e, but should be %.8e ( error is %e but tolerance is %e )",
                      tstep, number,
                      eig, value, fabs(eig-value), tolerance );
    }
    return check;
}


////////////////////////////////////////////////////////////////////////////////////////////////////


ErrorCheckingExportModule :: ErrorCheckingExportModule(int n, EngngModel *e) : ExportModule(n, e)
{
    allPassed = true;
    writeChecks = false;
}

IRResultType
ErrorCheckingExportModule :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    allPassed = true;
    this->errorCheckingRules.clear();

    filename = std::string("");

    if ( ir->hasField(_IFT_ErrorCheckingExportModule_filename) ) {
        IR_GIVE_FIELD(ir, this->filename, _IFT_ErrorCheckingExportModule_filename);
    }
    else {
        filename = emodel->giveReferenceFileName();
    }

    // Reads all the rules;
    std :: ifstream inputStream(this->filename);
    if ( !inputStream ) {
        OOFEM_WARNING("Couldn't open file '%s'\n", this->filename.c_str());
        return IRRT_BAD_FORMAT;
    }
    double tol = 0.;
    if ( this->scanToErrorChecks(inputStream,  tol) ) {
        for (;;) {
            std :: unique_ptr< ErrorCheckingRule > rule(this->giveErrorCheck(inputStream, tol));
            if ( !rule ) {
                break;
            }
            errorCheckingRules.push_back(std :: move(rule));
        }
    }

    this->writeIST.clear();
    writeChecks = ir->hasField(_IFT_ErrorCheckingExportModule_writeIST);
    if ( writeChecks ) {
        IR_GIVE_FIELD(ir, this->writeIST, _IFT_ErrorCheckingExportModule_writeIST);
    }

    if ( errorCheckingRules.size() == 0 && !writeChecks ) {
        OOFEM_WARNING("No rules found (possibly wrong file or syntax).");
    }

    return ExportModule :: initializeFrom(ir);
}

void
ErrorCheckingExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
#if 0
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }
#endif

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

ErrorCheckingRule *
ErrorCheckingExportModule :: giveErrorCheck(std :: ifstream &stream, double errorTolerance)
{
    std :: string line;
    while ( !stream.eof() ) {
        std :: getline(stream, line);
        if ( line.compare(0, 12, "#%END_CHECK%") == 0 ) {
            return NULL;
        }
        if ( line.size() > 2 && line[0] == '#' && line[1] != '#' ) {
            break;
        }
    }

    if ( line.compare(0, 5, "#NODE") == 0 ) {
        return new NodeErrorCheckingRule(line, errorTolerance);
    } else if ( line.compare(0, 8, "#ELEMENT") == 0 ) {
        return new ElementErrorCheckingRule(line, errorTolerance);
    } else if ( line.compare(0, 13, "#BEAM_ELEMENT") == 0 ) {
        return new BeamElementErrorCheckingRule(line, errorTolerance);
    } else if ( line.compare(0, 9, "#REACTION") == 0 ) {
        return new ReactionErrorCheckingRule(line, errorTolerance);
    } else if ( line.compare(0, 10, "#LOADLEVEL") == 0 ) {
        return new LoadLevelErrorCheckingRule(line, errorTolerance);
    } else if ( line.compare(0, 7, "#EIGVAL") == 0 ) {
        return new EigenValueErrorCheckingRule(line, errorTolerance);
    } else {
        OOFEM_ERROR("Unsupported rule '%s'", line.c_str());
        return NULL;
    }
}

} // end namespace oofem
