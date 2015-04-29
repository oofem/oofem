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

#include "../sm/EngineeringModels/xfemstatic.h"
#include "../sm/xfem/xfemstructuralelementinterface.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "../sm/mappers/primvarmapper.h"
#include "timestep.h"
#include "metastep.h"
#include "dictionary.h"
#include "classfactory.h"
#include "dofmanager.h"
#include "element.h"
#include "exportmodulemanager.h"
#include "vtkxmlexportmodule.h"
#include "gausspoint.h"
#include "matstatmapperint.h"
#include "unknownnumberingscheme.h"

namespace oofem {
REGISTER_EngngModel(XFEMStatic);

XFEMStatic :: XFEMStatic(int i, EngngModel *_master) :
    NonLinearStatic(i, _master),
    updateStructureFlag(false),
    mForceRemap(false),
    mSetValsFromDofMap(true)
{ }

XFEMStatic :: ~XFEMStatic() { }

void
XFEMStatic :: solveYourselfAt(TimeStep *tStep)
{
    // Initiates the total displacement to zero in the UnknownsDictionary at the first time step.
    // this must be done in order to support a dynamic equation system
    if ( tStep->isTheFirstStep() ) {
        printf("Initializing DofUnknownsDictionary... \n");
        this->initializeDofUnknownsDictionary(tStep);
    }

    // Initialization
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ); // 1 stands for domain?
    if ( this->needsStructureUpdate() ) {
        //      printf("Increasing size of displacement array.\n");

        if ( totalDisplacement.giveSize() != neq ) {
            if ( mSetValsFromDofMap ) {
                FloatArray totalDisplacementNew;
                setValsFromDofMap(totalDisplacementNew, totalDisplacement);
                totalDisplacement = totalDisplacementNew;
            } else {
                totalDisplacement.resize(neq);
                totalDisplacement.zero();
            }
        }


        if ( incrementOfDisplacement.giveSize() != neq ) {
            if ( mSetValsFromDofMap ) {
                FloatArray incrementOfDisplacementNew;
                setValsFromDofMap(incrementOfDisplacementNew, incrementOfDisplacement);
                incrementOfDisplacement = incrementOfDisplacementNew;
            } else {
                incrementOfDisplacement.resize(neq);
                incrementOfDisplacement.zero();
            }
        }

        if ( mSetValsFromDofMap ) {
            this->setTotalDisplacementFromUnknownsInDictionary(VM_Total, tStep);
        }

        if ( incrementalLoadVector.giveSize() != neq ) {
            if ( mSetValsFromDofMap ) {
                FloatArray incrementalLoadVectorNew;
                incrementalLoadVector.zero();                 // temp JB - load vector needs to be recomputed if xfem dofs are introduced
                setValsFromDofMap(incrementalLoadVectorNew, incrementalLoadVector);
                incrementalLoadVector = incrementalLoadVectorNew;
            } else {
                incrementalLoadVector.resize(neq);
                incrementalLoadVector.zero();
            }
        }

        ///////////////////////////////////////////////////////////////////
        // Map values in the old initialLoadVector to the new initialLoadVector
        if ( initialLoadVector.giveSize() != neq ) {
            if ( mSetValsFromDofMap ) {
                FloatArray initialLoadVectorNew;
                initialLoadVector.zero();                 // temp JB - load vector needs to be recomputed if xfem dofs are introduced
                setValsFromDofMap(initialLoadVectorNew, initialLoadVector);
                initialLoadVector = initialLoadVectorNew;
            } else {
                initialLoadVector.resize(neq);
                initialLoadVector.zero();
            }
        }
    }

    this->setUpdateStructureFlag(false);
    NonLinearStatic :: solveYourselfAt(tStep);
}

void
XFEMStatic :: terminate(TimeStep *tStep)
{
    this->doStepOutput(tStep);
    this->printReactionForces(tStep, 1);
    // update load vectors before storing context
    fflush( this->giveOutputStream() );
    this->updateLoadVectors(tStep);
    this->saveStepContext(tStep);

    // Propagate fronts
    for ( auto &domain: domainList ) {
        XfemManager *xMan = domain->giveXfemManager();
        bool frontsHavePropagated = false;
        xMan->propagateFronts(frontsHavePropagated);
    }


    // Update element subdivisions if necessary
    // (e.g. if a crack has moved and cut a new element)
    for ( int domInd = 1; domInd <= this->giveNumberOfDomains(); domInd++ ) {
        Domain *domain = this->giveDomain(domInd);

        // create a new set containing all elements
        Set elemSet(0, domain);
        elemSet.addAllElements();

        if ( domain->giveXfemManager()->hasPropagatingFronts() || mForceRemap ) {
            // If domain cloning is performed, there is no need to
            // set values from the dof map.
            mSetValsFromDofMap = false;

            // Take copy of the domain to allow mapping of state variables
            // to the new Gauss points.
            Domain *dNew = domain->Clone();

            bool deallocateOld = false;
            setDomain(1, dNew, deallocateOld);
            forceEquationNumbering();

            // Map primary variables
            LSPrimaryVariableMapper primMapper;
            FloatArray u;
            primMapper.mapPrimaryVariables(u, * domain, * dNew, VM_Total, * tStep);


            if ( totalDisplacement.giveSize() == u.giveSize() ) {
                FloatArray diff;
                diff.beDifferenceOf(totalDisplacement, u);

                printf( "diff norm: %e\n", diff.computeNorm() );
            }

            totalDisplacement = u;


            primMapper.mapPrimaryVariables(incrementOfDisplacement, * domain, * dNew, VM_Incremental, * tStep);


            int numEl = dNew->giveNumberOfElements();

            for ( int i = 1; i <= numEl; i++ ) {
                ////////////////////////////////////////////////////////
                // Map state variables for regular Gauss points
                StructuralElement *el = dynamic_cast< StructuralElement * >( dNew->giveElement(i) );
                el->createMaterialStatus();
                el->mapStateVariables(* domain, * tStep);


                ////////////////////////////////////////////////////////
                // Map state variables for cohesive zone if applicable
                XfemStructuralElementInterface *xFemEl = dynamic_cast< XfemStructuralElementInterface * >(el);
                if ( xFemEl != NULL ) {
                    if ( xFemEl->mpCZMat != NULL ) {
                        size_t numCzRules = xFemEl->mpCZIntegrationRules.size();

                        for ( size_t czIndex = 0; czIndex < numCzRules; czIndex++ ) {
                            if ( xFemEl->mpCZIntegrationRules [ czIndex ] != NULL ) {
                                for ( GaussPoint *gp: *xFemEl->mpCZIntegrationRules [ czIndex ] ) {

                                    MaterialStatus *ms = xFemEl->mpCZMat->giveStatus(gp);
                                    if ( ms == NULL ) {
                                        OOFEM_ERROR("Failed to fetch material status.");
                                    }

                                    MaterialStatusMapperInterface *interface = dynamic_cast< MaterialStatusMapperInterface * >
                                                                               ( xFemEl->mpCZMat->giveStatus(gp) );

                                    if ( interface == NULL ) {
                                        OOFEM_ERROR("Failed to fetch MaterialStatusMapperInterface.");
                                    }


                                    MaterialStatus *matStat = dynamic_cast< MaterialStatus * >( xFemEl->mpCZMat->giveStatus(gp) );
                                    StructuralInterfaceMaterialStatus *siMatStat = dynamic_cast< StructuralInterfaceMaterialStatus * >(matStat);
                                    if ( siMatStat == NULL ) {
                                        OOFEM_ERROR("Failed to cast to StructuralInterfaceMaterialStatus.");
                                    }
                                    interface->MSMI_map_cz(* gp, * domain, elemSet, * tStep, * siMatStat);
                                }
                            }
                        }
                    }
                }
            }


            delete domain;
            domain = this->giveDomain(1);

            // Set domain pointer to various components ...
            this->nMethod->setDomain(domain);

            int numExpModules = this->exportModuleManager->giveNumberOfModules();
            for ( int i = 1; i <= numExpModules; i++ ) {
                //  ... by diving deep into the hierarchies ... :-/
                VTKXMLExportModule *vtkxmlMod = dynamic_cast< VTKXMLExportModule * >( this->exportModuleManager->giveModule(i) );
                if ( vtkxmlMod != NULL ) {
                    vtkxmlMod->giveSmoother()->setDomain(domain);
                    vtkxmlMod->givePrimVarSmoother()->setDomain(domain);
                }
            }


            this->setUpdateStructureFlag(true);
        } // if( domain->giveXfemManager()->hasPropagatingFronts() )

        //#endif
    }

    // Fracture/failure mechanics evaluation
    for ( auto &domain: domainList ) {
        if ( domain->hasFractureManager() ) { // Will most likely fail if numDom > 1
            FractureManager *fracMan = domain->giveFractureManager();
            fracMan->evaluateYourself(tStep);
            fracMan->updateXFEM(tStep); // Update XFEM structure based on the fracture manager

            this->setUpdateStructureFlag( fracMan->giveUpdateFlag() ); // if the internal structure need to be updated
        }
    }
}

void
XFEMStatic :: updateLoadVectors(TimeStep *tStep)
{
    MetaStep *mstep = this->giveMetaStep( tStep->giveMetaStepNumber() );
    bool isLastMetaStep = ( tStep->giveNumber() == mstep->giveLastStepNumber() );

    if ( controlMode == nls_indirectControl ) { //todo@: not checked
        //if ((tStep->giveNumber() == mstep->giveLastStepNumber()) && ir->hasField("fixload")) {
        if ( isLastMetaStep ) {
            if ( !mstep->giveAttributesRecord()->hasField(_IFT_NonLinearStatic_donotfixload) ) {
                OOFEM_LOG_INFO("Fixed load level\n");

                //update initialLoadVector
                if ( initialLoadVector.isEmpty() ) {
                    initialLoadVector.resize( incrementalLoadVector.giveSize() );
                }

                incrementalLoadVector.times(loadLevel);
                initialLoadVector.add(incrementalLoadVector);

                incrementalLoadVectorOfPrescribed.times(loadLevel);
                initialLoadVectorOfPrescribed.add(incrementalLoadVectorOfPrescribed);

                incrementalLoadVector.zero();
                incrementalLoadVectorOfPrescribed.zero();

                this->loadInitFlag = 1;
            }

            //if (!mstep->giveAttributesRecord()->hasField("keepll")) this->loadLevelInitFlag = 1;
        }
    } else { // direct control
        //update initialLoadVector after each step of direct control
        //(here the loading is not proportional)

        /*if ( initialLoadVector.isEmpty() ) {
         *  initialLoadVector.resize( incrementalLoadVector.giveSize() );
         * }
         */
        OOFEM_LOG_DEBUG("Fixed load level\n");

        int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ); // 1 stands for domain?

        //        printf("Before: ");
        //        incrementalLoadVector.printYourself();
        if ( incrementalLoadVector.giveSize() != neq ) {
            //initialLoadVector.resize( incrementalLoadVector.giveSize() );
            //          initialLoadVector.printYourself();
            //            initialLoadVector.resize( 0 );

            FloatArray incrementalLoadVectorNew;
            setValsFromDofMap(incrementalLoadVectorNew, incrementalLoadVector);
            incrementalLoadVector = incrementalLoadVectorNew;
        }
        //        printf("After: ");
        //        incrementalLoadVector.printYourself();

        incrementalLoadVector.times(loadLevel);

        ///////////////////////////////////////////////////////////////////
        // Map values in the old initialLoadVector to the new initialLoadVector
        //        printf("Before: ");
        //        initialLoadVector.printYourself();
        if ( initialLoadVector.giveSize() != neq ) {
            //initialLoadVector.resize( incrementalLoadVector.giveSize() );
            //          initialLoadVector.printYourself();
            //            initialLoadVector.resize( 0 );

            FloatArray initialLoadVectorNew;
            setValsFromDofMap(initialLoadVectorNew, initialLoadVector);
            initialLoadVector = initialLoadVectorNew;
        }
        //        printf("After: ");
        //        initialLoadVector.printYourself();

        ///////////////////////////////////////////////////////////////////

        initialLoadVector.add(incrementalLoadVector);

        incrementalLoadVectorOfPrescribed.times(loadLevel);
        initialLoadVectorOfPrescribed.add(incrementalLoadVectorOfPrescribed);

        incrementalLoadVector.zero();
        incrementalLoadVectorOfPrescribed.zero();

        this->loadInitFlag = 1;
    }


    // if (isLastMetaStep) {
    if ( isLastMetaStep && !mstep->giveAttributesRecord()->hasField(_IFT_NonLinearStatic_donotfixload) ) {
#ifdef VERBOSE
        OOFEM_LOG_INFO("Reseting load level\n");
#endif
        if ( mstepCumulateLoadLevelFlag ) {
            cumulatedLoadLevel += loadLevel;
        } else {
            cumulatedLoadLevel = 0.0;
        }

        this->loadLevel = 0.0;
    }
}

void
XFEMStatic :: updateYourself(TimeStep *tStep)
{
    NonLinearStatic :: updateYourself(tStep);

    // TODO: Check if update is actually needed
    this->setUpdateStructureFlag(true);

    // Update the UnknownsDictionary if needed
    if ( this->needsStructureUpdate() ) {
        printf(" Updating DofUnknownsDictionary... \n");

        buildDofMap();

        for ( auto &domain: domainList ) {
            int nnodes = domain->giveNumberOfDofManagers();
            for ( int inode = 1; inode <= nnodes; inode++ ) {
                this->updateDofUnknownsDictionary(domain->giveDofManager(inode), tStep);
            }
        }
    }
}

double
XFEMStatic ::  giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    // Returns the unknown quantity corresponding to the dof
    if ( this->requiresUnknownsDictionaryUpdate() ) {
        int hash = this->giveUnknownDictHashIndx(mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else { // Value is not initiated in UnknownsDictionary
            return 0.0; ///@todo: how should one treat newly created dofs?
            // If we are not happy with setting them to zero,
            // I suggest that we use the primary variable mapper. /ES
            //OOFEM_ERROR("Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    } else {
        return NonLinearStatic ::  giveUnknownComponent(mode, tStep, d, dof);
    }
}

IRResultType XFEMStatic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = NonLinearStatic :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    int remapFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, remapFlag, _IFT_XFEMStatic_ForceRemap);

    if ( remapFlag == 1 ) {
        printf("Forcing remapping.\n");
        mForceRemap = true;
    }

    return IRRT_OK;
}

void
XFEMStatic :: initializeDofUnknownsDictionary(TimeStep *tStep)
{
    // Initializes all dof values to zero
    for ( auto &domain: domainList ) {
        int nnodes = domain->giveNumberOfDofManagers();
        for ( int inode = 1; inode <= nnodes; inode++ ) {
            DofManager *node = domain->giveDofManager(inode);
            for ( Dof *dof: *node ) {
                dof->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, 0.0);
            }
        }
    }
}

void
XFEMStatic :: setTotalDisplacementFromUnknownsInDictionary(ValueModeType mode, TimeStep *tStep)
{
    printf("Entering XFEMStatic :: setTotalDisplacementFromUnknownsInDictionary().\n");

    // Sets the values in the displacement vector based on stored values in the unknowns dictionaries.
    // Used in the beginning of each time step.
    for ( auto &domain: domainList ) {
        for ( int j = 1; j <= domain->giveNumberOfDofManagers(); j++ ) {
            DofManager *inode = domain->giveDofManager(j);
            int eqNum;
            for ( Dof *dof: *inode ) {
                eqNum = dof->giveEqn();
                if ( eqNum > 0 ) {
                    double val = dof->giveUnknown(mode, tStep);
                    totalDisplacement.at(eqNum) = val;
                }
            }
        }
    }
}

void
XFEMStatic :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{
    // update DoF unknowns dictionary.
    double val;
    for ( Dof *dof: *inode ) {
        int eqNum = dof->__giveEquationNumber();
        if ( dof->hasBc(tStep) ) {
            val = dof->giveBcValue(VM_Total, tStep);
        } else {
            if ( eqNum > 0 ) {
                val = totalDisplacement.at(eqNum);
            } else { // new eq number
                val = 0.0;
            }
        }

        dof->updateUnknownsDictionary(tStep, VM_Total, val);
    }
}

void XFEMStatic :: buildDofMap()
{
    printf("Building dof map.\n");
    mDofEqnNumMap.clear();

    for ( int domainIndex = 1; domainIndex <= this->giveNumberOfDomains(); domainIndex++ ) {
        Domain *domain = this->giveDomain(domainIndex);

        for ( int dManIndex = 1; dManIndex <= domain->giveNumberOfDofManagers(); dManIndex++ ) {
            DofManager *dMan = domain->giveDofManager(dManIndex);

            for ( Dof *dof: *dMan ) {
                int eqNum = dof->giveEqn();

                if ( eqNum > 0 ) {
                    std :: vector< int > key(3);
                    key [ 0 ] = domainIndex;
                    key [ 1 ] = dManIndex;
                    key [ 2 ] = dof->giveDofID();

                    mDofEqnNumMap [ key ] = eqNum;
                }
            }
        }
    }
}

void XFEMStatic :: setValsFromDofMap(FloatArray &oArray, const FloatArray &iArray)
{
    int neq = 0;
    for ( int domainIndex = 1; domainIndex <= this->giveNumberOfDomains(); domainIndex++ ) {
        neq += this->giveNumberOfDomainEquations( domainIndex, EModelDefaultEquationNumbering() );
    }

    int numEqOld = iArray.giveSize();
    printf("Setting values from dof map. neq: %d numEqOld: %d\n", neq, numEqOld);


    oArray.resize(neq);
    oArray.zero();

    for ( int domainIndex = 1; domainIndex <= this->giveNumberOfDomains(); domainIndex++ ) {
        Domain *domain = this->giveDomain(domainIndex);

        for ( int dManIndex = 1; dManIndex <= domain->giveNumberOfDofManagers(); dManIndex++ ) {
            DofManager *dMan = domain->giveDofManager(dManIndex);

            for ( Dof *dof: *dMan ) {
                int eqNumNew = dof->giveEqn();

                if ( eqNumNew > 0 ) {
                    std :: vector< int > key(3);
                    key [ 0 ] = domainIndex;
                    key [ 1 ] = dManIndex;
                    key [ 2 ] = dof->giveDofID();

                    if ( mDofEqnNumMap.find(key) != mDofEqnNumMap.end() ) {
                        int eqNumOld = mDofEqnNumMap [ key ];
                        //                      printf("eqNumNew: %d eqNumOld: %d\n", eqNumNew, eqNumOld);

                        if ( eqNumOld > 0 && eqNumOld <= numEqOld  ) {
                            oArray.at(eqNumNew) = iArray.at(eqNumOld);
                        }
                    }
                }
            }
        }
    }
}
} /* namespace oofem */
