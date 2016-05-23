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

#include "../sm/EngineeringModels/xfemsolverinterface.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "../sm/xfem/xfemstructuralelementinterface.h"
#include "../sm/mappers/primvarmapper.h"
#include "timestep.h"
#include "structengngmodel.h"
#include "staticstructural.h"
#include "primaryfield.h"
#include "domain.h"
#include "xfem/xfemmanager.h"
#include "element.h"
#include "matstatmapperint.h"
#include "nummet.h"
#include "floatarray.h"
#include "exportmodulemanager.h"
#include "vtkxmlexportmodule.h"

namespace oofem {

XfemSolverInterface::XfemSolverInterface():
mNeedsVariableMapping(false)
{

}

XfemSolverInterface::~XfemSolverInterface()
{

}

void XfemSolverInterface::propagateXfemInterfaces(TimeStep *tStep, StructuralEngngModel &ioEngngModel, bool iRecomputeStepAfterCrackProp)
{
    int domainInd = 1;
    Domain *domain = ioEngngModel.giveDomain(domainInd);

    if(domain->hasXfemManager()) {
        XfemManager *xMan = domain->giveXfemManager();
        bool frontsHavePropagated = false;
        if ( xMan->hasInitiationCriteria() ) {
            // TODO: generalise this?
            // Intitiate delaminations (only implemented for listbasedEI/delamination. Treated the same way as propagation)
            xMan->initiateFronts(frontsHavePropagated,tStep);
        }

        if( xMan->hasPropagatingFronts() ) {
            // Propagate crack tips
            xMan->propagateFronts(frontsHavePropagated);

        }

        bool eiWereNucleated = false;
        if( xMan->hasNucleationCriteria() ) {
        	xMan->nucleateEnrichmentItems(eiWereNucleated);
       }

        if(frontsHavePropagated || eiWereNucleated) {
            mNeedsVariableMapping = false;
            
            if ( mNeedsVariableMapping ) {
                mapVariables(tStep, ioEngngModel);
            } else {
                ioEngngModel.giveDomain(1)->postInitialize();
                ioEngngModel.forceEquationNumbering();
            }

            if(iRecomputeStepAfterCrackProp) {
                printf("Recomputing time step.\n");
                ioEngngModel.forceEquationNumbering();
                ioEngngModel.solveYourselfAt(tStep);
                ioEngngModel.updateYourself( tStep );
                ioEngngModel.terminate( tStep );

            }
        }
    }

}

void XfemSolverInterface::mapVariables(TimeStep *tStep, StructuralEngngModel &ioEngngModel)
{
    if(mNeedsVariableMapping) {

        int domainInd = 1;
        Domain *domain = ioEngngModel.giveDomain(domainInd);


        // create a new set containing all elements
        Set elemSet(0, domain);
        elemSet.addAllElements();

        // Clear stiffness matrix
//        int di = 1;
//        int neq = this->giveNumberOfDomainEquations( di, EModelDefaultEquationNumbering() );
//        printf("neq: %d\n", neq);
//
//        delete stiffnessMatrix;
//        stiffnessMatrix = NULL;


        // Create new domain and map variables
        Domain *dNew = domain->Clone();
        bool deallocateOld = false;
        ioEngngModel.setDomain(1, dNew, deallocateOld);

        ioEngngModel.forceEquationNumbering();




#if 1
        // Map primary variables ...
        LSPrimaryVariableMapper primMapper;
        FloatArray u;
        primMapper.mapPrimaryVariables(u, * domain, * dNew, VM_Total, * tStep);

        xfemUpdatePrimaryField(ioEngngModel, tStep, u);
        
        // ... and write to PrimaryField
//        field->update(VM_Total, tStep, u);

//        PrimaryField *oldField = field;
//        field = new PrimaryField(this, 1, FT_Displacements, 1);
//        delete oldField;


        // Map state variables
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

                                if(!siMatStat->giveNewlyInserted()) {
                                    interface->MSMI_map_cz(* gp, * domain, elemSet, * tStep, * siMatStat);
                                }
                            }
                        }
                    }
                }
            }
        }
#endif



        delete domain;
        domain = ioEngngModel.giveDomain(1);

        // Set domain pointer to various components ...
        ioEngngModel.giveNumericalMethod(NULL)->setDomain(domain);
//        ioEngngModel.nMethod->setDomain(domain);

        int numExpModules = ioEngngModel.giveExportModuleManager()->giveNumberOfModules();
        for ( int i = 1; i <= numExpModules; i++ ) {
            //  ... by diving deep into the hierarchies ... :-/
            VTKXMLExportModule *vtkxmlMod = dynamic_cast< VTKXMLExportModule * >( ioEngngModel.giveExportModuleManager()->giveModule(i) );
            if ( vtkxmlMod != NULL ) {
                vtkxmlMod->giveSmoother()->setDomain(domain);
                vtkxmlMod->givePrimVarSmoother()->setDomain(domain);
            }
        }



        mNeedsVariableMapping = false;
    }
}

#if 0
void XfemSolverInterface::mapVariables(TimeStep *tStep, StaticStructural &ioEngngModel)
{
    printf("Entering femSolverInterface::mapVariables(TimeStep *tStep, StaticStructural &ioEngngModel).\n");
    if(mNeedsVariableMapping) {
        int domainInd = 1;
        Domain *domain = ioEngngModel.giveDomain(domainInd);


        // create a new set containing all elements
        Set elemSet(0, domain);
        elemSet.addAllElements();

        // Clear stiffness matrix
//        int di = 1;
//        int neq = this->giveNumberOfDomainEquations( di, EModelDefaultEquationNumbering() );
//        printf("neq: %d\n", neq);
//
//        delete stiffnessMatrix;
//        stiffnessMatrix = NULL;


        // Create new domain and map variables
        Domain *dNew = domain->Clone();
        bool deallocateOld = false;
        ioEngngModel.setDomain(1, dNew, deallocateOld);

        ioEngngModel.forceEquationNumbering();




        // Map primary variables ...
        LSPrimaryVariableMapper primMapper;
        FloatArray u;
        primMapper.mapPrimaryVariables(u, * domain, * dNew, VM_Total, * tStep);

        // ... and write to PrimaryField
        printf("ioEngngModel.updatePrimaryField.\n");
        ioEngngModel.updatePrimaryField(VM_Total, tStep, u);

//        PrimaryField *oldField = field;
//        field = new PrimaryField(this, 1, FT_Displacements, 1);
//        delete oldField;


        // Map state variables
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
        domain = ioEngngModel.giveDomain(1);

        // Set domain pointer to various components ...
        ioEngngModel.giveNumericalMethod(NULL)->setDomain(domain);
//        ioEngngModel.nMethod->setDomain(domain);

        int numExpModules = ioEngngModel.giveExportModuleManager()->giveNumberOfModules();
        for ( int i = 1; i <= numExpModules; i++ ) {
            //  ... by diving deep into the hierarchies ... :-/
            VTKXMLExportModule *vtkxmlMod = dynamic_cast< VTKXMLExportModule * >( ioEngngModel.giveExportModuleManager()->giveModule(i) );
            if ( vtkxmlMod != NULL ) {
                vtkxmlMod->giveSmoother()->setDomain(domain);
                vtkxmlMod->givePrimVarSmoother()->setDomain(domain);
            }
        }



        mNeedsVariableMapping = false;

    }
}
#endif

void XfemSolverInterface::xfemUpdatePrimaryField(StaticStructural &ioEngngModel, TimeStep *tStep, const FloatArray &iNewSol)
{
    printf("XfemSolverInterface::updatePrimaryField(StaticStructural &ioEngngModel, TimeStep *tStep, const FloatArray &iNewSol)\n");
    ioEngngModel.setSolution(tStep, iNewSol);
}

} /* namespace oofem */
