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

#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iterator>

#include "matlabexportmodule.h"
#include "engngm.h"
#include "node.h"
#include "mathfem.h"
#include "gausspoint.h"
#include "weakperiodicbc.h"
#include "solutionbasedshapefunction.h"
#include "timestep.h"
#include "classfactory.h"
#include "unknownnumberingscheme.h"

#ifdef __FM_MODULE
 #include "../fm/tr21stokes.h"
 #include "../fm/tet21stokes.h"
 #include "../fm/stokesflow.h"
#endif

#ifdef __SM_MODULE
 #include "../sm/structengngmodel.h"
#endif


namespace oofem {
REGISTER_ExportModule(MatlabExportModule)

MatlabExportModule :: MatlabExportModule(int n, EngngModel *e) : ExportModule(n, e), internalVarsToExport(), primaryVarsToExport()
{
    exportMesh = false;
    exportData = false;
    exportArea = false;
    exportSpecials = false;
    exportReactionForces = false;
    reactionForcesDofManList.clear();
    exportIntegrationPointFields = false;
    elList.clear();
}


MatlabExportModule :: ~MatlabExportModule()
{ }


IRResultType
MatlabExportModule :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    ExportModule :: initializeFrom(ir);

    exportMesh = ir->hasField(_IFT_MatlabExportModule_mesh);
    exportData = ir->hasField(_IFT_MatlabExportModule_data);
    exportArea = ir->hasField(_IFT_MatlabExportModule_area);
    exportSpecials = ir->hasField(_IFT_MatlabExportModule_specials);

    exportReactionForces = ir->hasField(_IFT_MatlabExportModule_ReactionForces);
    if ( exportReactionForces ) {
        IR_GIVE_OPTIONAL_FIELD(ir, reactionForcesDofManList, _IFT_MatlabExportModule_DofManList);
    }

    exportIntegrationPointFields = ir->hasField(_IFT_MatlabExportModule_IntegrationPoints);
    if ( exportIntegrationPointFields ) {
        IR_GIVE_FIELD(ir, internalVarsToExport, _IFT_MatlabExportModule_internalVarsToExport);
        IR_GIVE_OPTIONAL_FIELD(ir, elList, _IFT_MatlabExportModule_ElementList);
    }

    return IRRT_OK;
}


void
MatlabExportModule :: computeArea()
{
    Domain *domain = emodel->giveDomain(1);

    smax.clear();
    smin.clear();

    for ( int i = 1; i <= domain->giveNumberOfSpatialDimensions(); i++ ) {
        smax.push_back( domain->giveDofManager(1)->giveCoordinate(i) );
        smin.push_back( domain->giveDofManager(1)->giveCoordinate(i) );
    }

    for ( int i = 0; i < domain->giveNumberOfDofManagers(); i++ ) {
        for ( int j = 0; j < domain->giveNumberOfSpatialDimensions(); j++ ) {
            smax.at(j) = max( smax.at(j), domain->giveDofManager(i + 1)->giveCoordinate(j + 1) );
            smin.at(j) = min( smin.at(j), domain->giveDofManager(i + 1)->giveCoordinate(j + 1) );
        }
    }


    Area = 0;
    Volume = 0;

    if ( domain->giveNumberOfSpatialDimensions() == 2 ) {
        for ( int i = 1; i <= domain->giveNumberOfElements(); i++ ) {
            Area = Area + domain->giveElement(i)->computeArea();
        }
    } else {
        for ( int i = 1; i <= domain->giveNumberOfElements(); i++ ) {
            Volume = Volume + domain->giveElement(i)->computeVolume();
        }
    }
}


void
MatlabExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    FILE *FID;
    FID = giveOutputStream(tStep);
    Domain *domain  = emodel->giveDomain(1);
    ndim = domain->giveNumberOfSpatialDimensions();

    // Output header
    fprintf(FID, "%%%% OOFEM generated export file \n");
    fprintf( FID, "%% Output for time %f\n", tStep->giveTargetTime() );


    fprintf( FID, "function [mesh area data specials ReactionForces IntegrationPointFields]=%s\n\n", functionname.c_str() );

    if ( exportMesh ) {
        doOutputMesh(tStep, FID);
    } else {
        fprintf(FID, "\tmesh=[];\n");
    }

    if ( exportData ) {
        doOutputData(tStep, FID);
    } else {
        fprintf(FID, "\tdata=[];\n");
    }

    if ( exportArea ) {
        computeArea();
        fprintf( FID, "\tarea.xmax=%f;\n", smax.at(0) );
        fprintf( FID, "\tarea.xmin=%f;\n", smin.at(0) );
        fprintf( FID, "\tarea.ymax=%f;\n", smax.at(1) );
        fprintf( FID, "\tarea.ymin=%f;\n", smin.at(1) );
        if ( ndim == 2 ) {
            fprintf(FID, "\tarea.area=%f;\n", Area);
            fprintf(FID, "\tvolume=[];\n");
        } else {
            fprintf( FID, "\tarea.zmax=%f;\n", smax.at(2) );
            fprintf( FID, "\tarea.zmin=%f;\n", smin.at(2) );
            fprintf(FID, "\tarea.area=[];\n");
            fprintf(FID, "\tarea.volume=%f;\n", Volume);
        }
    } else {
        fprintf(FID, "\tarea.area=[];\n");
        fprintf(FID, "\tarea.volume=[];\n");
    }

    if ( exportSpecials ) {
        if ( !exportArea ) {
            computeArea();
        }

        doOutputSpecials(tStep, FID);
    } else {
        fprintf(FID, "\tspecials=[];\n");
    }

    // Reaction forces
    if ( exportReactionForces ) {
        doOutputReactionForces(tStep, FID);
    } else {
        fprintf(FID, "\tReactionForces=[];\n");
    }

    // Internal variables in integration points
    if ( exportIntegrationPointFields ) {
        doOutputIntegrationPointFields(tStep, FID);
    } else {
        fprintf(FID, "\tIntegrationPointFields=[];\n");
    }

    fprintf(FID, "\nend\n");
    fclose(FID);
}


void
MatlabExportModule :: doOutputMesh(TimeStep *tStep, FILE *FID)
{
    Domain *domain  = emodel->giveDomain(1);

    fprintf(FID, "\tmesh.p=[");
    for ( int i = 1; i <= domain->giveNumberOfDofManagers(); i++ ) {
        for ( int j = 1; j <= domain->giveNumberOfSpatialDimensions(); j++ ) {
            //double x = domain->giveDofManager(i)->giveCoordinate(1), y = domain->giveDofManager(i)->giveCoordinate(2);
            double c = domain->giveDofManager(i)->giveCoordinate(j);
            fprintf(FID, "%f, ", c);
        }
        fprintf(FID, "; ");
    }

    fprintf(FID, "]';\n");

    int numberOfDofMans = domain->giveElement(1)->giveNumberOfDofManagers();

    fprintf(FID, "\tmesh.t=[");
    for ( int i = 1; i <= domain->giveNumberOfElements(); i++ ) {
        if ( domain->giveElement(i)->giveNumberOfDofManagers() == numberOfDofMans ) {
            for ( int j = 1; j <= domain->giveElement(i)->giveNumberOfDofManagers(); j++ ) {
                fprintf( FID, "%d,", domain->giveElement(i)->giveDofManagerNumber(j) );
            }
        }
        fprintf(FID, ";");
    }

    fprintf(FID, "]';\n");
}


void
MatlabExportModule :: doOutputData(TimeStep *tStep, FILE *FID)
{
    Domain *domain  = emodel->giveDomain(1);
    std :: vector< int >DofIDList;
    std :: vector< int > :: iterator it;
    std :: vector< std :: vector< double > * >valuesList;
    std :: vector< double > *values;

    for ( int i = 1; i <= domain->giveNumberOfDofManagers(); i++ ) {
        for ( int j = 1; j <= domain->giveDofManager(i)->giveNumberOfDofs(); j++ ) {
            Dof *thisDof;
            thisDof = domain->giveDofManager(i)->giveDof(j);
            it = std :: find( DofIDList.begin(), DofIDList.end(), thisDof->giveDofID() );

            if ( it == DofIDList.end() ) {
                DofIDList.push_back( thisDof->giveDofID() );
                values = new(std :: vector< double > );
                valuesList.push_back(values);
            } else {
                std :: size_t pos = it - DofIDList.begin();
                values = valuesList.at(pos);
            }

            double value = thisDof->giveUnknown(VM_Total, tStep);
            values->push_back(value);
        }
    }

    fprintf(FID, "\tdata.DofIDs=[");
    for ( size_t i = 0; i < DofIDList.size(); i++ ) {
        fprintf( FID, "%d, ", DofIDList.at(i) );
    }

    fprintf(FID, "];\n");

    for ( size_t i = 0; i < valuesList.size(); i++ ) {
        fprintf(FID, "\tdata.a{%lu}=[", static_cast< long unsigned int >(i) + 1);
        for ( size_t j = 0; j < valuesList.at(i)->size(); j++ ) {
            fprintf( FID, "%f,", valuesList.at(i)->at(j) );
        }

        fprintf(FID, "];\n");
    }
}


void
MatlabExportModule :: doOutputSpecials(TimeStep *tStep,    FILE *FID)
{
    FloatMatrix v_hat, GradPTemp, v_hatTemp;

    Domain *domain  = emodel->giveDomain(1);

    v_hat.resize(ndim, 1);
    v_hat.zero();

    for ( int i = 1; i <= domain->giveNumberOfElements(); i++ ) {
#ifdef __FM_MODULE

        if ( Tr21Stokes * T = dynamic_cast< Tr21Stokes * >( domain->giveElement(i) ) ) {
            T->giveIntegratedVelocity(v_hatTemp, tStep);
            v_hat.add(v_hatTemp);
        } else if ( Tet21Stokes * T = dynamic_cast< Tet21Stokes * >( domain->giveElement(i) ) ) {
            T->giveIntegratedVelocity(v_hatTemp, tStep);
            v_hat.add(v_hatTemp);
        }

#endif
    }

    // Compute intrinsic area/volume
    double intrinsicSize = 1.0;

    std :: vector< double >V;

    for ( int i = 0; i < ndim; i++ ) {
        intrinsicSize = intrinsicSize * ( smax.at(i) - smin.at(i) );
    }

    for ( int i = 1; i <= ndim; i++ ) {
        V.push_back(v_hat.at(i, 1) / intrinsicSize);
    }

    fprintf(FID, "\tspecials.velocitymean=[");
    for ( int i = 0; i < ndim; i++ ) {
        fprintf( FID, "%e", V.at(i) );
        if ( i != ( ndim - 1 ) ) {
            fprintf(FID, ", ");
        }
    }
    fprintf(FID, "];\n");

    // Output weak periodic boundary conditions
    unsigned int wpbccount = 1, sbsfcount = 1;

    for ( int i = 1; i <= domain->giveNumberOfBoundaryConditions(); i++ ) {
        WeakPeriodicBoundaryCondition *wpbc = dynamic_cast< WeakPeriodicBoundaryCondition * >( domain->giveBc(i) );
        if ( wpbc ) {
            for ( int j = 1; j <= wpbc->giveNumberOfInternalDofManagers(); j++ ) {
                fprintf( FID, "\tspecials.weakperiodic{%u}.descType=%u;\n", wpbccount, wpbc->giveBasisType() );
                fprintf(FID, "\tspecials.weakperiodic{%u}.coefficients=[", wpbccount);
                for ( int k = 1; k <= wpbc->giveInternalDofManager(j)->giveNumberOfDofs(); k++ ) {
                    FloatArray unknowns;
                    IntArray DofMask;
                    double X = wpbc->giveInternalDofManager(j)->giveDof(k)->giveUnknown(VM_Total, tStep);
                    fprintf(FID, "%e\t", X);
                }

                fprintf(FID, "];\n");
                wpbccount++;
            }
        }
        SolutionbasedShapeFunction *sbsf = dynamic_cast< SolutionbasedShapeFunction * >( domain->giveBc(i) );
        if ( sbsf ) {
            fprintf(FID, "\tspecials.solutionbasedsf{%u}.values=[", sbsfcount);
            for ( int k = 1; k <= sbsf->giveInternalDofManager(1)->giveNumberOfDofs(); k++ ) {                  // Only one internal dof manager
                FloatArray unknowns;
                IntArray DofMask;
                double X = sbsf->giveInternalDofManager(1)->giveDof(k)->giveUnknown(VM_Total, tStep);
                fprintf(FID, "%e\t", X);
            }
            fprintf(FID, "];\n");
            sbsfcount++;
        }
    }
}


void
MatlabExportModule :: doOutputReactionForces(TimeStep *tStep,    FILE *FID)
{
    int domainIndex = 1;
    Domain *domain  = emodel->giveDomain(domainIndex);

    FloatArray reactions;
    IntArray dofManMap, dofMap, eqnMap;
#ifdef __SM_MODULE
    StructuralEngngModel *strEngMod = dynamic_cast< StructuralEngngModel * >(emodel);
    if ( strEngMod ) {
        strEngMod->buildReactionTable(dofManMap, dofMap, eqnMap, tStep, domainIndex);
        strEngMod->computeReaction(reactions, tStep, 1);
    } else
#endif
    {
        OOFEM_ERROR("Cannot export reaction forces - only implemented for structural problems.");
    }

    int numDofManToExport = this->reactionForcesDofManList.giveSize();
    if ( numDofManToExport == 0 ) { // No dofMan's given - export every dMan with reaction forces
        for ( int i = 1; i <= domain->giveNumberOfDofManagers(); i++ ) {
            if ( dofManMap.contains(i) ) {
                this->reactionForcesDofManList.followedBy(i);
            }
        }
        numDofManToExport = this->reactionForcesDofManList.giveSize();
    }


    // Output header
    fprintf(FID, "\n %%%% Export of reaction forces \n\n");

    // Output the dofMan numbers that are exported
    fprintf(FID, "\tReactionForces.DofManNumbers = [");
    for ( int i = 1; i <= numDofManToExport; i++ ) {
        fprintf( FID, "%i ", this->reactionForcesDofManList.at(i) );
    }
    fprintf(FID, "];\n");


    // Define the reaction forces as a cell object
    fprintf(FID, "\tReactionForces.ReactionForces = cell(%i,1); \n", numDofManToExport);
    fprintf(FID, "\tReactionForces.DofIDs = cell(%i,1); \n", numDofManToExport);


    // Output the reaction forces for each dofMan. If a certain dof is not prescribed zero is exported.
    IntArray dofIDs;
    for ( int i = 1; i <= numDofManToExport; i++ ) {
        int dManNum = this->reactionForcesDofManList.at(i);

        fprintf(FID, "\tReactionForces.ReactionForces{%i} = [", i);
        if ( dofManMap.contains(dManNum) ) {
            DofManager *dofMan = domain->giveDofManager(dManNum);
            dofIDs.resize( dofMan->giveNumberOfDofs() );
            dofIDs.zero();

            for ( int j = 1; j <= dofMan->giveNumberOfDofs(); j++ ) {
                Dof *dof = dofMan->giveDof(j);
                int num = dof->giveEquationNumber( EModelDefaultPrescribedEquationNumbering() );
                int pos = eqnMap.findFirstIndexOf(num);
                dofIDs.at(j) = ( int ) dof->giveDofID();
                if ( pos > 0 ) {
                    fprintf( FID, "%e ", reactions.at(pos) );
                } else {
                    fprintf(FID, "%e ", 0.0);   // if not prescibed output zero
                }
            }
        }
        fprintf(FID, "];\n");

        // Output dof ID's

        fprintf(FID, "\tReactionForces.DofIDs{%i} = [", i);
        if ( dofManMap.contains(dManNum) ) {
            for ( int j = 1; j <= dofIDs.giveSize(); j++ ) {
                fprintf( FID, "%i ", dofIDs.at(j) );
            }
        }
        fprintf(FID, "];\n");
    }
}


void
MatlabExportModule :: doOutputIntegrationPointFields(TimeStep *tStep,    FILE *FID)
{
    //if ( !testTimeStepOutput(tStep) ) {
    //    return;
    //}

    int domainIndex = 1;
    Domain *domain  = emodel->giveDomain(domainIndex);

    // Output header
    fprintf(FID, "\n %%%% Export of internal variables in integration points \n\n");
    fprintf(FID, "\n %% for interpretation of internal var. numbers see internalstatetype.h\n");


    int numVars = this->internalVarsToExport.giveSize();
    // Output the internalVarsToExport-list
    fprintf(FID, "\tIntegrationPointFields.InternalVarsToExport = [");
    for ( int i = 1; i <= numVars; i++ ) {
        fprintf( FID, "%i ", this->internalVarsToExport.at(i) );
    }
    fprintf(FID, "];\n");





    FloatArray valueArray;

    int nelem = this->elList.giveSize();
    if ( nelem == 0 ) { // no list given, export all elements
        nelem = domain->giveNumberOfElements();
        this->elList.resize(nelem);
        for ( int i = 1; i <= nelem; i++ ) {
            this->elList.at(i) = i;
        }
    }

    fprintf(FID, "\tIntegrationPointFields.Elements = cell(%i,1); \n", nelem);

    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        Element *el = domain->giveElement( this->elList.at(ielem) );
        fprintf( FID, "\tIntegrationPointFields.Elements{%i}.elementNumber = %i; \n", ielem, el->giveNumber() );

        int numIntRules = el->giveNumberOfIntegrationRules();
        fprintf(FID, "\tIntegrationPointFields.Elements{%i}.integrationRule = cell(%i,1); \n", ielem, numIntRules);
        for ( int i = 1; i <= numIntRules; i++ ) {
            IntegrationRule *iRule = el->giveIntegrationRule(i - 1);

            fprintf( FID, "\tIntegrationPointFields.Elements{%i}.integrationRule{%i}.ip = cell(%i,1); \n ",
                    ielem, i, iRule->giveNumberOfIntegrationPoints() );

            // Loop over integration points
            for ( int j = 1; j <= iRule->giveNumberOfIntegrationPoints(); j++ ) {
                IntegrationPoint *ip = iRule->getIntegrationPoint(j - 1);

                double weight = ip->giveWeight();

                fprintf(FID, "\tIntegrationPointFields.Elements{%i}.integrationRule{%i}.ip{%i}.ipWeight = %e; \n ",
                        ielem, i, j, weight);


                // export Gauss point coordinates
                fprintf(FID, "\tIntegrationPointFields.Elements{%i}.integrationRule{%i}.ip{%i}.coords = [",
                        ielem, i, j);

                FloatArray coords;
                el->computeGlobalCoordinates( coords, * ( ip->giveCoordinates() ) );
                for ( int ic = 1; ic <= coords.giveSize(); ic++ ) {
                    fprintf( FID, "%e ", coords.at(ic) );
                }
                fprintf(FID, "]; \n");

                // export internal variables
                fprintf(FID, "\tIntegrationPointFields.Elements{%i}.integrationRule{%i}.ip{%i}.valArray = cell(%i,1); \n",
                        ielem, i, j, numVars);

                for ( int iv = 1; iv <= numVars; iv++ ) {
                    fprintf(FID, "\tIntegrationPointFields.Elements{%i}.integrationRule{%i}.ip{%i}.valArray{%i} = [",
                            ielem, i, j, iv);
                    InternalStateType vartype = ( InternalStateType ) this->internalVarsToExport.at(iv);
                    el->giveIPValue(valueArray, ip, vartype, tStep);
                    int nv = valueArray.giveSize();
                    for ( int ic = 1; ic <= nv; ic++ ) {
                        fprintf( FID, "%.6e ", valueArray.at(ic) );
                    }
                    fprintf(FID, "]; \n");
                }
            }
        }
    }
}


void
MatlabExportModule :: initialize()
{ }


void
MatlabExportModule :: terminate()
{ }


FILE *
MatlabExportModule :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;
    std :: ostringstream baseFileName;
    std :: string fileName;

    fileName = this->emodel->giveOutputBaseFileName();

    size_t foundDot;
    foundDot = fileName.rfind(".");
    fileName.replace(foundDot, 1, "_");

    char fext [ 100 ];
    if ( this->testSubStepOutput() ) {
        // include tStep version in output file name
#ifdef __PARALLEL_MODE
        if ( this->emodel->isParallel() && this->emodel->giveNumberOfProcesses() > 1 ) {
            sprintf( fext, "_%03d_m%d_%d_%d", emodel->giveRank(), this->number, tStep->giveNumber(), tStep->giveSubStepNumber() );
        } else
#endif
        sprintf( fext, "_m%d_%d_%d", this->number, tStep->giveNumber(), tStep->giveSubStepNumber() );
    } else   {
#ifdef __PARALLEL_MODE
        if ( this->emodel->isParallel() && this->emodel->giveNumberOfProcesses() > 1 ) {
            sprintf( fext, "_%03d_m%d_%d", emodel->giveRank(), this->number, tStep->giveNumber() );
        } else
#endif
        sprintf( fext, "_m%d_%d", this->number, tStep->giveNumber() );
    }

    fileName += fext;

    functionname = fileName;
    fileName += ".m";

    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }

    return answer;
}
} // end namespace oofem
