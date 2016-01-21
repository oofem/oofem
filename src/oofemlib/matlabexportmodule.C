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
#include "set.h"
#include "unknownnumberingscheme.h"
#include "prescribedmean.h"
#include "feinterpol.h"
    
#ifdef __FM_MODULE
#include "../fm/tr21stokes.h"
#include "../fm/tet21stokes.h"
#include "../fm/stokesflow.h"
#endif

#ifdef __SM_MODULE
#include "../sm/Elements/nlstructuralelement.h"
#include "../sm/EngineeringModels/structengngmodel.h"
#endif


namespace oofem {

REGISTER_ExportModule( MatlabExportModule )

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
    reactionForcesNodeSet = 0;
    IPFieldsElSet = 0;
    noscaling = false;
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
    exportHomogenizeIST = ir->hasField(_IFT_MatlabExportModule_homogenizeInternalVars);


    exportReactionForces = ir->hasField(_IFT_MatlabExportModule_ReactionForces);
    reactionForcesDofManList.resize(0);
    if ( exportReactionForces ) {
        IR_GIVE_OPTIONAL_FIELD(ir, reactionForcesDofManList, _IFT_MatlabExportModule_DofManList);
        this->reactionForcesNodeSet = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, this->reactionForcesNodeSet, _IFT_MatlabExportModule_ReactionForcesNodeSet);
    }

    exportIntegrationPointFields = ir->hasField(_IFT_MatlabExportModule_IntegrationPoints);
    elList.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, _IFT_MatlabExportModule_internalVarsToExport);
    if ( exportIntegrationPointFields ) {
        IR_GIVE_OPTIONAL_FIELD(ir, elList, _IFT_MatlabExportModule_ElementList);
        IPFieldsElSet = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, IPFieldsElSet, _IFT_MatlabExportModule_IPFieldsElSet);
    }

    noscaling = ir->hasField(_IFT_MatlabExportModule_noScaledHomogenization);

    return ExportModule :: initializeFrom(ir);
}



void
MatlabExportModule :: computeArea(TimeStep *tStep)
{
    Domain *domain = emodel->giveDomain(1);

    smax.clear();
    smin.clear();

    for (int i = 1; i <= domain->giveNumberOfSpatialDimensions(); i++) {
        smax.push_back(domain->giveDofManager(1)->giveCoordinate(i));
        smin.push_back(domain->giveDofManager(1)->giveCoordinate(i));
    }

    for ( int i = 0; i < domain->giveNumberOfDofManagers(); i++ ) {
        for (int j = 0; j < domain->giveNumberOfSpatialDimensions(); j++) {
            smax.at(j)=max(smax.at(j), domain->giveDofManager(i+1)->giveCoordinate(j+1));
            smin.at(j)=min(smin.at(j), domain->giveDofManager(i+1)->giveCoordinate(j+1));
        }
    }

    Area = 0;
    Volume = 0;

    if ( domain->giveNumberOfSpatialDimensions() == 2 ) {
        for ( auto &elem : domain->giveElements() ) {
            Area += elem->computeArea();
        }
    } else {

        for (size_t i = 0; i < partVolume.size(); i++ ) {
            partVolume.at(i) = 0.0;
        }

        for ( auto &elem : domain->giveElements() ) {
            //
            double v;
#ifdef __SM_MODULE
            if ( NLStructuralElement *e = dynamic_cast< NLStructuralElement *>( elem.get() ) ) {
                v = e->computeCurrentVolume(tStep);
            } else {
#endif
                v = elem->computeVolume();
#ifdef __SM_MODULE
            }
#endif

            std :: string eName ( elem->giveClassName() );
            int j = -1;

            //printf("%s\n", eName.c_str());

            for ( size_t k = 0; k < partName.size(); k++ ) {
                //printf("partName.at(%u) = %s\n", k, partName.at(k).c_str() );
                if ( eName.compare(partName.at(k)) == 0 ) {
                    j = k;
                    break;
                }
            }

            if ( j == -1 ) {
                partName.push_back( elem->giveClassName() );
                partVolume.push_back( v );
            } else {
                partVolume.at(j) += v;
            }

            Volume += v;

        }
    }

}


void
MatlabExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }


    int nelem = this->elList.giveSize();
    if ( nelem == 0 ) { // no list given, export all elements
        this->elList.enumerate(this->emodel->giveDomain(1)->giveNumberOfElements());
    }

    FILE *FID;
    FID = giveOutputStream(tStep);
    Domain *domain  = emodel->giveDomain(1);
    ndim=domain->giveNumberOfSpatialDimensions();

    // Output header
    fprintf( FID, "%%%% OOFEM generated export file \n");
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
        computeArea(tStep);
        fprintf(FID, "\tarea.xmax=%f;\n", smax.at(0));
        fprintf(FID, "\tarea.xmin=%f;\n", smin.at(0));
        fprintf(FID, "\tarea.ymax=%f;\n", smax.at(1));
        fprintf(FID, "\tarea.ymin=%f;\n", smin.at(1));
        if ( ndim == 2 ) {
            fprintf(FID, "\tarea.area=%f;\n", Area);
            fprintf(FID, "\tvolume=[];\n");
        } else {
            fprintf(FID, "\tarea.zmax=%f;\n", smax.at(2));
            fprintf(FID, "\tarea.zmin=%f;\n", smin.at(2));
            fprintf(FID, "\tarea.area=[];\n");
            fprintf(FID, "\tarea.volume=%f;\n", Volume);
            for (size_t i=0; i<this->partName.size(); i++) {
                fprintf(FID, "\tarea.volume_%s=%f;\n", partName.at(i).c_str(), partVolume.at(i));
            }
        }
    } else {
        fprintf(FID, "\tarea.area=[];\n");
        fprintf(FID, "\tarea.volume=[];\n");
    }

    if ( exportSpecials ) {
        if ( !exportArea ) {
            computeArea(tStep);
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

    // Homogenized quantities
    if ( exportHomogenizeIST ) {
        doOutputHomogenizeDofIDs(tStep, FID);
    }

    fprintf(FID, "\nend\n");
    fclose(FID);
}


void
MatlabExportModule :: doOutputMesh(TimeStep *tStep, FILE *FID)
{
    Domain *domain  = emodel->giveDomain(1);

    fprintf(FID, "\tmesh.p=[");
    for ( auto &dman : domain->giveDofManagers() ) {
        for ( int j = 1; j <= domain->giveNumberOfSpatialDimensions(); j++) {
            double c = dman->giveCoordinate(j);
            fprintf(FID, "%f, ", c);
        }
        fprintf(FID, "; ");
    }

    fprintf(FID, "]';\n");

    int numberOfDofMans=domain->giveElement(1)->giveNumberOfDofManagers();

    fprintf(FID, "\tmesh.t=[");
    for ( auto &elem : domain->giveElements() ) {
        if ( elem->giveNumberOfDofManagers() == numberOfDofMans ) {
            for ( int j = 1; j <= elem->giveNumberOfDofManagers(); j++ ) {
                fprintf( FID, "%d,", elem->giveDofManagerNumber(j) );
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
    std :: vector< std :: vector< double > >valuesList;

    for ( auto &dman : domain->giveDofManagers() ) {
        for ( Dof *thisDof: *dman ) {
            it = std :: find( DofIDList.begin(), DofIDList.end(), thisDof->giveDofID() );

            double value = thisDof->giveUnknown(VM_Total, tStep);
            if ( it == DofIDList.end() ) {
                DofIDList.push_back( thisDof->giveDofID() );
                valuesList.push_back({value});
            } else {
                std::size_t pos = it - DofIDList.begin();
                valuesList[pos].push_back(value);
            }
        }
    }

    fprintf(FID, "\tdata.DofIDs=[");
    for ( auto &dofid : DofIDList ) {
        fprintf( FID, "%d, ", dofid );
    }

    fprintf(FID, "];\n");

    for ( size_t i = 0; i < valuesList.size(); i++ ) {
        fprintf(FID, "\tdata.a{%lu}=[", static_cast< long unsigned int >(i) + 1);
        for ( double val: valuesList[i] ) {
            fprintf( FID, "%f,", val );
        }

        fprintf(FID, "];\n");
    }

}


void
MatlabExportModule :: doOutputSpecials(TimeStep *tStep,    FILE *FID)
{
//    FloatArray v_hat, GradPTemp, v_hatTemp;

    Domain *domain  = emodel->giveDomain(1);
/*
    v_hat.clear();

        ///@todo Sort out this hack in a nicer (modular) way / Mikael
#if 0
    for ( auto &elem : domain->giveElements() ) {
#ifdef __FM_MODULE

        if ( Tr21Stokes *Tr = dynamic_cast< Tr21Stokes * >( elem.get() ) ) {
            Tr->giveIntegratedVelocity(v_hatTemp, tStep);
            v_hat.add(v_hatTemp);
        } else if ( Tet21Stokes *Tet = dynamic_cast< Tet21Stokes * >( elem.get() ) ) {
            Tet->giveIntegratedVelocity(v_hatTemp, tStep);
            v_hat.add(v_hatTemp);
        }

#endif
    }
#endif

    // Compute intrinsic area/volume
    double intrinsicSize = 1.0;

    std :: vector <double> V;

    for ( int i = 0; i < (int)smax.size(); i++ ) {
        intrinsicSize *= ( smax.at(i) - smin.at(i) );
    }

    for ( double vh: v_hat ) {
        V.push_back(vh);
    }

    fprintf(FID, "\tspecials.velocitymean=[");
    if (V.size()>0) {
        for (int i=0; i<ndim; i++) {
            fprintf(FID, "%e", V.at(i));
            if (i!=(ndim-1)) fprintf (FID, ", ");
        }
        fprintf(FID, "];\n");
    } else {
        fprintf(FID, "]; %% No velocities\n");
    }

    */

    // Output weak periodic boundary conditions
    unsigned int wpbccount = 1, sbsfcount = 1, mcount = 1;

    for ( auto &gbc : domain->giveBcs() ) {
        WeakPeriodicBoundaryCondition *wpbc = dynamic_cast< WeakPeriodicBoundaryCondition * >( gbc.get() );
        if ( wpbc ) {
            for ( int j = 1; j <= wpbc->giveNumberOfInternalDofManagers(); j++ ) {
                fprintf(FID, "\tspecials.weakperiodic{%u}.descType=%u;\n", wpbccount, wpbc->giveBasisType() );
                fprintf(FID, "\tspecials.weakperiodic{%u}.coefficients=[", wpbccount);
                for ( Dof *dof: *wpbc->giveInternalDofManager(j) ) {
                    double X = dof->giveUnknown(VM_Total, tStep);
                    fprintf(FID, "%e\t", X);
                }

                fprintf(FID, "];\n");
                wpbccount++;
            }
        }
        SolutionbasedShapeFunction *sbsf = dynamic_cast< SolutionbasedShapeFunction *>( gbc.get());
        if (sbsf) {
            fprintf(FID, "\tspecials.solutionbasedsf{%u}.values=[", sbsfcount);
            for ( Dof *dof: *sbsf->giveInternalDofManager(1) ) {                  // Only one internal dof manager
                double X = dof->giveUnknown(VM_Total, tStep);
                fprintf(FID, "%e\t", X);
            }
            fprintf(FID, "];\n");
            sbsfcount++;
        }
        PrescribedMean *m = dynamic_cast<PrescribedMean *> ( gbc.get() );
        if (m) {
            fprintf(FID, "\tspecials.prescribedmean{%u}.value=[", mcount);
            for ( Dof *dof: *m->giveInternalDofManager(1)) {
                double X = dof->giveUnknown(VM_Total, tStep);
                fprintf(FID, "%e\t", X);
            }
            fprintf(FID, "];\n");
            mcount++;
        }
    }
}


void
MatlabExportModule :: doOutputReactionForces(TimeStep *tStep,    FILE *FID)
{

    int domainIndex = 1;
    Domain *domain  = emodel->giveDomain( domainIndex );

    FloatArray reactions;
    IntArray dofManMap, dofidMap, eqnMap;
#ifdef __SM_MODULE
    StructuralEngngModel *strEngMod = dynamic_cast< StructuralEngngModel * >(emodel);
    if ( strEngMod ) {
        strEngMod->buildReactionTable(dofManMap, dofidMap, eqnMap, tStep, domainIndex);
        strEngMod->computeReaction(reactions, tStep, 1);
    } else
#endif
    {
        OOFEM_ERROR("Cannot export reaction forces - only implemented for structural problems.");
    }

    // Set the nodes and elements to export based on sets
    if ( this->reactionForcesNodeSet > 0 ) {
        Set *set = domain->giveSet( this->reactionForcesNodeSet );
        reactionForcesDofManList = set->giveNodeList();
    }


    int numDofManToExport = this->reactionForcesDofManList.giveSize();
    if ( numDofManToExport == 0 ) { // No dofMan's given - export every dMan with reaction forces

        for (int i = 1; i <= domain->giveNumberOfDofManagers(); i++) {
            if ( dofManMap.contains(i) ) {
                this->reactionForcesDofManList.followedBy(i);
            }
        }
        numDofManToExport = this->reactionForcesDofManList.giveSize();
    }


    // Output header
    fprintf( FID, "\n %%%% Export of reaction forces \n\n" );

    // Output the dofMan numbers that are exported
    fprintf( FID, "\tReactionForces.DofManNumbers = [" );
    for ( int i = 1; i <= numDofManToExport; i++ ) {
        fprintf( FID, "%i ", this->reactionForcesDofManList.at(i) );
    }
    fprintf( FID, "];\n" );


    // Define the reaction forces as a cell object
    fprintf( FID, "\tReactionForces.ReactionForces = cell(%i,1); \n", numDofManToExport );
    fprintf( FID, "\tReactionForces.DofIDs = cell(%i,1); \n", numDofManToExport );


    // Output the reaction forces for each dofMan. If a certain dof is not prescribed zero is exported.
    IntArray dofIDs;
    for ( int i = 1; i <= numDofManToExport; i++ ) {
        int dManNum = this->reactionForcesDofManList.at(i);

        fprintf(FID, "\tReactionForces.ReactionForces{%i} = [", i);
        if ( dofManMap.contains( dManNum ) ) {

            DofManager *dofMan = domain->giveDofManager( dManNum );
            dofIDs.clear();

            for ( Dof *dof: *dofMan ) {
                int num = dof->giveEquationNumber( EModelDefaultPrescribedEquationNumbering() );
                int pos = eqnMap.findFirstIndexOf( num );
                dofIDs.followedBy(dof->giveDofID());
                if ( pos > 0 ) {
                    fprintf(FID, "%e ", reactions.at(pos));
                } else {
                    fprintf( FID, "%e ", 0.0 ); // if not prescibed output zero
                }
            }
        }
        fprintf(FID, "];\n");

        // Output dof ID's

        fprintf( FID, "\tReactionForces.DofIDs{%i} = [", i);
        if ( dofManMap.contains( dManNum ) ) {
            for ( int id: dofIDs ) {
                fprintf( FID, "%i ", id );
            }
        }
        fprintf(FID, "];\n");
    }
}


void
MatlabExportModule :: doOutputIntegrationPointFields(TimeStep *tStep,    FILE *FID)
{

    int domainIndex = 1;
    Domain *domain  = emodel->giveDomain( domainIndex );

    // Output header
    fprintf( FID, "\n %%%% Export of internal variables in integration points \n\n" );
    fprintf( FID, "\n %% for interpretation of internal var. numbers see internalstatetype.h\n");


    int numVars = this->internalVarsToExport.giveSize();
    // Output the internalVarsToExport-list
    fprintf( FID, "\tIntegrationPointFields.InternalVarsToExport = [" );
    for ( int i = 1; i <= numVars; i++ ) {
        fprintf( FID, "%i ", this->internalVarsToExport.at(i) );
    }
    fprintf( FID, "];\n" );




    if ( this->IPFieldsElSet > 0 ) {
        Set *set = domain->giveSet( this->IPFieldsElSet );
        elList = set->giveElementList();
    }


    FloatArray valueArray;

    int nelem = this->elList.giveSize();

    fprintf( FID, "\tIntegrationPointFields.Elements = cell(%i,1); \n", nelem );

    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        Element *el = domain->giveElement( this->elList.at(ielem) );
        fprintf( FID, "\tIntegrationPointFields.Elements{%i}.elementNumber = %i; \n", ielem, el->giveNumber());

        int numIntRules = el->giveNumberOfIntegrationRules();
        fprintf( FID, "\tIntegrationPointFields.Elements{%i}.integrationRule = cell(%i,1); \n", ielem, numIntRules);
        for ( int i = 1; i <= numIntRules; i++ ) {
            IntegrationRule *iRule = el->giveIntegrationRule(i-1);

            fprintf( FID, "\tIntegrationPointFields.Elements{%i}.integrationRule{%i}.ip = cell(%i,1); \n ",
                     ielem, i, iRule->giveNumberOfIntegrationPoints() );

            // Loop over integration points
            for ( GaussPoint *ip: *iRule ) {

                double weight = ip->giveWeight();

                fprintf( FID, "\tIntegrationPointFields.Elements{%i}.integrationRule{%i}.ip{%i}.ipWeight = %e; \n ",
                         ielem, i, ip->giveNumber(), weight);


                // export Gauss point coordinates
                fprintf( FID, "\tIntegrationPointFields.Elements{%i}.integrationRule{%i}.ip{%i}.coords = [",
                         ielem, i, ip->giveNumber());

                FloatArray coords;
                el->computeGlobalCoordinates( coords, ip->giveNaturalCoordinates() );
                for ( int ic = 1; ic <= coords.giveSize(); ic++ ) {
                    fprintf( FID, "%e ", coords.at(ic) );
                }
                fprintf( FID, "]; \n" );

                // export internal variables
                fprintf( FID, "\tIntegrationPointFields.Elements{%i}.integrationRule{%i}.ip{%i}.valArray = cell(%i,1); \n",
                         ielem, i, ip->giveNumber(), numVars);

                for ( int iv = 1; iv <= numVars; iv++ ) {
                    fprintf( FID, "\tIntegrationPointFields.Elements{%i}.integrationRule{%i}.ip{%i}.valArray{%i} = [",
                             ielem, i, ip->giveNumber(), iv);
                    InternalStateType vartype = ( InternalStateType ) this->internalVarsToExport.at(iv);
                    el->giveIPValue(valueArray, ip, vartype, tStep);
                    int nv = valueArray.giveSize();
                    for ( int ic = 1; ic <= nv; ic++ ) {
                        fprintf( FID, "%.6e ", valueArray.at(ic) );
                    }
                    fprintf( FID, "]; \n" );
                }
            }

        }
    }

}


void
MatlabExportModule :: initialize()
{ 
    ExportModule :: initialize();
}


void
MatlabExportModule :: terminate()
{ }


FILE *
MatlabExportModule :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;

    char fext[100];
    sprintf( fext, "_m%d_%d", this->number, tStep->giveNumber() );

    if ( this->testSubStepOutput() ) {
        // include tStep version in output file name
        if ( this->emodel->isParallel() && this->emodel->giveNumberOfProcesses() > 1 ) {
            sprintf( fext, "_%03d_m%d_%d_%d", emodel->giveRank(), this->number, tStep->giveNumber(), tStep->giveSubStepNumber() );
        } else {
            sprintf( fext, "_m%d_%d_%d", this->number, tStep->giveNumber(), tStep->giveSubStepNumber() );
        }
    } else   {
        if ( this->emodel->isParallel() && this->emodel->giveNumberOfProcesses() > 1 ) {
            sprintf( fext, "_%03d_m%d_%d", emodel->giveRank(), this->number, tStep->giveNumber() );
        } else {
            sprintf( fext, "_m%d_%d", this->number, tStep->giveNumber() );
        }
    }

    std :: string fileName;

    fileName = this->emodel->giveOutputBaseFileName();

    size_t foundDot;
    foundDot = fileName.rfind(".");

    // CARL
    while (foundDot != std :: string :: npos) {
        fileName.replace(foundDot, 1, "_");
        foundDot = fileName.rfind(".");
//        printf("%s\n", fileName.c_str());
    }

    fileName += fext;

    std :: string temp;
    temp = fileName;
    size_t backslash = temp.rfind("/");

    if (backslash != std :: string :: npos ) {
        functionname = temp.substr(backslash+1, std :: string :: npos);
    } else {
        functionname = temp;
    }

    fileName += ".m";

    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }

    return answer;
}

void
MatlabExportModule :: doOutputHomogenizeDofIDs(TimeStep *tStep,    FILE *FID)
{

    std :: vector <FloatArray*> HomQuantities;
    double Vol = 0.0;

    // Initialize vector of arrays constaining homogenized quantities
    HomQuantities.resize(internalVarsToExport.giveSize());

    for (int j=0; j<internalVarsToExport.giveSize(); j++) {
        HomQuantities.at(j) = new FloatArray;
    }

    int nelem = this->elList.giveSize();
    for (int i = 1; i<=nelem; i++) {
        Element *e = this->emodel->giveDomain(1)->giveElement(elList.at(i));
        FEInterpolation *Interpolation = e->giveInterpolation();

        Vol = Vol + e->computeVolumeAreaOrLength();

        for ( GaussPoint *gp: *e->giveDefaultIntegrationRulePtr() ) {

            for (int j=0; j<internalVarsToExport.giveSize(); j++) {
                FloatArray elementValues;
                e->giveIPValue(elementValues, gp, (InternalStateType) internalVarsToExport(j), tStep);
                double detJ=fabs(Interpolation->giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(e)));

                elementValues.times(gp->giveWeight()*detJ);
                if (HomQuantities.at(j)->giveSize() == 0) {
                    HomQuantities.at(j)->resize(elementValues.giveSize());
                    HomQuantities.at(j)->zero();
                };
                HomQuantities.at(j)->add(elementValues);
            }
        }
    }


    if (noscaling) Vol=1.0;

    for ( std :: size_t i = 0; i < HomQuantities.size(); i ++) {
        FloatArray *thisIS;
        thisIS = HomQuantities.at(i);
        thisIS->times(1.0/Vol);
        fprintf(FID, "\tspecials.%s = [", __InternalStateTypeToString ( InternalStateType (internalVarsToExport(i)) ) );

        for (int j = 0; j<thisIS->giveSize(); j++) {
            fprintf(FID, "%e", thisIS->at(j+1));
            if (j!=(thisIS->giveSize()-1) ) {
                fprintf(FID, ", ");
            }
        }
        fprintf(FID, "];\n");
        delete HomQuantities.at(i);
    }

}

} // end namespace oofem
