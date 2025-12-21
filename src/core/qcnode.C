
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


#include "qcnode.h"
#include "sm/EngineeringModels/qclinearstatic.h"
#include "hangingnode.h"
#include "slavedof.h"
#include "floatarray.h"
#include "intarray.h"
#include "element.h"
#include "feinterpol.h"
#include "spatiallocalizer.h"
#include "classfactory.h"
#include "engngm.h"
#include "paramkey.h"

#include "outputmanager.h"

namespace oofem {
REGISTER_DofManager(qcNode);

ParamKey qcNode::IPK_qcNode_masterElement("masterElement");
ParamKey qcNode::IPK_qcNode_masterRegion("masterRegion");

qcNode :: qcNode(int n, Domain *aDomain) : Node(n, aDomain)
{
#ifdef __OOFEG
    initialized = false;
#endif
}

void qcNode :: initializeFrom(InputRecord &ir, int priority)
{
    Node :: initializeFrom(ir, priority);

    ParameterManager &ppm =  this->giveDomain()->dofmanPPM;

    PM_UPDATE_PARAMETER(masterElement, ppm, ir, this->number, IPK_qcNode_masterElement, priority) ;
    PM_UPDATE_PARAMETER(masterRegion, ppm, ir, this->number, IPK_qcNode_masterRegion, priority) ;


}

int qcNode :: checkConsistency()
{
    // for all nodes
    int result = Node :: checkConsistency();

    // for hanging nodes only
    if ( this->qcNodeTypeLabel == 2 ) {
        Element *e = this->domain->giveElement(this->masterElement);

#if 0
        // Check if master is in same mode
        if ( parallel_mode != DofManager_local ) {
            for ( int i = 1; i <= countOfMasterNodes; i++ ) {
                if ( e->giveNode(i)->giveParallelMode() != parallel_mode ) {
                    OOFEM_WARNING("Mismatch in parallel mode of qcNode and master");
                    return false;
                }
            }
        }
#endif

        // Check local coordinate systems
        for ( int i = 1; i <= e->giveNumberOfNodes(); ++i ) {
            if ( !this->hasSameLCS( e->giveNode(i) ) ) {
                OOFEM_WARNING("Different lcs for master/slave nodes.");
                result = false;
            }
        }
    }
    return result;
}

void qcNode :: postInitialize()
{
    #ifdef __SM_MODULE
    QClinearStatic *em = dynamic_cast< QClinearStatic * >( this->giveDomain()->giveEngngModel() );
    if ( em ) {
        if ( em->giveQcApproachNumber() == 0 ) {
            this->setAsRepnode();
        } else {
            // set node type according to fullsolveddomain
            if ( this->initializeAsRepnode() ) {
                this->setAsRepnode();
            } else {
                this->setAsHanging();
            }
        }
    } else {
        OOFEM_ERROR("\"qcNode\" can be used only in \"QClinearStatic\" EngngModel");
    }
#else
    OOFEM_ERROR("\"qcNode\" can be used only in \"QClinearStatic\" EngngModel");
#endif


    Node :: postInitialize();
    if ( this->qcNodeTypeLabel == 2 ) {
        this->postInitializeAsHangingNode();
    }
}


void qcNode :: postInitializeAsHangingNode()
{
    // Node :: postInitialize();

    Element *e;
    FEInterpolation *fei;
    FloatArray lcoords, masterContribution;

#ifdef __OOFEG
    if ( initialized ) {
        return;
    }
    initialized = true;
#endif

    // First check element and interpolation
    if ( masterElement == -1 ) { // Then we find it by taking the closest (probably containing element)
        FloatArray closest;
        SpatialLocalizer *sp = this->domain->giveSpatialLocalizer();
        sp->init();
        // Closest point or containing point? It should be contained, but with numerical errors it might be slightly outside
        // so the closest point is more robust.
        if ( !( e = sp->giveElementClosestToPoint(lcoords, closest, coordinates, this->masterRegion) ) ) {
            OOFEM_ERROR("Couldn't find closest element (automatically).");
        }
        this->masterElement = e->giveNumber();
    } else if ( !( e = this->giveDomain()->giveElement(this->masterElement) ) ) {
        OOFEM_ERROR("Requested element %d doesn't exist.", this->masterElement);
    }
    if ( !( fei = e->giveInterpolation() ) ) {
        OOFEM_ERROR("Requested element %d doesn't have a interpolator.", this->masterElement);
    }

    if ( lcoords.giveSize() == 0 ) { // we don't need to do this again if the spatial localizer was used.
        fei->global2local( lcoords, coordinates, FEIElementGeometryWrapper(e) );
    }

    // Initialize slave dofs (inside check of consistency of receiver and master dof)
    const IntArray &masterNodes = e->giveDofManArray();
    for ( Dof *dof : *this ) {
        SlaveDof *sdof = dynamic_cast< SlaveDof * >(dof);
        if ( sdof ) {
            DofIDItem id = sdof->giveDofID();
            fei = e->giveInterpolation(id);
            if ( !fei ) {
                OOFEM_ERROR("Requested interpolation for dof id %d doesn't exist in element %d.",
                            id, this->masterElement);
            }
#if 0 // This won't work (yet), as it requires some more general FEI classes, or something similar.
            if ( fei->hasMultiField() ) {
                FloatMatrix multiContribution;
                IntArray masterDofIDs, masterNodesDup, dofids;
                fei->evalMultiN(multiContribution, dofids, lcoords, FEIElementGeometryWrapper(e), 0.0);
                masterContribution.flatten(multiContribution);
                masterDofIDs.clear();
                for ( int i = 0; i <= multiContribution.giveNumberOfColumns(); ++i ) {
                    masterDofIDs.followedBy(dofids);
                    masterNodesDup.followedBy(masterNodes);
                }
                sdof->initialize(masterNodesDup, & masterDofIDs, masterContribution);
            } else { }
#else
            // Note: There can be more masterNodes than masterContributions, since all the
            // FEI classes are based on that the first nodes correspond to the simpler/linear interpolation.
            // If this assumption is changed in FEIElementGeometryWrapper + friends,
            // masterNode will also need to be modified for each dof accordingly.
            fei->evalN( masterContribution, lcoords, FEIElementGeometryWrapper(e) );
            sdof->initialize(masterNodes, IntArray(), masterContribution);
#endif
        }
    }
}

/* Test if node will be initialized as repnode  i.e. if node is in fullsolved domain or if has BC or IC
 */
bool
qcNode :: initializeAsRepnode()
{
    // Nodes with prescribed BC and IC can be set as automatically repnodes here
    /*
     * // if node has BC
     * if (dofBCmap){
     * if (dofBCmap->size()!=0) {
     *  return true;
     * }
     * }
     * // if node has IC
     * if (dofICmap){
     * if (dofICmap->size()!=0) {
     *  return true;
     * }
     * }
     */

    // if node is in fullsolved domain
  
#ifdef __SM_MODULE
    QClinearStatic *em = dynamic_cast<  QClinearStatic * >( this->giveDomain()->giveEngngModel() );
    if ( !em ) {
        OOFEM_ERROR("qcNode is used in unsupported Engineering Models");
    }
    
    if ( em->nodeInFullSolvedDomainTest(this) ) {
        return true;
    }

    //else
    return false;

#else
    OOFEM_ERROR("qcNode is used in unsupported Engineering Models");
    return false;
#endif

}

void qcNode :: setAsRepnode()
{
    // delete content of DofTypeMap (if exist) // (set all doftype=0 is not enough)
    if ( this->giveDofTypeMap() !=  NULL ) {
        this->giveDofTypeMap()->clear();
    }
    this->qcNodeTypeLabel = 1;
}
void qcNode :: setAsHanging()
{
    // set all doftype=2 in dofTypemap

    if ( this->giveDofTypeMap() == NULL ) {
        this->dofTypemap.clear();
    }

    int DofTypeMapSize = this->giveDofTypeMap()->size();
    // insert "2" into new (empty) dofTypemap
    if ( DofTypeMapSize == 0 ) {
        for ( int i = 1; i <= this->giveDomain()->giveDefaultNodeDofIDArry().giveSize(); i++ ) {
            this->giveDofTypeMap()->insert( std :: pair< int, int >(i, 2) );
        }
    }
    // rewrite old dofTypemap by "2"
    else {
        for ( int i = 1; i <= DofTypeMapSize; i++ ) {
            this->giveDofTypeMap()->at(i) = 2;
        }
    }

    this->qcNodeTypeLabel = 2;
}


void qcNode :: printOutputAt(FILE *stream, TimeStep *tStep)
{
    EngngModel *emodel = this->giveDomain()->giveEngngModel();

    if ( this->giveQcNodeType() == 1 ) {
        fprintf( stream, "%-8s R%8d (%8d):\n", this->giveClassName(), this->giveLabel(), this->giveNumber() );
        for ( Dof *dof : *this ) {
            emodel->printDofOutputAt(stream, dof, tStep);
        }
    } else if ( this->giveQcNodeType() == 2 ) {
        fprintf( stream, "%-8s H%8d (%8d): el. %8d\n", this->giveClassName(), this->giveLabel(), this->giveNumber(), this->giveMasterElementNumber() );
        for ( Dof *dof : *this ) {
            emodel->printDofOutputAt(stream, dof, tStep);
        }
    } else {
        OOFEM_WARNING( "Node %d cannot be printed out: unknown QcNodeType", this->giveGlobalNumber() );
    }
}
} // end namespace oofem
