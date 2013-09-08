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


#include "xfemstatic.h"
#include "timestep.h"
#include "metastep.h"
#include "dictionary.h"
#include "classfactory.h"
#include "dofmanager.h"
#include "xfemelementinterface.h"
#include "element.h"

namespace oofem {

REGISTER_EngngModel( XFEMStatic );

XFEMStatic::XFEMStatic(int i, EngngModel *_master):
NonLinearStatic(i, _master),
updateStructureFlag(false)
{
	printf("Entering XFEMStatic::XFEMStatic(int i, EngngModel *_master).\n");
}

XFEMStatic::~XFEMStatic() {

}

void
XFEMStatic :: solveYourselfAt(TimeStep *tStep)
{

    // Initiates the total displacement to zero in the UnknownsDictionary at the first time step.
    // this must be done in order to support a dynamic equation system
    if ( tStep->isTheFirstStep() ) {
        printf("Initializing DofUnknownsDictionary... \n");
        this->initializeDofUnknownsDictionary(tStep);
    }


//    printf("Before: ");
//    initialLoadVector.printYourself();
    // Initialization
    int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()); // 1 stands for domain?
    if ( totalDisplacement.giveSize() != neq ) {

    	printf("Increasing size of displacement array.\n");

//        totalDisplacement.resize(neq);
//        totalDisplacement.zero();
        if ( totalDisplacement.giveSize() != neq ) {
        	FloatArray totalDisplacementNew;
            setValsFromDofMap(totalDisplacementNew, totalDisplacement);
            totalDisplacement = totalDisplacementNew;
        }


        if ( incrementOfDisplacement.giveSize() != neq ) {
        	FloatArray incrementOfDisplacementNew;
            setValsFromDofMap(incrementOfDisplacementNew, incrementOfDisplacement);
            incrementOfDisplacement = incrementOfDisplacementNew;
        }

//        incrementOfDisplacement.resize(neq);
//        incrementOfDisplacement.zero();
        this->setTotalDisplacementFromUnknownsInDictionary(EID_MomentumBalance, VM_Total, tStep);


        if ( incrementalLoadVector.giveSize() != neq ) {
        	FloatArray incrementalLoadVectorNew;
            setValsFromDofMap(incrementalLoadVectorNew, incrementalLoadVector);
            incrementalLoadVector = incrementalLoadVectorNew;
        }

        ///////////////////////////////////////////////////////////////////
        // Map values in the old initialLoadVector to the new initialLoadVector
        if ( initialLoadVector.giveSize() != neq ) {
        	FloatArray initialLoadVectorNew;
            setValsFromDofMap(initialLoadVectorNew, initialLoadVector);
            initialLoadVector = initialLoadVectorNew;
        }

    }

//    printf("After: ");
//    initialLoadVector.printYourself();

#ifdef USE_FRACTURE_MANAGER
    // Instanciate fracture manager
    // should be made in a more proper way with input and the like
    if ( tStep->isTheFirstStep() ) {

        this->fMan = new FractureManager( this->giveDomain(1) );
        this->fMan->failureCriterias = new AList< FailureCriteria >(1); // list of all the criterias to evaluate
        //FailureCriteria *fc = new FailureCriteria(FC_MaxShearStress);
        //FailureCriteria *fc = new FailureCriteria(FC_DamagedNeighborCZ);

        FailureCriteria *fc = new FailureCriteria(FC_DamagedNeighborCZ, this->fMan);


        fc->thresholds.resize(1);
        fc->thresholds.at(1) = -10.0;
        this->fMan->failureCriterias->put(1, fc);
    }


#endif

    this->setUpdateStructureFlag(false);
    NonLinearStatic :: solveYourselfAt(tStep);

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
            initialLoadVector.resize( incrementalLoadVector.giveSize() );
        }
        */
        OOFEM_LOG_DEBUG("Fixed load level\n");

        int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()); // 1 stands for domain?

//        printf("Before: ");
//        incrementalLoadVector.printYourself();
        if ( incrementalLoadVector.giveSize() != neq ) {
            //initialLoadVector.resize( incrementalLoadVector.giveSize() );
//        	initialLoadVector.printYourself();
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
//        	initialLoadVector.printYourself();
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

    #ifdef USE_FRACTURE_MANAGER

    // Fracture/failure mechanics evaluation
    // Upadate components like the XFEM manager and its sub-components
    this->fMan->update(tStep);



    this->setUpdateStructureFlag( this->fMan->giveUpdateFlag() ); // if the internal structure need to be updated
#endif


    // TODO: Check if update is actually needed
    this->setUpdateStructureFlag( true );

    // Update the UnknownsDictionary if needed
    if ( this->needsStructureUpdate() ) {
        printf(" Updating DofUnknownsDictionary... \n");

        buildDofMap();

        for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
            Domain *domain = this->giveDomain(idomain);
            int nnodes = domain->giveNumberOfDofManagers();
            for ( int inode = 1; inode <= nnodes; inode++ ) {
                this->updateDofUnknownsDictionary(domain->giveDofManager(inode), tStep);
            }
        }
    }


    // Update element subdivisions if necessary
    // (e.g. if a crack has moved and cut a new element)
    for( int domInd = 1; domInd <= this->giveNumberOfDomains(); domInd++ ) {

        Domain *domain = this->giveDomain(domInd);
    	int numEl = domain->giveNumberOfElements();

    	for(int i = 1; i <= numEl; i++) {
    		Element *el = domain->giveElement(i);

    		XfemElementInterface *xfemEl = dynamic_cast<XfemElementInterface*> (el);

    		if( xfemEl != NULL ) {
    			xfemEl->recomputeGaussPoints();
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
            //OOFEM_ERROR2( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    } else {
        return NonLinearStatic ::  giveUnknownComponent(mode, tStep, d, dof);
    }

}

void
XFEMStatic :: initializeDofUnknownsDictionary(TimeStep *tStep)
{
    // Initializes all dof values to zero

    Domain *domain;
    Dof *iDof;
    DofManager *node;

    int nDofs;
    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        int nnodes = domain->giveNumberOfDofManagers();
        for ( int inode = 1; inode <= nnodes; inode++ ) {
            node = domain->giveDofManager(inode);
            nDofs = node->giveNumberOfDofs();
            for ( int i = 1; i <= nDofs; i++ ) {
                iDof = node->giveDof(i);
                iDof->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, 0.0);
            }
        }
    }
}

void
XFEMStatic :: setTotalDisplacementFromUnknownsInDictionary(EquationID type, ValueModeType mode, TimeStep *tStep)
{
	printf("Entering XFEMStatic :: setTotalDisplacementFromUnknownsInDictionary().\n");

    // Sets the values in the displacement vector based on stored values in the unknowns dictionaries.
    // Used in the beginning of each time step.
    Domain *domain;
    DofManager *inode;
    Dof *iDof;
    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        for ( int j = 1; j <= domain->giveNumberOfDofManagers(); j++ ) {
            inode = domain->giveDofManager(j);
            int eqNum;
            for ( int i = 1; i <= inode->giveNumberOfDofs(); i++ ) {
                iDof = inode->giveDof(i);
                eqNum = iDof->giveEqn();
                if ( eqNum > 0 ) {
                    double val = iDof->giveUnknown(mode, tStep);
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
    Dof *iDof;
    double val;
    for ( int i = 1; i <= inode->giveNumberOfDofs(); i++ ) {
        iDof = inode->giveDof(i);
        int eqNum = iDof->__giveEquationNumber();
        if ( iDof->hasBc(tStep) ) {
            val = iDof->giveBcValue(VM_Total, tStep);
        } else {
            if ( eqNum > 0 ) {
                val = totalDisplacement.at(eqNum);
            } else { // new eq number
                val = 0.0;
            }
        }

        iDof->updateUnknownsDictionary(tStep, VM_Total, val);
    }
}

void XFEMStatic :: buildDofMap() {

	printf("Building dof map.\n");
	mDofEqnNumMap.clear();

    for ( int domainIndex = 1; domainIndex <= this->giveNumberOfDomains(); domainIndex++ ) {
    	Domain *domain = this->giveDomain(domainIndex);

        for ( int dManIndex = 1; dManIndex <= domain->giveNumberOfDofManagers(); dManIndex++ ) {
        	DofManager *dMan = domain->giveDofManager(dManIndex);

            for ( int k = 1; k <= dMan->giveNumberOfDofs(); k++ ) {

            	Dof *dof = dMan->giveDof(k);
            	int eqNum = dof->giveEqn();

                if ( eqNum > 0 ) {

                	std::vector<int> key(3);
                	key[0] = domainIndex;
                	key[1] = dManIndex;
                	key[2] = k;

                	mDofEqnNumMap[ key ] = eqNum;
                }

            }
        }
    }

}

void XFEMStatic :: setValsFromDofMap(FloatArray &oArray, const FloatArray &iArray) {


    int neq = 0;
    for ( int domainIndex = 1; domainIndex <= this->giveNumberOfDomains(); domainIndex++ ) {
    	neq += this->giveNumberOfDomainEquations(domainIndex, EModelDefaultEquationNumbering());
    }

	int numEqOld = iArray.giveSize();
	printf("Setting values from dof map. neq: %d numEqOld: %d\n", neq, numEqOld);


	oArray.resize( neq );
	oArray.zero();

    for ( int domainIndex = 1; domainIndex <= this->giveNumberOfDomains(); domainIndex++ ) {
    	Domain *domain = this->giveDomain(domainIndex);

        for ( int dManIndex = 1; dManIndex <= domain->giveNumberOfDofManagers(); dManIndex++ ) {
        	DofManager *dMan = domain->giveDofManager(dManIndex);

            for ( int k = 1; k <= dMan->giveNumberOfDofs(); k++ ) {

            	Dof *dof = dMan->giveDof(k);
            	int eqNumNew = dof->giveEqn();

                if ( eqNumNew > 0 ) {

                	std::vector<int> key(3);
                	key[0] = domainIndex;
                	key[1] = dManIndex;
                	key[2] = k;

                	if( mDofEqnNumMap.find(key) != mDofEqnNumMap.end() ){
                		int eqNumOld = mDofEqnNumMap[ key ];
//                		printf("eqNumNew: %d eqNumOld: %d\n", eqNumNew, eqNumOld);

                		if( eqNumOld > 0 && eqNumOld <= numEqOld  ) {
                			oArray.at(eqNumNew) = iArray.at(eqNumOld);
                		}
                	}
                }

            }
        }
    }


}


} /* namespace oofem */
