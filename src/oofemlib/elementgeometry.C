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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#include "elementgeometry.h"
#include "crosssection.h"
#include "integrationrule.h"
#include "errorestimator.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "primaryfield.h"
#include "verbose.h"
#include "entityrenumberingscheme.h"
#include "error.h"
#include "classfactory.h"
#include "datastream.h"
#include "materialmapperinterface.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "feinterpol.h"
#include "feinterpol1d.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"
#include "function.h"
#include "dofmanager.h"
#include "node.h"
#include "dynamicinputrecord.h"
#include "matstatmapperint.h"

#include <cstdio>

namespace oofem {
ElementGeometry :: ElementGeometry(int n, Domain *aDomain) :
    FEMComponent(n, aDomain), dofManArray()
{
    material           = 0;
    numberOfDofMans    = 0;
    numberOfIntegrationRules = 0;
	activityTimeFunction = 0;
    integrationRulesArray  = NULL;
}


ElementGeometry :: ~ElementGeometry()
{
    if ( integrationRulesArray ) {
        for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
            delete integrationRulesArray [ i ];
        }

        delete[] integrationRulesArray;
    }
}

Material *ElementGeometry :: giveMaterial()
// Returns the material of the receiver.
{
#ifdef DEBUG
    if ( !material ) {
        // material = this -> readInteger("mat") ;
        OOFEM_ERROR("giveMaterial: material not defined");
    }
#endif
    return domain->giveMaterial(material);
}


CrossSection *ElementGeometry :: giveCrossSection()
// Returns the crossSection of the receiver.
{
#ifdef DEBUG
    if ( !crossSection ) {
      OOFEM_ERROR("giveCrossSection: crossSection not defined");
    }
#endif
    return domain->giveCrossSection(crossSection);
}


int
ElementGeometry :: giveRegionNumber()
{
    return this->giveCrossSection()->giveNumber();
}


DofManager *
ElementGeometry :: giveDofManager(int i) const
// Returns the i-th node of the receiver.
{
#ifdef DEBUG
    if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
        OOFEM_ERROR("giveNode: Node %i is not defined", i);
    }
#endif
    return domain->giveDofManager( dofManArray.at(i) );
}


Node *
ElementGeometry :: giveNode(int i) const
// Returns the i-th node of the receiver.
{
#ifdef DEBUG
    if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
        OOFEM_ERROR("giveNode: Node is not defined");
    }
#endif
    return domain->giveNode( dofManArray.at(i) );
}


ElementSide *
ElementGeometry :: giveSide(int i) const
// Returns the i-th side of the receiver.
{
#ifdef DEBUG
    if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
        OOFEM_ERROR("giveNode: Side is not defined");
    }
#endif
    return domain->giveSide( dofManArray.at(i) );
}


void
ElementGeometry :: setDofManagers(const IntArray &_dmans)
{
    this->dofManArray = _dmans;
}


void
ElementGeometry :: setIntegrationRules(AList< IntegrationRule > *irlist)
{
    if ( integrationRulesArray ) {
        for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
            delete integrationRulesArray [ i ];
        }

        delete[] integrationRulesArray;
    }

    numberOfIntegrationRules = irlist->giveSize();
    integrationRulesArray = new IntegrationRule * [ irlist->giveSize() ];

    for ( int j = 0; j < irlist->giveSize(); j++ ) {
        integrationRulesArray [ j ] =  irlist->at(j + 1);
        irlist->unlink(j + 1);
    }
}





IRResultType
ElementGeometry :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                          // Required by IR_GIVE_FIELD macro

#  ifdef VERBOSE
    // VERBOSE_PRINT1("Instanciating element ",number);
#  endif
    //IR_GIVE_FIELD(ir, material, _IFT_Element_mat);
    IR_GIVE_OPTIONAL_FIELD(ir, material, _IFT_ElementGeometry_mat);

    //IR_GIVE_FIELD(ir, crossSection, _IFT_Element_crosssect);
	IR_GIVE_OPTIONAL_FIELD(ir, crossSection, _IFT_ElementGeometry_crosssect);

	IR_GIVE_FIELD(ir, dofManArray, _IFT_ElementGeometry_nodes);

    elemLocalCS.resize(0, 0);

	if (ir->hasField(_IFT_ElementGeometry_lcs)) { //local coordinate system
        double n1 = 0.0, n2 = 0.0;
        FloatArray triplets;
        triplets.resize(0);
		IR_GIVE_OPTIONAL_FIELD(ir, triplets, _IFT_ElementGeometry_lcs);
        elemLocalCS.resize(3, 3);
        for ( int j = 1; j <= 3; j++ ) {
            elemLocalCS.at(j, 1) = triplets.at(j);
            n1 += triplets.at(j) * triplets.at(j);
            elemLocalCS.at(j, 2) = triplets.at(j + 3);
            n2 += triplets.at(j + 3) * triplets.at(j + 3);
        }

        n1 = sqrt(n1);
        n2 = sqrt(n2);
        for ( int j = 1; j <= 3; j++ ) { // normalize e1' e2'
            elemLocalCS.at(j, 1) /= n1;
            elemLocalCS.at(j, 2) /= n2;
        }

        // vector e3' computed from vector product of e1', e2'
        elemLocalCS.at(1, 3) = ( elemLocalCS.at(2, 1) * elemLocalCS.at(3, 2) - elemLocalCS.at(3, 1) * elemLocalCS.at(2, 2) );
        elemLocalCS.at(2, 3) = ( elemLocalCS.at(3, 1) * elemLocalCS.at(1, 2) - elemLocalCS.at(1, 1) * elemLocalCS.at(3, 2) );
        elemLocalCS.at(3, 3) = ( elemLocalCS.at(1, 1) * elemLocalCS.at(2, 2) - elemLocalCS.at(2, 1) * elemLocalCS.at(1, 2) );
        //elemLocalCS.printYourself();
    }

#ifdef __PARALLEL_MODE
    partitions.resize(0);
	IR_GIVE_OPTIONAL_FIELD(ir, partitions, _IFT_ElementGeometry_partitions);
    if ( ir->hasField(_IFT_Element_remote) ) {
        parallel_mode = Element_remote;
    } else {
        parallel_mode = Element_local;
    }

#endif

 

	IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, _IFT_ElementGeometry_nip);

	IR_GIVE_OPTIONAL_FIELD(ir, activityTimeFunction, _IFT_ElementGeometry_activityTimeFunction);


    return IRRT_OK;
}


void
ElementGeometry :: giveInputRecord(DynamicInputRecord &input)
{
    FEMComponent :: giveInputRecord(input);

	input.setField(material, _IFT_ElementGeometry_mat);

	input.setField(crossSection, _IFT_ElementGeometry_crosssect);

	input.setField(dofManArray, _IFT_ElementGeometry_nodes);



    if ( elemLocalCS.giveNumberOfRows() > 0 ) {
        FloatArray triplets(6);
        for ( int j = 1; j <= 3; j++ ) {
            triplets.at(j) = elemLocalCS.at(j, 1);
            triplets.at(j + 3) = elemLocalCS.at(j, 2);
        }
		input.setField(triplets, _IFT_ElementGeometry_lcs);
    }


#ifdef __PARALLEL_MODE
    if ( this->giveDomain()->giveEngngModel()->isParallel() && partitions.giveSize() > 0 ) {
		input.setField(this->partitions, _IFT_ElementGeometry_partitions);
		if ( this->parallel_mode == ElementGeometry_remote ) {
			input.setField(_IFT_ElementGeometry_remote);
        }
    }
#endif


	input.setField(numberOfGaussPoints, _IFT_ElementGeometry_nip);

	if (activityTimeFunction > 0) {
		input.setField(activityTimeFunction, _IFT_ElementGeometry_activityTimeFunction);
	}


}


void
ElementGeometry :: postInitialize()
{
    this->computeGaussPoints();
}


void
ElementGeometry :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->printOutputAt(file, tStep);
    }
}


void
ElementGeometry :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
#  ifdef VERBOSE
    // VERBOSE_PRINT1("Updating Element ",number)
#  endif

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->updateYourself(tStep);
    }
}



void
ElementGeometry :: initForNewStep()
// initializes receiver to new time step or can be used
// if current time step must be restarted
{
    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->initForNewStep();
    }
}


contextIOResultType ElementGeometry :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    int _val;

    if ( ( iores = FEMComponent :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( mode & CM_Definition ) ) {
        if ( !stream->write(& numberOfDofMans, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->write(& material, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->write(& crossSection, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

#ifdef __PARALLEL_MODE
        if ( mode & CM_DefinitionGlobal ) {
            // send global numbers instead of local ones
            int s = dofManArray.giveSize();
            IntArray globDN(s);
            for ( int i = 1; i <= s; i++ ) {
                globDN.at(i) = this->giveDofManager(i)->giveGlobalNumber();
            }

            if ( ( iores = globDN.storeYourself(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        } else {
            if ( ( iores = dofManArray.storeYourself(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        }

#else
        if ( ( iores = dofManArray.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

#endif
       

        if ( !stream->write(& numberOfIntegrationRules, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
            _val = integrationRulesArray [ i ]->giveIntegrationRuleType();
            if ( !stream->write(& _val, 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

#ifdef __PARALLEL_MODE
        int _mode;
        if ( !stream->write(& globalNumber, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        _mode = parallel_mode;
        if ( !stream->write(& _mode, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = partitions.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

#endif
    }

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        if ( ( iores = integrationRulesArray [ i ]->saveContext(stream, mode, obj) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType ElementGeometry :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    int _nrules;

    if ( ( iores = FEMComponent :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        if ( !stream->read(& numberOfDofMans, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->read(& material, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->read(& crossSection, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = dofManArray.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream->read(& _nrules, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        // restore integration rules
        IntArray dtypes(_nrules);
        for ( int i = 1; i <= _nrules; i++ ) {
            if ( !stream->read(& dtypes.at(i), 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        if ( _nrules != numberOfIntegrationRules ) {
            // delete old int rule array
            if ( integrationRulesArray ) {
                for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
                    delete integrationRulesArray [ i ];
                }

                delete[] integrationRulesArray;
            }

            // AND ALLOCATE NEW ONE
            integrationRulesArray = new IntegrationRule * [ _nrules ];
            for ( int i = 0; i < _nrules; i++ ) {
                integrationRulesArray [ i ] = classFactory.createIRule( ( IntegrationRuleType ) dtypes(i), i + 1, this );
            }

            numberOfIntegrationRules = _nrules;
        } else {
            for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
                if ( integrationRulesArray [ i ]->giveIntegrationRuleType() != dtypes(i) ) {
                    delete integrationRulesArray [ i ];
                    integrationRulesArray [ i ] = classFactory.createIRule( ( IntegrationRuleType ) dtypes(i), i + 1, this );
                }
            }
        }

#ifdef __PARALLEL_MODE
        int _mode;
        if ( !stream->read(& globalNumber, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->read(& _mode, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        parallel_mode = ( elementParallelMode ) _mode;
        if ( ( iores = partitions.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

#endif
    }


    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        if ( ( iores = integrationRulesArray [ i ]->restoreContext(stream, mode, this) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


double
ElementGeometry :: computeVolumeAreaOrLength()
// the element computes its volume, area or length
// (depending on the spatial dimension of that element)
{
    GaussPoint *gp;
    double answer = 0.;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    if ( iRule ) {
        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule->getIntegrationPoint(i);
            answer += this->computeVolumeAround(gp);
        }

        return answer;
    }

    return -1.; // means "cannot be evaluated"
}


double
ElementGeometry :: computeMeanSize()
// Computes the size of the element defined as its length,
// square root of area or cube root of volume (depending on spatial dimension)
{
    // if the method "giveArea" is properly implemented
    // for the particular type of element, the mean size is the square root of area
    // 8 July 2010 - does not seem to work any more (Element does not inherit from ElementGeometry)
    // double area = this->giveArea();
    // if (area>0.)
    //  return sqrt(area);

    // if "giveArea" is not implemented (default value 0.),
    // then the contributing areas or volumes are collected from Gauss points
    double volume = this->computeVolumeAreaOrLength();
    if ( volume < 0. ) {
        return -1.; // means "cannot be evaluated"
    }

    int dim = this->giveSpatialDimension();
    switch ( dim ) {
    case 1: return volume;

    case 2: return sqrt(volume);

    case 3: return pow(volume, 1. / 3.);
    }

    return -1.; // means "cannot be evaluated"
}


double
ElementGeometry :: computeVolume()
{
    FEInterpolation3d *fei = dynamic_cast< FEInterpolation3d * >( this->giveInterpolation() );
#ifdef DEBUG
    if ( !fei ) {
        OOFEM_ERROR("Element :: computeVolume - Function not overloaded and necessary interpolator isn't available");
        return 0.0;
    }
#endif
    return fei->giveVolume( FEIElementGeometryWrapper(this) );
}


double
ElementGeometry :: computeArea()
{
    FEInterpolation2d *fei = dynamic_cast< FEInterpolation2d * >( this->giveInterpolation() );
#ifdef DEBUG
    if ( !fei ) {
        OOFEM_ERROR("Element :: computeArea - Function not overloaded and necessary interpolator isn't available");
        return 0.0;
    }
#endif
    return fei->giveArea( FEIElementGeometryWrapper(this) );
}


double
ElementGeometry :: computeLength()
{
    FEInterpolation1d *fei = dynamic_cast< FEInterpolation1d * >( this->giveInterpolation() );
#ifdef DEBUG
    if ( !fei ) {
        OOFEM_ERROR("Element :: computeLength - Function not overloaded and necessary interpolator isn't available");
        return 0.0;
    }
#endif
    return fei->giveLength( FEIElementGeometryWrapper(this) );
}


double
ElementGeometry :: giveLenghtInDir(const FloatArray &normalToCrackPlane)
//
// returns receivers projection length (for some material models)
// to direction given by normalToCrackPlane;
//
{
    FloatArray *coords;
    double maxDis, minDis, dis;
    int nnode = giveNumberOfNodes();

    coords = this->giveNode(1)->giveCoordinates();
    minDis = maxDis = normalToCrackPlane.dotProduct( * coords, coords->giveSize() );

    for ( int i = 2; i <= nnode; i++ ) {
        coords = this->giveNode(i)->giveCoordinates();
        dis = normalToCrackPlane.dotProduct( * coords, coords->giveSize() );
        if ( dis > maxDis ) {
            maxDis = dis;
        } else if ( dis < minDis ) {
            minDis = dis;
        }
    }

    return maxDis - minDis;
}


int
ElementGeometry :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    FEInterpolation *fei = this->giveInterpolation();
#ifdef DEBUG
    if ( !fei ) {
        answer.resize(0);
        return false;
    }
#endif
    fei->local2global( answer, lcoords, FEIElementGeometryWrapper(this) );
    return true;
}


bool
ElementGeometry :: computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords)
{
    FEInterpolation *fei = this->giveInterpolation();
    if ( fei ) {
        return fei->global2local( answer, gcoords, FEIElementGeometryWrapper(this) );
    } else {
        return false;
    }
}


int
ElementGeometry :: giveLocalCoordinateSystem(FloatMatrix &answer)
{
    if ( elemLocalCS.isNotEmpty() ) {
        answer = elemLocalCS;
        return 1;
    } else {
        answer.clear();
    }

    return 0;
}


void
ElementGeometry :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *)
// valid only for plane elements (shells, plates, ....)
// computes mid-plane normal at gaussPoint - for materials with orthotrophy
{
    OOFEM_ERROR("Unable to compute mid-plane normal, not supported");
}


int
ElementGeometry :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_ErrorIndicatorLevel ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            answer.at(1) = ee->giveElementError(indicatorET, this, tStep);
        } else {
            answer.resize(0);
            return 0;
        }

        return 1;
    } else if ( type == IST_InternalStressError ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            answer.at(1) = ee->giveElementError(internalStressET, this, tStep);
        } else {
            answer.resize(0);
            return 0;
        }
        return 1;
    } else if ( type == IST_PrimaryUnknownError ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            answer.at(1) = ee->giveElementError(primaryUnknownET, this, tStep);
        } else {
            answer.resize(0);
            return 0;
        }
        return 1;
    } else if ( type == IST_CrossSectionNumber ) {
        answer.resize(1);
        answer.at(1) = gp->giveCrossSection()->giveNumber();
        return 1;
    } else if ( type == IST_ElementNumber ) {
        answer.resize(1);
        answer.at(1) = this->giveNumber();
        return 1;
    } else {
        return this->giveCrossSection()->giveIPValue(answer, gp, type, tStep);
    }
}

int
ElementGeometry :: giveSpatialDimension()
{
    ///@todo Just ask the interpolator instead?
    switch ( this->giveGeometryType() ) {
    case EGT_point:
        return 0;

    case EGT_line_1:
    case EGT_line_2:
        return 1;

    case EGT_triangle_1:
    case EGT_triangle_2:
    case EGT_quad_1:
    case EGT_quad_2:
    case EGT_quad9_2:
    case EGT_quad_1_interface:
    case EGT_quad_21_interface:
        return 2;

    case EGT_tetra_1:
    case EGT_tetra_2:
    case EGT_hexa_1:
    case EGT_hexa_2:
    case EGT_hexa_27:
    case EGT_wedge_1:
    case EGT_wedge_2:
        return 3;

    case EGT_Composite:
    case EGT_unknown:
        break;
    }

    OOFEM_ERROR("giveSpatialDimension: failure (maybe new element type was registered)");
    return 0; //to make compiler happy
}


int
ElementGeometry :: giveNumberOfBoundarySides()
{
    ///@todo Just ask the interpolator instead?
    switch ( this->giveGeometryType() ) {
    case EGT_point:
        return 0;

    case EGT_line_1:
    case EGT_line_2:
    case EGT_quad_1_interface:
    case EGT_quad_21_interface:
        return 2;

    case EGT_triangle_1:
    case EGT_triangle_2:
        return 3;

    case EGT_quad_1:
    case EGT_quad_2:
    case EGT_quad9_2:
        return 4;

    case EGT_tetra_1:
    case EGT_tetra_2:
        return 4;

    case EGT_wedge_1:
    case EGT_wedge_2:
        return 5;

    case EGT_hexa_1:
    case EGT_hexa_2:
    case EGT_hexa_27:
        return 6;

    case EGT_Composite:
    case EGT_unknown:
        break;
    }

    OOFEM_ERROR( "giveSpatialDimension: failure, unsupported geometry type (%s)",
             __Element_Geometry_TypeToString( this->giveGeometryType() ) );
    return 0; // to make compiler happy
}


int
ElementGeometry :: adaptiveMap(Domain *oldd, TimeStep *tStep)
{
    int result = 1;
    IntegrationRule *iRule;
    MaterialModelMapperInterface *interface = static_cast< MaterialModelMapperInterface * >
                                              ( this->giveMaterial()->giveInterface(MaterialModelMapperInterfaceType) );

    if ( !interface ) {
        return 0;
    }

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            result &= interface->MMI_map(iRule->getIntegrationPoint(j), oldd, tStep);
        }
    }

    return result;
}

int
ElementGeometry :: mapStateVariables(Domain &iOldDom, const TimeStep &iTStep)
{
    int result = 1;
    // create source set (this is quite inefficient here as done for each element geometry.
    // the alternative MaterialModelMapperInterface approach allows to cache sets on material model
    Set sourceElemSet = Set(0, & iOldDom);
    IntArray el;
    // compile source list to contain all elements on old odmain with the same material id
    for ( int i = 1; i <= iOldDom.giveNumberOfElements(); i++ ) {
      if ( iOldDom.giveElementGeometry(i)->giveMaterial()->giveNumber() == this->giveMaterial()->giveNumber() ) {
	// add oldd domain element to source list
	el.followedBy(i, 10);
      }
    }
    sourceElemSet.setElementList(el);

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        IntegrationRule *iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint &gp = * ( iRule->getIntegrationPoint(j) );

            MaterialStatus *ms = dynamic_cast< MaterialStatus * >( gp.giveMaterialStatus() );
            if ( ms == NULL ) {
                OOFEM_ERROR("In Element :: mapStateVariables(): failed to fetch MaterialStatus.\n");
            }

            MaterialStatusMapperInterface *interface = dynamic_cast< MaterialStatusMapperInterface * >( ms );
            if ( interface == NULL ) {
                OOFEM_ERROR("In Element :: mapStateVariables(): Failed to fetch MaterialStatusMapperInterface.\n");
            }

            result &= interface->MSMI_map( gp, iOldDom, sourceElemSet, iTStep, * ( ms ) );
        }
    }

    return result;
}


int
ElementGeometry :: adaptiveFinish(TimeStep *tStep)
{
    int result = 1;
    IntegrationRule *iRule;
    MaterialModelMapperInterface *interface = static_cast< MaterialModelMapperInterface * >
                                              ( this->giveMaterial()->giveInterface(MaterialModelMapperInterfaceType) );

    if ( !interface ) {
        return 0;
    }

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            result &= interface->MMI_finish(tStep);
        }
    }

    return result;
}


void
ElementGeometry :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        dofManArray.at(i) = f(dofManArray.at(i), ERS_DofManager);
    }
}


integrationDomain
ElementGeometry :: giveIntegrationDomain() const
{
    FEInterpolation *fei = this->giveInterpolation();
    return fei ? fei->giveIntegrationDomain() : _Unknown_integrationDomain;
}


Element_Geometry_Type
ElementGeometry :: giveGeometryType() const
{
    FEInterpolation *fei = this->giveInterpolation();
    return fei ? fei->giveGeometryType() : EGT_unknown;
}


bool
ElementGeometry :: isActivated(TimeStep *tStep)
{
	if ( activityTimeFunction ) {
		if ( tStep ) {
			return ( domain->giveFunction(activityTimeFunction)->evaluateAtTime( tStep->giveIntrinsicTime() ) > 1.e-3 );
		} else {
			return false;
		}
	} else {
		return true;
	}
}


#ifdef __PARALLEL_MODE
int
ElementGeometry :: packUnknowns(CommunicationBuffer &buff, TimeStep *tStep)
{
    int result = 1;
    IntegrationRule *iRule;

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            result &= this->giveCrossSection()->packUnknowns( buff, tStep, iRule->getIntegrationPoint(j) );
        }
    }

    return result;
}


int
ElementGeometry :: unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *tStep)
{
    int result = 1;
    IntegrationRule *iRule;

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            result &= this->giveCrossSection()->unpackAndUpdateUnknowns( buff, tStep, iRule->getIntegrationPoint(j) );
        }
    }

    return result;
}


int
ElementGeometry :: estimatePackSize(CommunicationBuffer &buff)
{
    int result = 0;
    IntegrationRule *iRule;

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            result += this->giveCrossSection()->estimatePackSize( buff, iRule->getIntegrationPoint(j) );
        }
    }

    return result;
}


double
ElementGeometry :: predictRelativeComputationalCost()
{
    int nip;
    double wgt = 0;
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    nip = iRule->giveNumberOfIntegrationPoints();
    for ( int j = 0; j < nip; j++ ) {
        wgt += this->giveCrossSection()->predictRelativeComputationalCost( iRule->getIntegrationPoint(j) );
    }

    return ( this->giveRelativeSelfComputationalCost() * wgt );
}
#endif


#ifdef __OOFEG
void
ElementGeometry :: drawYourself(oofegGraphicContext &gc)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc);
    } else if ( mode == OGC_elementAnnotation ) {
        this->drawAnnotation(gc);
    } else if ( mode == OGC_deformedGeometry ) {
        this->drawDeformedGeometry(gc, DisplacementVector);
    } else if ( mode == OGC_eigenVectorGeometry ) {
        this->drawDeformedGeometry(gc, EigenVector);
    } else if ( mode == OGC_scalarPlot ) {
        this->drawScalar(gc);
    } else if ( mode == OGC_elemSpecial ) {
        this->drawSpecial(gc);
    } else {
        _error("drawYourself : unsupported mode");
    }
}


void
ElementGeometry :: drawAnnotation(oofegGraphicContext &gc)
{
    int i, count = 0;
    Node *node;
    WCRec p [ 1 ]; /* point */
    GraphicObj *go;
    char num [ 30 ];

    p [ 0 ].x = p [ 0 ].y = p [ 0 ].z = 0.0;
    // compute element center
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        if ( ( node = this->giveNode(i) ) ) {
            p [ 0 ].x += node->giveCoordinate(1);
            p [ 0 ].y += node->giveCoordinate(2);
            p [ 0 ].z += node->giveCoordinate(3);
            count++;
        }
    }

    p [ 0 ].x /= count;
    p [ 0 ].y /= count;
    p [ 0 ].z /= count;

    EASValsSetLayer(OOFEG_ELEMENT_ANNOTATION_LAYER);
    EASValsSetColor( gc.getElementColor() );
 #ifdef __PARALLEL_MODE
    sprintf( num, "%d(%d)", this->giveNumber(), this->giveGlobalNumber() );
 #else
    sprintf( num, "%d", this->giveNumber() );
 #endif
    go = CreateAnnText3D(p, num);
    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


int
ElementGeometry :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                   int node, TimeStep *tStep)
{
    if ( type == IST_RelMeshDensity ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            answer.at(1) = this->giveDomain()->giveErrorEstimator()->giveRemeshingCrit()->
                           giveRequiredDofManDensity(this->giveNode(node)->giveNumber(), tStep, 1);
            return 1;
        } else {
            answer.resize(0);
            return 0;
        }
    } else {
        if ( mode == ISM_recovered ) {
            const FloatArray *nodval;
            NodalRecoveryModel *smoother = this->giveDomain()->giveSmoother();
            int result = smoother->giveNodalVector( nodval, this->giveNode(node)->giveNumber(),
                                                    smoother->giveElementVirtualRegionNumber(this->number) );
            if ( nodval ) {
                answer = * nodval;
            } else {
                answer.resize(0);
            }

            return result;
        } else {
            return 0;
        }
    }
}


#endif
} // end namespace oofem
