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

#include "dofmanvalfield.h"
#include "domain.h"
#include "spatiallocalizer.h"
#include "element.h"
#include "crosssection.h"
#include "timestep.h"
#include "util.h"
#include "classfactory.h"
#include "dofmanager.h"
#include "feinterpol.h"
#include "problemmode.h"
#include "nodalrecoverymodel.h"
#include "vtkbaseexportmodule.h"
#include "octreelocalizer.h"


namespace oofem {
DofManValueField::DofManValueField(FieldType ft, Domain *d) : Field(ft), dmanvallist()
{
    int ndofman = d->giveNumberOfDofManagers();
    this->domain = d;
    this->dmanvallist.resize(ndofman);
}

DofManValueField::DofManValueField(FieldType ft, int nNodes, int nElements, const std::string engngModel = "transienttransport", const std::string domainDofsDefaults = "heattransfer") : Field(ft), dmanvallist()
{
    //ft FieldType::FT_Unknown
    this->eModel = classFactory.createEngngModel(engngModel.c_str(), 1, NULL);
    if ( !this->eModel ) {
        OOFEM_ERROR("Couldn't create EngngModel\n");
    }
    eModel->Instanciate_init();//create one domain
    this->domain = eModel->giveDomain(1);
    this->domain->resolveDomainDofsDefaults(domainDofsDefaults.c_str() );
    this->domain->dofManagerList.clear();
    this->domain->dofManagerList.resize(nNodes);
    dmanvallist.resize(nNodes);
    this->domain->elementList.clear();
    this->domain->elementList.resize(nElements);
    this->crossSect = classFactory.createCrossSection("emptycs", 1, this->domain);//create one dummy cross-section
    this->domain->crossSectionList.clear();
    this->domain->crossSectionList.resize(1);
    this->domain->crossSectionList[0] = std::move(this->crossSect);
    std::unique_ptr< SpatialLocalizer > spatialLocalizer;
    spatialLocalizer = std::make_unique<OctreeSpatialLocalizer>(this->domain);
    spatialLocalizer->init();    
    this->domain->setSpatialLocalizer(std::move(spatialLocalizer));
}


void DofManValueField::addNode(int i, const FloatArray &coords) {
    std::unique_ptr< DofManager >dman(classFactory.createDofManager("node", i, eModel->giveDomain(1) ) );
    if ( !dman ) {
        OOFEM_ERROR("Couldn't create node %d\n", i);
    }
    dman->setCoordinates(coords);
    dman->setGlobalNumber(i);
    this->domain->dofManagerList [ i - 1 ] = std::move(dman);
}


void 
DofManValueField::addElement(int i, const char *name, const IntArray &nodes) {
    std::unique_ptr< Element >elem(classFactory.createElement(name, i, this->domain) );
    if ( !elem ) {
        OOFEM_ERROR("Couldn't create element %d: %s\n", i, name);
    }

    elem->setDofManagers(nodes);
    elem->setGlobalNumber(i);
    elem->setCrossSection(1);
    this->domain->elementList [ i - 1 ] = std::move(elem);
}


int
DofManValueField::evaluateAt(FloatArray &answer, const FloatArray &coords, ValueModeType mode, TimeStep *tStep)
{
    int result = 0; // assume ok
    FloatArray lc, n;
    answer.resize(0);
    
    // request element containing target point
    Element *elem = this->domain->giveSpatialLocalizer()->giveElementContainingPoint(coords);
    if ( elem ) { // ok element containing target point found
        FEInterpolation *interp = elem->giveInterpolation();
        if ( interp ) {
            // map target point to element local coordinates
            if ( interp->global2local(lc, coords, FEIElementGeometryWrapper(elem) ) ) {
                // evaluate interpolation functions at target point
                interp->evalN(n, lc, FEIElementGeometryWrapper(elem) );
                // loop over element nodes
                for ( int i = 1; i <= n.giveSize(); i++ ) {
                    // multiply nodal value by value of corresponding shape function and add this to answer
                    answer.add(n.at(i), this->dmanvallist [ elem->giveDofManagerNumber(i) - 1 ]);
                }
            } else { // mapping from global to local coordinates failed
                OOFEM_ERROR("Error in mapping from global to local coordinates\n");
            }
        } else {  // element without interpolation
            OOFEM_ERROR("Error in element without interpolation\n");
        }
    } else { // no element containing given point found
        OOFEM_ERROR("Point out of defined elements\n");
    }
    return result;
}


int
DofManValueField::evaluateAt(FloatArray &answer, DofManager *dman, ValueModeType mode, TimeStep *tStep)
{
    answer = this->dmanvallist [ dman->giveNumber() - 1 ];
    return 1;
}

void
DofManValueField::setDofManValue(int dofMan, FloatArray value)
{
    this->dmanvallist [ dofMan - 1 ] = std::move(value);
}

const FloatArray &DofManValueField::getNodeCoordinates(int i)
{
    return this->domain->dofManagerList [ i - 1 ]->giveCoordinates();
}


void
DofManValueField::saveContext(DataStream &stream)
{}

void
DofManValueField::restoreContext(DataStream &stream)
{}
} // end namespace oofem
