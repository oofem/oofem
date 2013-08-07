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

#include "patchintegrationrule.h"
#include "xfemelementinterface.h"
//#include "patch.h"
#include "integrationrule.h"
#include "gaussintegrationrule.h"
#include "geometry.h"
#include "classfactory.h"
#include "contextioerr.h"
#include "datastream.h"
#include "gausspoint.h"
#include "fei2dtrlin.h"

#include "XFEMDebugTools.h"


namespace oofem {
PatchIntegrationRule :: PatchIntegrationRule(int n, Element *e, const std::vector<Triangle> &iTriangles) :
GaussIntegrationRule(n, e),
mTriangles(iTriangles)
{

}

PatchIntegrationRule :: ~PatchIntegrationRule()
{

}

FEI2dTrLin PatchIntegrationRule :: mTriInterp(1, 2);

int
PatchIntegrationRule :: SetUpPointsOnTriangle(int nPoints, MaterialMode mode)
{
	int pointsPassed = 0;

	// TODO: set properly
    firstLocalStrainIndx = 1;
    lastLocalStrainIndx = 3;


    ////////////////////////////////////////////
    // Allocate Gauss point array
    int nPointsTot = nPoints*mTriangles.size();
    FloatArray coords_xi1, coords_xi2, weights;
    this->giveTriCoordsAndWeights(nPoints, coords_xi1, coords_xi2, weights);
    this->gaussPointArray = new GaussPoint * [ nPointsTot ];
    ////////////////////////////////////////////


    std::vector<FloatArray> newGPCoord;

	// Loop over triangles
	for(int i = 0; i < mTriangles.size(); i++)
	{
        // TODO: Probably unnecessary to allocate here
        const FloatArray **coords = new const FloatArray * [ mTriangles[i].giveNrVertices() ];
        // this we should put into the function before
        for ( int k = 0; k < mTriangles[i].giveNrVertices(); k++ ) {
            coords [ k ] = new FloatArray( *(mTriangles[i].giveVertex(k+1)) );
        }

		// Can not be used because it writes to the start of the array instead of appending.
//		int nPointsTri = GaussIntegrationRule :: SetUpPointsOnTriangle(nPoints, mode);

		for(int j = 0; j < nPoints; j++)
		{
		    FloatArray global;
	        GaussPoint *&gp = this->gaussPointArray[ pointsPassed ];

	        FloatArray *coord = new FloatArray(2);
	        coord->at(1) = coords_xi1.at(j+1);
	        coord->at(2) = coords_xi2.at(j+1);
	        gp = new GaussPoint(this, pointsPassed+1, coord, weights.at(j+1), mode);



	        mTriInterp.local2global(global, * gp->giveCoordinates(),
	                                         FEIVertexListGeometryWrapper(mTriangles[i].giveNrVertices(), coords));

	        newGPCoord.push_back(global);


	        FloatArray local;
	        this->elem->computeLocalCoordinates(local, global);

	        gp->setCoordinates(local);



	        double parentArea = this->elem->computeArea();

	        double refElArea = this->elem->giveParentElSize();

	        gp->setWeight(2.0 *refElArea* gp->giveWeight() * mTriangles[i].getArea() / parentArea); // update integration weight


	        pointsPassed++;
		}


        for ( int k = 0; k < mTriangles[i].giveNrVertices(); k++ ) {
            delete coords [ k ];
        }

        delete [] coords;

	}

#if PATCH_INT_DEBUG > 0
	int elIndex = this->elem->giveGlobalNumber();
    std::stringstream str;
    str << "GaussPoints" << elIndex << ".vtk";
    std::string name = str.str();

	XFEMDebugTools::WritePointsToVTK( name, newGPCoord);
#endif

	numberOfIntegrationPoints = pointsPassed;


	return numberOfIntegrationPoints;

}

contextIOResultType
PatchIntegrationRule :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
	// TODO: Implement

    //
    // saves full  context (saves state variables, that completely describe
    // current state)
    //

    // save parent data
    contextIOResultType iores;

    if ( ( iores = IntegrationRule :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
/*
    // save patch data
    if ( this->patch ) {
        // store patch type
        int _type = this->patch->givePatchType();
        if ( !stream->write(& _type, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        patch->saveContext(stream, mode, obj);
    } else {
        OOFEM_ERROR("saveContex : can't store NULL patch");
    }
*/
    return CIO_OK;
}

contextIOResultType
PatchIntegrationRule :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
	// TODO: Implement

    //
    // restores full element context (saves state variables, that completely describe
    // current state)
    //

    contextIOResultType iores;

    if ( stream == NULL ) {
        OOFEM_ERROR("restoreContex : can't write into NULL stream");
    }

    if ( ( iores = IntegrationRule :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
/*
    // restore patch data
    if ( this->patch ) {
        delete this->patch;
    }
*/
    int _ptype;
    if ( !stream->read(& _ptype, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
/*
    // create new patch
    this->patch = classFactory.createPatch( ( Patch :: PatchType ) _ptype, this->giveElement() );
    this->patch->restoreContext(stream, mode, obj);
*/
    return CIO_OK;
}
} // end namespace oofem

