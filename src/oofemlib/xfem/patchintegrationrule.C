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
#include "fei3dtrquad.h"
#include "XFEMDebugTools.h"

#include "timestep.h"
#include "engngm.h"

namespace oofem {
PatchIntegrationRule :: PatchIntegrationRule(int n, Element *e, const std :: vector< Triangle > &iTriangles) :
    GaussIntegrationRule(n, e),
    mTriangles(iTriangles)
{ }

PatchIntegrationRule :: ~PatchIntegrationRule()
{ }

FEI2dTrLin PatchIntegrationRule :: mTriInterp(1, 2);
FEI3dTrQuad PatchIntegrationRule :: mTriInterpQuad;

int
PatchIntegrationRule :: SetUpPointsOnTriangle(int nPoints, MaterialMode mode)
{
    int pointsPassed = 0;

    // TODO: set properly
    firstLocalStrainIndx = 1;
    lastLocalStrainIndx = 3;


    ////////////////////////////////////////////
    // Allocate Gauss point array


    // It may happen that the patch contains triangles with
    // zero area. This does no harm, since their weights in
    // the quadrature will be zero. However, they invoke additional
    // computational cost and therefore we want to avoid them.
    // Thus, count the number of triangles with finite area
    // and keep only those triangles.

    double totArea = 0.0;
    for ( size_t i = 0; i < mTriangles.size(); i++ ) {
        totArea += mTriangles [ i ].getArea();
    }

    std :: vector< int >triToKeep;
    const double triTol = ( 1.0e-6 ) * totArea;

    for ( size_t i = 0; i < mTriangles.size(); i++ ) {
        if ( mTriangles [ i ].getArea() > triTol ) {
            triToKeep.push_back(i);
        }
    }

    int nPointsTot = nPoints * triToKeep.size();
    FloatArray coords_xi1, coords_xi2, weights;
    this->giveTriCoordsAndWeights(nPoints, coords_xi1, coords_xi2, weights);
    this->gaussPoints.resize( nPointsTot );
    ////////////////////////////////////////////


    std :: vector< FloatArray >newGPCoord;

    double parentArea = this->elem->computeArea();

    // Loop over triangles
    for ( int tri: triToKeep ) {
        // TODO: Probably unnecessary to allocate here
        std::vector< FloatArray > coords( mTriangles [ tri ].giveNrVertices() );
        // this we should put into the function before
        for ( int k = 1; k <= mTriangles [ tri ].giveNrVertices(); k++ ) {
            coords[ k - 1 ] = mTriangles [ tri ].giveVertex( k );
        }

        // Can not be used because it writes to the start of the array instead of appending.
        //		int nPointsTri = GaussIntegrationRule :: SetUpPointsOnTriangle(nPoints, mode);

        for ( int j = 0; j < nPoints; j++ ) {
            FloatArray global;
            GaussPoint * &gp = this->gaussPoints [ pointsPassed ];

            FloatArray *coord = new FloatArray(2);
            coord->at(1) = coords_xi1.at(j + 1);
            coord->at(2) = coords_xi2.at(j + 1);
            gp = new GaussPoint(this, pointsPassed + 1, coord, weights.at(j + 1), mode);



            mTriInterp.local2global( global, * gp->giveNaturalCoordinates(),
                                    FEIVertexListGeometryWrapper(coords) );

            newGPCoord.push_back(global);


            FloatArray local;
            this->elem->computeLocalCoordinates(local, global);

            gp->setGlobalCoordinates(global);
            gp->setNaturalCoordinates(local);
            gp->setSubPatchCoordinates(local);




            double refElArea = this->elem->giveParentElSize();

            gp->setWeight(2.0 * refElArea * gp->giveWeight() * mTriangles [ tri ].getArea() / parentArea); // update integration weight


            pointsPassed++;
        }
    }

    XfemManager *xMan = elem->giveDomain()->giveXfemManager();
    if ( xMan != NULL ) {
        if ( xMan->giveVtkDebug() ) {
            double time = 0.0;

            Element *el = this->elem;
            if ( el != NULL ) {
                Domain *dom = el->giveDomain();
                if ( dom != NULL ) {
                    EngngModel *em = dom->giveEngngModel();
                    if ( em != NULL ) {
                        TimeStep *ts = em->giveCurrentStep();
                        if ( ts != NULL ) {
                            time = ts->giveTargetTime();
                        }
                    }
                }
            }

            int elIndex = this->elem->giveGlobalNumber();
            std :: stringstream str;
            str << "GaussPointsTime" << time << "El" << elIndex << ".vtk";
            std :: string name = str.str();

            XFEMDebugTools :: WritePointsToVTK(name, newGPCoord);
        }
    }

    return this->giveNumberOfIntegrationPoints();
}




int
PatchIntegrationRule :: SetUpPointsOnWedge(int nPointsTri, int nPointsDepth, MaterialMode mode)
{
    int pointsPassed = 0;

    // TODO: set properly
    firstLocalStrainIndx = 1;
    lastLocalStrainIndx = 3;

    double totArea = 0.0;
    for ( size_t i = 0; i < mTriangles.size(); i++ ) {
        totArea += mTriangles [ i ].getArea();
    }

    std :: vector< int >triToKeep;
    const double triTol = ( 1.0e-6 ) * totArea;

    for ( size_t i = 0; i < mTriangles.size(); i++ ) {
        if ( mTriangles [ i ].getArea() > triTol ) {
            triToKeep.push_back(i);
        }
    }

    int nPointsTot = nPointsTri * nPointsDepth * triToKeep.size();
    FloatArray coords_xi1, coords_xi2, coords_xi3, weightsTri, weightsDepth;
    this->giveTriCoordsAndWeights(nPointsTri, coords_xi1, coords_xi2, weightsTri);
    this->giveLineCoordsAndWeights(nPointsDepth, coords_xi3, weightsDepth);

    //this->gaussPointArray = new GaussPoint * [ nPointsTot ];
    this->gaussPoints.resize(nPointsTot);

    std :: vector< FloatArray >newGPCoord;

    double parentArea = this->elem->computeArea();
    int count = 0;

    // Loop over triangles
    for ( int i = 0; i < int( triToKeep.size() ); i++ ) {

        Triangle triangle = mTriangles [ triToKeep [ i ] ];
        //const FloatArray **gCoords = new const FloatArray * [ triangle.giveNrVertices() ];
        std::vector< FloatArray > gCoords( triangle.giveNrVertices() );

        // global coords of the the triangle verticies
        for ( int j = 0; j < triangle.giveNrVertices(); j++ ) {
            //gCoords [ j ] = new FloatArray( ( triangle.giveVertex(j + 1) ) );
            gCoords[j] = (triangle.giveVertex(j + 1));
        }
        

        for ( int k = 1; k <= nPointsTri; k++ ) {
            for ( int m = 1; m <= nPointsDepth; m++ ) {
                FloatArray *lCoords = new FloatArray(3);
                lCoords->at(1) = coords_xi1.at(k);
                lCoords->at(2) = coords_xi2.at(k);
                lCoords->at(3) = coords_xi3.at(m);

                // Compute global gp coordinate in the element from local gp coord in patch triangle
                FloatArray global;
                //mTriInterp.local2global( global, *lCoords, FEIVertexListGeometryWrapper( triangle.giveNrVertices(), gCoords) );
                mTriInterp.local2global(global, *lCoords, FEIVertexListGeometryWrapper(gCoords));
                

                // Compute local gp coordinate in the element from global gp coord in the element
                FloatArray local;
                this->elem->computeLocalCoordinates(local, global);
                local.at(3) = coords_xi3.at(m); // manually set third coordinate
    
                double refElArea = 0.5;
                //double refElArea = this->elem->giveParentElSize();
                double oldWeight = weightsTri.at(k) * weightsDepth.at(m);
                double newWeight = 2.0 * refElArea * oldWeight * triangle.getArea() / parentArea; 


                FloatArray *coord = new FloatArray(3);
                coord->at(1) = local.at(1);
                coord->at(2) = local.at(2);
                coord->at(3) = local.at(3);

                GaussPoint *gp = new GaussPoint(this, count + 1, coord, newWeight, mode);

                //GaussPoint *gp = new GaussPoint(this, count + 1, &local, newWeight, mode);
                //this->gaussPointArray [ count ] = gp;
                this->gaussPoints[count] = gp;
                count++;

                // Store new global gp coord for vtk output
                this->elem->computeGlobalCoordinates(global, *gp->giveNaturalCoordinates() );
                newGPCoord.push_back(global);
            }
        }


        //for ( int k = 0; k < mTriangles [ triToKeep [ i ] ].giveNrVertices(); k++ ) {
        //    delete gCoords [ k ];
        //}

        //delete [] gCoords;
    }

    XfemManager *xMan = elem->giveDomain()->giveXfemManager();
    if ( xMan != NULL ) {
        if ( xMan->giveVtkDebug() ) {
            double time = 0.0;

            Element *el = this->elem;
            if ( el != NULL ) {
                Domain *dom = el->giveDomain();
                if ( dom != NULL ) {
                    EngngModel *em = dom->giveEngngModel();
                    if ( em != NULL ) {
                        TimeStep *ts = em->giveCurrentStep();
                        if ( ts != NULL ) {
                            time = ts->giveTargetTime();
                        }
                    }
                }
            }

            int elIndex = this->elem->giveGlobalNumber();
            std :: stringstream str;
            str << "GaussPointsTime" << time << "El" << elIndex << ".vtk";
            std :: string name = str.str();

            XFEMDebugTools :: WritePointsToVTK(name, newGPCoord);
        }
    }

    
    return this->giveNumberOfIntegrationPoints();
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
     *  // save patch data
     *  if ( this->patch ) {
     *      // store patch type
     *      int _type = this->patch->givePatchType();
     *      if ( !stream->write(& _type, 1) ) {
     *          THROW_CIOERR(CIO_IOERR);
     *      }
     *
     *      patch->saveContext(stream, mode, obj);
     *  } else {
     *      OOFEM_ERROR("can't store NULL patch");
     *  }
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
        OOFEM_ERROR("can't write into NULL stream");
    }

    if ( ( iores = IntegrationRule :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    /*
     *  // restore patch data
     *  if ( this->patch ) {
     *      delete this->patch;
     *  }
     */
    int _ptype;
    if ( !stream->read(& _ptype, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    /*
     *  // create new patch
     *  this->patch = classFactory.createPatch( ( Patch :: PatchType ) _ptype, this->giveElement() );
     *  this->patch->restoreContext(stream, mode, obj);
     */
    return CIO_OK;
}
} // end namespace oofem
