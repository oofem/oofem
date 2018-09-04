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

#include "delaunaytriangulator.h"
#include "octreelocalizert.h"
#include "delaunaytriangle.h"
#include "tr1_2d_pfem.h"
#include "intarray.h"
#include "contextioerr.h"
#include "verbose.h"
#include "timer.h"
#include "mathfem.h"
#include "pfemparticle.h"
#include "pfem.h"

#define _USING_OCTREE

namespace oofem {
DelaunayTriangulator :: DelaunayTriangulator(Domain *d, double setAlpha) :
    alphaShapeEdgeList(0), triangleOctree()
{
    domain = d;
    alphaValue = setAlpha;
    nnode = domain->giveNumberOfDofManagers();
    // Option 2: setting bounds vor computed Alpha
    //minAlpha = 0.08;
    //maxAlpha = 0.2;
}

DelaunayTriangulator :: ~DelaunayTriangulator()
{
    for ( auto tri : generalTriangleList ) {
        delete tri;
    }

    for ( auto edge : edgeList ) {
        delete edge;
    }
}



void DelaunayTriangulator :: addUniqueEdgeToPolygon(Edge2D edge, std :: list< Edge2D > &polygon)
{
    int do_add = 1;

    for ( auto pos = polygon.begin(); pos != polygon.end(); ) {
        if ( ( * pos ) == edge ) {
            pos = polygon.erase(pos);
            do_add = 0;
            break;
        } else {
            ++pos;
        }
    }

    if ( do_add ) {
        polygon.push_back(edge);
    }
}


void DelaunayTriangulator :: generateMesh()
{
    InsertTriangleBasedOnCircumcircle tInsert(domain);

    buildInitialBBXMesh(tInsert);

    initializeTimers();

    for ( int insertedNode = 1; insertedNode <= nnode; insertedNode++ ) {
        PFEMParticle *particle = dynamic_cast< PFEMParticle * >( domain->giveDofManager(insertedNode) );
        if ( particle->isActive() ) {
            std :: list< Edge2D >polygon;

            findNonDelaunayTriangles(insertedNode, tInsert, polygon);

            meshPolygon(insertedNode, tInsert, polygon);
        }
    }

    //giveTimeReport();

    cleanUpTriangleList();

    //deleting BBX-corner particles
    domain->resizeDofManagers(nnode);

    int size = generalTriangleList.size();

    VERBOSE_PRINT0("Number of generated elements", size);
    triangleOctree.giveReport();

    if ( alphaValue > 0.001 ) {
        computeAlphaComplex();
        giveAlphaShape();
    }

    cleanUpTriangleList();

    writeMesh();
}

void DelaunayTriangulator :: writeMesh()
{
    int num = 1;
    int nelem = generalTriangleList.size();
    IntArray dmans(3);
    DofManager *dman;
    DofIDItem type;
    bool hasNoBcOnItself;

    int mat = 0;
    int cs = 0;
    int pressureBC = 0;
    PFEM *pfemEngngModel = dynamic_cast< PFEM * >( domain->giveEngngModel() );
    if ( pfemEngngModel ) {
        mat = pfemEngngModel->giveAssociatedMaterialNumber();
        cs = pfemEngngModel->giveAssociatedCrossSectionNumber();
        pressureBC = pfemEngngModel->giveAssociatedPressureBC();
    }

    domain->resizeElements(nelem);
    //from domain.C
    for ( auto gen : generalTriangleList ) {
        auto elem = std::make_unique<TR1_2D_PFEM>(num, domain, gen->giveNode(1), gen->giveNode(2), gen->giveNode(3), mat, cs);

        domain->setElement(num, std::move(elem));
        num++;
    }

    if ( alphaValue > 0.001 ) {
        // first reset all pressure boundary conditions and alphaShapeProperty
        for ( int i = 1; i <= domain->giveNumberOfDofManagers(); i++ ) {
            Dof *jDof = domain->giveDofManager(i)->giveDofWithID(P_f);
            jDof->setBcId(0);

            dynamic_cast< PFEMParticle * >( domain->giveDofManager(i) )->setOnAlphaShape(false);
        }

        // and then prescribe zero pressure on the free surface
        for ( auto el : alphaShapeEdgeList ) {
            bool oneIsFree = false;
            hasNoBcOnItself = true;
            dman = domain->giveDofManager( el->giveFirstNodeNumber() );
            for ( Dof *dof: *dman ) {
                type = dof->giveDofID();
                if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                    if ( dof->giveBcId() ) {
                        hasNoBcOnItself = false;
                    }
                }
            }

            if ( hasNoBcOnItself ) {
                oneIsFree = true;
            }
            dynamic_cast< PFEMParticle * >( dman )->setOnAlphaShape();

            hasNoBcOnItself = true;
            dman = domain->giveDofManager( el->giveSecondNodeNumber() );
            for ( Dof *dof: *dman ) {
                type = dof->giveDofID();
                if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                    if ( dof->giveBcId() ) {
                        hasNoBcOnItself = false;
                    }
                }
            }

            if ( hasNoBcOnItself ) {
                oneIsFree = true;
            }
            dynamic_cast< PFEMParticle * >( dman )->setOnAlphaShape();

            if ( oneIsFree ) {
                Dof *dofOnNode1 = domain->giveDofManager( el->giveFirstNodeNumber() )->giveDofWithID(P_f);
                dofOnNode1->setBcId(pressureBC);

                Dof *dofOnNode2 = domain->giveDofManager( el->giveSecondNodeNumber() )->giveDofWithID(P_f);
                dofOnNode2->setBcId(pressureBC);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////
void DelaunayTriangulator :: computeAlphaComplex()
{
    bool alphaLengthInitialized = false;
    double minimalLength = 1.e6;

    for ( auto &gen : generalTriangleList ) {
        if ( alphaLengthInitialized ) {
            minimalLength = min( gen->giveShortestEdgeLength(), minimalLength );
        } else {
            minimalLength = gen->giveShortestEdgeLength();
            alphaLengthInitialized = true;
        }

        int par1 = gen->giveNode(1);
        int par2 = gen->giveNode(2);
        int par3 = gen->giveNode(3);
        double ccRadius = gen->giveCircumRadius();

        AlphaEdge2D *containedEdge;

        AlphaEdge2D *edge1 = new AlphaEdge2D( par1, par2, gen->giveEdgeLength(1, 2) );

        containedEdge = giveBackEdgeIfAlreadyContainedInList(*edge1);

        if ( containedEdge ) {
            delete edge1;
            containedEdge->setSharing( 2, gen );
            double outAlph = containedEdge->giveOuterAlphaBound();
            if ( ccRadius < outAlph ) {
                containedEdge->setOuterAlphaBound(ccRadius);
            }

            double innAlph = containedEdge->giveInnerAlphaBound();
            if ( ccRadius > innAlph ) {
                containedEdge->setInnerAlphaBound(ccRadius);
            }

            containedEdge->setHullFlag(false);

            edgeList.push_back(containedEdge);
        } else {
            edge1->setOuterAlphaBound(ccRadius);
            edge1->setInnerAlphaBound(ccRadius);
            edge1->setHullFlag(true);

            edge1->setSharing( 1, gen );
            edgeList.push_back(edge1);
        }

        AlphaEdge2D *edge2 = new AlphaEdge2D( par2, par3, gen->giveEdgeLength(2, 3) );

        containedEdge = giveBackEdgeIfAlreadyContainedInList(*edge2);

        if ( containedEdge ) {
            delete edge2;
            containedEdge->setSharing( 2, gen );

            double outAlph = containedEdge->giveOuterAlphaBound();
            if ( ccRadius < outAlph ) {
                containedEdge->setOuterAlphaBound(ccRadius);
            }

            double innAlph = containedEdge->giveInnerAlphaBound();
            if ( ccRadius > innAlph ) {
                containedEdge->setInnerAlphaBound(ccRadius);
            }

            containedEdge->setHullFlag(false);

            edgeList.push_back(containedEdge);
        } else {
            edge2->setOuterAlphaBound(ccRadius);
            edge2->setInnerAlphaBound(ccRadius);
            edge2->setHullFlag(true);

            edge2->setSharing( 1, gen );
            edgeList.push_back(edge2);
        }

        AlphaEdge2D *edge3 = new AlphaEdge2D( par3, par1, gen->giveEdgeLength(3, 1) );

        containedEdge = giveBackEdgeIfAlreadyContainedInList(*edge3);

        if ( containedEdge ) {
            delete edge3;
            containedEdge->setSharing( 2, gen );

            double outAlph = containedEdge->giveOuterAlphaBound();
            if ( ccRadius < outAlph ) {
                containedEdge->setOuterAlphaBound(ccRadius);
            }

            double innAlph = containedEdge->giveInnerAlphaBound();
            if ( ccRadius > innAlph ) {
                containedEdge->setInnerAlphaBound(ccRadius);
            }

            containedEdge->setHullFlag(false);

            edgeList.push_back(containedEdge);
        } else {
            edge3->setOuterAlphaBound(ccRadius);
            edge3->setInnerAlphaBound(ccRadius);
            edge3->setHullFlag(true);

            edge3->setSharing( 1, gen );
            edgeList.push_back(edge3);
        }
    }

    //alphaValue *= minimalLength;
}

//gives back pointer from edgeList
//////////////////////////////////////////////////////////////////////////
AlphaEdge2D *DelaunayTriangulator :: giveBackEdgeIfAlreadyContainedInList(AlphaEdge2D &alphaEdge)
{
    for ( auto elIT = edgeList.begin(); elIT != edgeList.end(); ++elIT ) {
        if ( ( * * elIT ) == alphaEdge ) {
            AlphaEdge2D *foundEdge = * elIT;
            edgeList.erase(elIT);
            return foundEdge;
        }
    }

    return nullptr;
}

//////////////////////////////////////////////////////////////////////////
void DelaunayTriangulator :: giveAlphaShape()
{
    for ( auto el : edgeList ) {
        // Option 2 : setting bounds vor computed Alpha
        //double alpha = max(min(alphaValue * el->giveLength(), maxAlpha), minAlpha);
        double alpha = alphaValue;
        double outBound = el->giveOuterAlphaBound();

        //innerBound = infinity
        if ( el->giveHullFlag() ) {
            if ( alpha > outBound ) {
                alphaShapeEdgeList.push_back(el);
            } else {
                //invalidating element
                el->giveShared(1)->setValidFlag(false);
            }
        } else {
            double innBound = el->giveInnerAlphaBound();
            if ( alpha > outBound && alpha < innBound ) {
                alphaShapeEdgeList.push_back(el);
                if ( el->giveShared(1)->giveCircumRadius() > alpha ) {
                    el->giveShared(1)->setValidFlag(false);
                } else {
                    el->giveShared(2)->setValidFlag(false);
                }
            }

            if ( alpha < outBound ) {
                for ( int i = 1; i <= 2; i++ ) {
                    el->giveShared(i)->setValidFlag(false);
                }
            }
        }
    }
}

void
DelaunayTriangulator :: initializeTimers()
{
    this->meshingTimer.startTimer();
    this->searchingTimer.startTimer();
    this->searchingTimer.pauseTimer();
    this->polygonTimer.startTimer();
    this->polygonTimer.pauseTimer();
    this->creativeTimer.startTimer();
    this->creativeTimer.pauseTimer();
}

void
DelaunayTriangulator :: findNonDelaunayTriangles(int insertedNode, InsertTriangleBasedOnCircumcircle &tInsert, std :: list< Edge2D > &polygon)
{
    DofManager *node = domain->giveNode(insertedNode);
    const auto &nodeCoords = *node->giveCoordinates();

    ElementCircumCirclesContainingNode findElements(nodeCoords, domain);
    std :: list< DelaunayTriangle * > nonDelaunayTriangles;

    this->searchingTimer.resumeTimer();
    triangleOctree.proceedDataOnFilterAndRemoveFromOctree(nonDelaunayTriangles, findElements, tInsert, searchingTimer);
    this->searchingTimer.pauseTimer();

    this->polygonTimer.resumeTimer();
    for ( auto triangle : nonDelaunayTriangles ) {
        addUniqueEdgeToPolygon({triangle->giveNode(1), triangle->giveNode(2)}, polygon);
        addUniqueEdgeToPolygon({triangle->giveNode(2), triangle->giveNode(3)}, polygon);
        addUniqueEdgeToPolygon({triangle->giveNode(3), triangle->giveNode(1)}, polygon);
        triangle->setValidFlag(false);
    }

    this->polygonTimer.pauseTimer();
}

void
DelaunayTriangulator :: meshPolygon(int insertedNode, InsertTriangleBasedOnCircumcircle &tInsert, std :: list< Edge2D > &polygon)
{
    for ( auto polygonIT = polygon.begin(); polygonIT != polygon.end(); polygonIT++ ) {
        DelaunayTriangle *newTriangle = new DelaunayTriangle(domain, ( * polygonIT ).giveFirstNodeNumber(), ( * polygonIT ).giveSecondNodeNumber(), insertedNode);

        this->creativeTimer.resumeTimer();
        triangleOctree.insertMemberIntoOctree(newTriangle, tInsert);
        this->creativeTimer.pauseTimer();

        generalTriangleList.push_back(newTriangle);
    }

    polygon.clear();
}

void
DelaunayTriangulator :: giveTimeReport()
{
    this->meshingTimer.stopTimer();
    double _utime = this->meshingTimer.getUtime();

    this->searchingTimer.stopTimer();
    double _searchutime = this->searchingTimer.getUtime();

    this->polygonTimer.stopTimer();
    double _polygonutime = this->polygonTimer.getUtime();

    this->creativeTimer.stopTimer();
    double _creativeutime = this->creativeTimer.getUtime();

    printf("\nUser time consumed by searching: %.3f [s]\n\n", _searchutime);
    printf("\nUser time consumed by building insertion polygon: %.3f [s]\n\n", _polygonutime);
    printf("\nUser time consumed by creating elements: %.3f [s]\n\n", _creativeutime);
    printf("\nUser time consumed by meshing: %.3f [s]\n\n", _utime);
}

void
DelaunayTriangulator :: cleanUpTriangleList()
{
    for ( auto genIT = generalTriangleList.begin(); genIT != generalTriangleList.end(); ) {
        int node1 = ( * genIT )->giveNode(1);
        int node2 = ( * genIT )->giveNode(2);
        int node3 = ( * genIT )->giveNode(3);

        if ( node1 == nnode + 1 || node1 == nnode + 2 || node1 == nnode + 3 || node1 == nnode + 4 ||
             node2 == nnode + 1 || node2 == nnode + 2 || node2 == nnode + 3 || node2 == nnode + 4 ||
             node3 == nnode + 1 || node3 == nnode + 2 || node3 == nnode + 3 || node3 == nnode + 4 ||
             !( ( * genIT )->giveValidFlag() ) ) {
            delete( * genIT );
            genIT = generalTriangleList.erase(genIT);
        } else {
            genIT++;
        }
    }
}

void DelaunayTriangulator :: buildInitialBBXMesh(InsertTriangleBasedOnCircumcircle &tInsert)
{
    BoundingBox BBX;
    this->computeBBXBasedOnNodeData(BBX);

    //int initialDivision = 5;
    triangleOctree.init(BBX);

    FloatArray coords;

    auto bottomLeftNode = std::make_unique<Node>(nnode + 1, domain);
    BBX.giveOrigin(coords);
    double diff = BBX.giveSize();
    for ( int ci = 1; ci <= 2; ci++ ) {
        coords.at(ci) -= diff;
    }

    bottomLeftNode->setCoordinates(coords);

    auto bottomRightNode = std::make_unique<Node>(nnode + 2, domain);
    coords.at(1) += BBX.giveSize() + 2.0 * diff;
    bottomRightNode->setCoordinates(coords);

    auto topRightNode = std::make_unique<Node>(nnode + 3, domain);
    coords.at(2) += BBX.giveSize() + 2.0 * diff;
    topRightNode->setCoordinates(coords);

    auto topLeftNode = std::make_unique<Node>(nnode + 4, domain);
    coords.at(1) -= BBX.giveSize() + 2.0 * diff;
    topLeftNode->setCoordinates(coords);

    domain->resizeDofManagers(nnode + 4);
    domain->setDofManager(nnode + 1, std::move(bottomLeftNode));
    domain->setDofManager(nnode + 2, std::move(bottomRightNode));
    domain->setDofManager(nnode + 3, std::move(topRightNode));
    domain->setDofManager(nnode + 4, std::move(topLeftNode));

    DelaunayTriangle *firstTriangle = new DelaunayTriangle(domain, nnode + 1, nnode + 2, nnode + 3);
    generalTriangleList.push_back(firstTriangle);

    DelaunayTriangle *secondTriangle = new DelaunayTriangle(domain, nnode + 3, nnode + 4, nnode + 1);
    generalTriangleList.push_back(secondTriangle);

    triangleOctree.insertMemberIntoOctree(firstTriangle, tInsert);
    triangleOctree.insertMemberIntoOctree(secondTriangle, tInsert);
}

void
DelaunayTriangulator:: computeBBXBasedOnNodeData(BoundingBox &BBX)
{
    int init = 1, nnode = this->domain->giveNumberOfDofManagers();
    FloatArray minc(3), maxc(3);

    // first determine domain extends (bounding box), and check for degenerated domain type
    for ( int i = 1; i <= nnode; i++ ) {
        auto node = domain->giveNode(i);
        const auto &coords = node->giveCoordinates();
        if ( init ) {
            init = 0;
            for ( int j = 1; j <= coords->giveSize(); j++ ) {
                minc.at(j) = maxc.at(j) = coords->at(j);
            }
        } else {
            for ( int j = 1; j <= coords->giveSize(); j++ ) {
                if ( coords->at(j) < minc.at(j) ) {
                    minc.at(j) = coords->at(j);
                }
                if ( coords->at(j) > maxc.at(j) ) {
                    maxc.at(j) = coords->at(j);
                }
            }
        }
    }                 // end loop over nodes

    BBX.setOrigin(minc);

    // determine root size
    double rootSize = 0.0;
    for ( int i = 1; i <= 3; i++ ) {
        rootSize = 1.000001 * max( rootSize, maxc.at(i) - minc.at(i) );
    }

    BBX.setSize(rootSize);

    // check for degenerated domain
    double resolutionLimit = min(1.e-3, rootSize / 1.e6);
    for ( int i = 1; i <= 3; i++ ) {
        if ( ( maxc.at(i) - minc.at(i) ) > resolutionLimit ) {
            BBX.setMask(i, 1);
        } else {
            BBX.setMask(i, 0);
        }
    }
}


} // end namespace oofem
