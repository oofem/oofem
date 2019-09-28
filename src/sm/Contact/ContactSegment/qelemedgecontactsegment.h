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


#pragma once
#include "contactsegment.h"
#include "node.h"
#include "inputrecord.h"
#include "Elements/structuralelement.h"
#include "set.h"
#include "feinterpol.h"
#include "classfactory.h"

#define _IFT_QElementEdgeContactSegment_Name "qelementedgecontactsegment"
#define _IFT_QElementEdgeContactSegment_edgeSet "edgeset"

namespace oofem {
class QElementEdgeContactSegment : public ContactSegment
{
public:
    QElementEdgeContactSegment( int n, Domain *aDomain ) :
        ContactSegment( n, aDomain ) { ; }
    ~QElementEdgeContactSegment(){};

    IRResultType initializeFrom( InputRecord *ir );

    //returns normalized n, which is an normal vector of contact
    void computeNormal( FloatArray &answer, Node *node, TimeStep *tstep ) override;

    //returns normalized n, which is an normal vector of contact
    void computeTangent( FloatArray &answer, Node *node, TimeStep *tstep ) override;

    //returns an extended N (aka A) matrix, integrated at point of contact of given node
    void computeSegmentNMatrix( FloatMatrix &answer, Node *node, TimeStep *tStep ) override;

    //computes the penetration of node given
    double computePenetration( Node *node, TimeStep *tStep ) override;

    //returns a xi derivative of the extended N matrix at point of contact of given node
    virtual void computeSegmentBMatrix( FloatMatrix &answer, Node *node, TimeStep *tStep ) override{};

	//returns the metric tensor m
    virtual void computeMetricTensor( FloatMatrix &answer, Node *node, TimeStep *tStep ) override{};

    //determines whether geometrical non-linearity should be considered
    virtual bool hasNonLinearGeometry( Node *node, TimeStep *tStep ) override { return false; };

    void giveLocationArray( const IntArray &dofIdArray, IntArray &s_loc, const UnknownNumberingScheme &c_s ) const override;
    void giveLocationArrays( const IntArray &dofIdArray, IntArray &s_loc, const UnknownNumberingScheme &c_s ) override;

    void updateYourself( TimeStep *tStep ) override;

    const char *giveClassName() const override { return "Qelemenetedgecontactsegment"; }
    const char *giveInputRecordName() const override { return _IFT_QElementEdgeContactSegment_Name; }

    void postInitialize() override;

private:
    IntArray edges;
    IntArray lastEdge; //last edge worked with
    std::vector<Node *> knownNodes;
    std::vector<IntArray> knownClosestEdges;
    int setnum;

    //gives the closest edge to a given node in the form of an IntArray(elempos,edgepos)
    //only computes it again if it wasn't determined for this node in this solution step yet
    //(i.e. the array of known closest edges is reset after step convergence)
    void giveClosestEdge( IntArray &answer, Node *node, TimeStep *tStep );

    //computes distance of a point given by a set of coordinates to a line given by two sets of coordinates
    //returns whether point of intersection is inbetween the line points
    //bool computeContactPoint(FloatArray& answer, const FloatArray& externalPoint, const FloatArray& linePoint1, const FloatArray& linePoint2, FloatArray * contactPointCoords = nullptr);

    //searches the array of nodes to which a closest edge was already computed
    //returns -1 if unsuccessful
    inline int giveIndexOfKnownNode( const Node *node )
    {
        for ( int i = 0; i < knownNodes.size(); i++ ) {
            if ( node == knownNodes.at( i ) ) return i;
        }
        return -1;
    }
};

} // namespace oofem
