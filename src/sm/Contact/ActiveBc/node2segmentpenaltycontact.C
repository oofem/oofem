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


#include "sm/Contact/ActiveBc/node2segmentpenaltycontact.h"
#include "set.h"
#include "domain.h"
#include "node.h"
#include "floatmatrix.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"
#include "classfactory.h"
#include "mathfem.h"

namespace oofem {
REGISTER_BoundaryCondition( Node2SegmentPenaltyContact );



void Node2SegmentPenaltyContact::computeTangentFromContact( FloatMatrix &answer, Node *node, ContactSegment *segment, TimeStep *tStep )
{
    double gap;
    this->computeGap( gap, node, segment, tStep );

    //assembling the first part of the tangent
    //considering Nv = N^T * n/||n||
    //and         K1 = p * (Nv*Nv^T)
    // (this part is always assembled as it is the simplest way to know the required dimensions of answer)
    // (if there is no contact, the answer, already with the right dimensions, is later zeroed, see below)
    FloatArray Nv;
    this->computeNvMatrixAt( Nv, node, segment, tStep );
    answer.beDyadicProductOf( Nv, Nv );
    answer.times( this->penalty );

    if ( gap < 0.0 && segment->hasNonLinearGeometry(node, tStep)) {
        //assembling the rest of the tangent (for large deformations)

        //K2 = - p g/l (BvTv + TvBv + g/l BvBv)
        FloatMatrix k2, k3, k4;
        FloatArray Bv, Tv;
        double l;

        this->computeBvMatrixAt( Bv, node, segment, tStep );
        this->computeTvMatrixAt( Tv, node, segment, tStep );

		FloatMatrix m;
        segment->computeMetricTensor( m, node, tStep );
        //placeholder for 2D only
        l = sqrt( m.at( 1, 1 ) );


		k2.beDyadicProductOf( Bv, Tv );
        k2.times( -penalty * gap  / l );
	
		k3.beDyadicProductOf( Tv, Bv );
        k3.times( -penalty * gap / l );
	
		k4.beDyadicProductOf( Bv, Bv );
        k4.times( -penalty * gap * gap / (l*l) );

		answer.add( k2 );
        answer.add( k3 );
        answer.add( k4 );
    }
    if (gap >= 0.0){
		//zero in the case of no contact occuring
        answer.zero();
    }
}


void Node2SegmentPenaltyContact::computeExternalForcesFromContact( FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *tStep )
{
    double gap;
    this->computeGap( gap, node, segment, tStep );
    this->computeNvMatrixAt( answer, node, segment, tStep );
    if ( gap < 0.0 ) {
        answer.times( -penalty * gap );
    } else {
        answer.times( 0 );
    }
}


void Node2SegmentPenaltyContact::computeGap( double &answer, Node *node, ContactSegment *segment, TimeStep *tStep )
{
    answer = segment->computePenetration( node, tStep );
}


void Node2SegmentPenaltyContact::computeNvMatrixAt( FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *tStep )
{

    FloatArray normal;
    FloatMatrix segmentN, extendedN;
	int ndof = node->giveNumberOfDofs();

    if ( prescribedNormal.giveSize() == node->giveNumberOfDofs() ) {
        normal = prescribedNormal;
    } else {
        segment->computeNormal( normal, node, tStep );
    }
    int norm = normal.computeNorm();
    if(norm) {
      normal.normalize();
    }

    segment->computeSegmentNMatrix( segmentN, node, tStep );

	if (segmentN.giveNumberOfRows() != ndof) {
		OOFEM_ERROR("Dimension mismatch between node and contact segment");
	}

	//append extesion
	extendedN.resize(segmentN.giveNumberOfRows(), segmentN.giveNumberOfColumns() + 2);
	extendedN.zero();
	FloatMatrix extension(ndof, ndof);
	extension.beUnitMatrix();
	segmentN.times(-1);
	extendedN.setSubMatrix(extension, 1, 1);
	extendedN.setSubMatrix(segmentN, 1, ndof + 1);

    //Nv should be given just as N^t * n;
    answer.beTProductOf(extendedN, normal );
}

void Node2SegmentPenaltyContact::computeTvMatrixAt( FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *tStep )
{
    FloatArray tangent;
	FloatMatrix segmentN, extendedN;
	int ndof = node->giveNumberOfDofs();

    if ( prescribedNormal.giveSize() == node->giveNumberOfDofs() ) {
        OOFEM_WARNING( "Prescribed normal inapplicable for use with large strains" );
    }
    segment->computeTangent( tangent, node, tStep );
    tangent.normalize();

    segment->computeSegmentNMatrix(segmentN, node, tStep );

	if (segmentN.giveNumberOfRows() != ndof) {
		OOFEM_ERROR("Dimension mismatch between node and contact segment");
	}

	extendedN.resize(segmentN.giveNumberOfRows(), segmentN.giveNumberOfColumns() + 2);
	extendedN.zero();
	FloatMatrix extension(ndof, ndof);
	extension.beUnitMatrix();
	segmentN.times(-1);
	extendedN.setSubMatrix(extension, 1, 1);
	extendedN.setSubMatrix(segmentN, 1, ndof+1);

    //Tv should be given just as N^t * t;
    answer.beTProductOf( extendedN, tangent );
}

void Node2SegmentPenaltyContact::computeBvMatrixAt( FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *tStep )
{
    FloatArray normal;
    FloatMatrix segmentB, extendedB;
	int ndof = node->giveNumberOfDofs();

    if ( prescribedNormal.giveSize() == node->giveNumberOfDofs() ) {
        normal = prescribedNormal;
    } else {
        segment->computeNormal( normal, node, tStep );
    }
    normal.normalize();

    segment->computeSegmentBMatrix(segmentB, node, tStep );

	if (segmentB.giveNumberOfRows() != ndof) {
		OOFEM_ERROR("Dimension mismatch between node and contact segment");
	}

	extendedB.resize(segmentB.giveNumberOfRows(), segmentB.giveNumberOfColumns() + 2);
	extendedB.zero();
	FloatMatrix extension(ndof, ndof);
	segmentB.times(-1);
	extendedB.setSubMatrix(extension, 1, 1);
	extendedB.setSubMatrix(segmentB, 1, ndof + 1);

    //Bv should be given just as B^t * n;
    answer.beTProductOf( extendedB, normal );
}




IRResultType Node2SegmentPenaltyContact::initializeFrom( InputRecord *ir )
{
    IRResultType result;
    IR_GIVE_FIELD( ir, this->penalty, _IFT_Node2SegmentPenaltyContact_penalty );
    this->useTangent = ir->hasField( _IFT_Node2SegmentPenaltyContact_useTangent );


    IR_GIVE_FIELD( ir, this->segmentSet, _IFT_Node2SegmentPenaltyContact_segmentSet );
    IR_GIVE_FIELD( ir, this->nodeSet, _IFT_Node2SegmentPenaltyContact_nodeSet );

    IR_GIVE_OPTIONAL_FIELD( ir, this->prescribedNormal, _IFT_Node2SegmentPenaltyContact_prescribedNormal );

    return ActiveBoundaryCondition::initializeFrom( ir );
}


void Node2SegmentPenaltyContact::assemble( SparseMtrx &answer, TimeStep *tStep,
    CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale )
{
    if ( !this->useTangent || type != TangentStiffnessMatrix ) {
        return;
    }

    FloatMatrix K;
    IntArray loc, node_loc;

    IntArray dofIdArray = giveDomain()->giveDefaultNodeDofIDArry();

    //iterate over all pairs of nodes and segments
    for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); ++nodePos ) {
        for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {

            Node *node              = this->giveDomain()->giveNode( nodeSet.at( nodePos ) );
            ContactSegment *segment = ( this->giveDomain()->giveContactSegment( segmentSet.at( segmentPos ) ) );
            this->computeTangentFromContact( K, node, segment, tStep );
	    this->giveLocationArray(loc, r_s, node, segment);
            answer.assemble( loc, K );
        }
    }
}


void Node2SegmentPenaltyContact::assembleVector( FloatArray &answer, TimeStep *tStep,
    CharType type, ValueModeType mode,
    const UnknownNumberingScheme &s, FloatArray *eNorms )
{
    if ( type != ExternalForcesVector ) {
        return;
    }


    IntArray loc;
    FloatArray fext;

    for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); ++nodePos ) {
        for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {
            Node *node              = this->giveDomain()->giveNode( nodeSet.at( nodePos ) );
            ContactSegment *segment = ( this->giveDomain()->giveContactSegment( segmentSet.at( segmentPos ) ) );
            this->computeExternalForcesFromContact( fext, node, segment, tStep );
	    this->giveLocationArray(loc, s, node, segment);
            answer.assemble( fext, loc );
        }
    }
}

void Node2SegmentPenaltyContact:: giveLocationArray(IntArray &loc, const UnknownNumberingScheme &ns, const Node *node, const ContactSegment *segment)
{
  IntArray seg_loc;
  IntArray dofIdArray = giveDomain()->giveDefaultNodeDofIDArry();
  //assembling for both node and segment
  node->giveLocationArray( dofIdArray, loc, ns );
  segment->giveLocationArray( dofIdArray, seg_loc, ns);
  loc.followedBy( seg_loc );
}

void Node2SegmentPenaltyContact::giveLocationArrays( std::vector<IntArray> &rows, std::vector<IntArray> &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s )
{
    //returns all possible combinations of dof that can theoretically be triggered by contact
    //of any segment with any node. Room for optimization aplenty...
    IntArray n_loc, s_loc;

    int ncombinations = nodeSet.giveSize() * segmentSet.giveSize();
    rows.resize( ncombinations * 2 );
    cols.resize( ncombinations * 2 );

    IntArray dofIdArray = giveDomain()->giveDefaultNodeDofIDArry();

    int pos = 0;

    for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); nodePos++ ) {
        for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {
            Node *node              = this->giveDomain()->giveNode( nodeSet.at( nodePos ) );
            ContactSegment *segment = ( this->giveDomain()->giveContactSegment( segmentSet.at( segmentPos ) ) );

            node->giveLocationArray( dofIdArray, n_loc, r_s );
            segment->giveLocationArrays( dofIdArray, s_loc, c_s );

            // insert location arrays into the answer arrays
            rows[pos]     = n_loc;
            cols[pos]     = s_loc;
            rows[pos + 1] = s_loc;
            cols[pos + 1] = n_loc;
            pos += 2;
        }
    }
}

  

} // namespace oofem
