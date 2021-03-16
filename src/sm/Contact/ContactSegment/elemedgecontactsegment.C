#include "elemedgecontactsegment.h"
#include "classfactory.h"

namespace oofem {

    REGISTER_ContactSegment(Linear2dElementEdgeContactSegment);


	IRResultType Linear2dElementEdgeContactSegment::initializeFrom(InputRecord * ir)
	{
		IRResultType result;
		////IR_GIVE_FIELD(ir, this->elemSet, _IFT_ElementEdgeContactSegment_elemSet);
		//IR_GIVE_FIELD(ir, setnum, _IFT_Linear2dElementEdgeContactSegment_edgeSet);
		//
		//updateMode = UM_EachStep;
        //int updateModeInt = 0;
		//    IR_GIVE_OPTIONAL_FIELD(ir, updateModeInt, _IFT_Linear2dElementEdgeContactSegment_pairUpdateMode);
		//    if ( result == IRRT_OK ) {
		//        updateMode = (UpdateMode)updateModeInt;
		//        if ( updateModeInt < 0 || updateModeInt > 2 ) OOFEM_ERROR("Contact segment pairing update mode can be only 0, 1 or 2");
		//   }

        return BoundaryContactSegment::initializeFrom(ir);
    }

    void Linear2dElementEdgeContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep* tStep)
    {
		if (hasNonLinearGeometry(node, tStep)) {
			//determine normal from tangent as t vector cross product e3
			FloatArray tangent, tangent3D, e3, normal3D;
			computeTangent(tangent, node, tStep);
			tangent3D.resize(3);
			tangent3D.addSubVector(tangent, 1);
			e3.resize(3);
			e3.at(3) = 1.;
			normal3D.beVectorProductOf(e3, tangent3D);
			answer = normal3D;
			answer.times(-1./normal3D.computeNorm());
			answer.resizeWithValues(2);
		}
		else {
			//determine normal from edge
			IntArray closestEdge;
			giveClosestEdge(closestEdge, node, tStep);
			if (closestEdge.giveSize() != 2) {
				//no closest edge means no contact
				//return zeros
				answer.resize(2);
				return;
			}

			StructuralElement* elem = (StructuralElement*)this->giveDomain()->giveElement(closestEdge.at(1));
			int edgePos = closestEdge.at(2);

			FloatArray cPoint, cPointLocal, normal;

			bool inbetween = computeContactPoint(cPointLocal,node, elem, edgePos, tStep);
			//no need to care here whether distance is negative or not

			//retrieve edge normal from element interpolation
			FEInterpolation2d* interpolation = dynamic_cast<FEInterpolation2d*>(elem->giveInterpolation());
			if (interpolation == nullptr) {
				OOFEM_ERROR("Non-2D element encountered in Linear2dElementEdgeContactSegment");
			}
			interpolation->edgeEvalNormal(normal, edgePos, cPointLocal, FEIElementGeometryWrapper(elem));

			answer = normal;
		}
    }

    void Linear2dElementEdgeContactSegment::computeTangent( FloatArray &answer, Node *node, TimeStep *tStep )
    {
        IntArray closestEdge;
        giveClosestEdge( closestEdge, node, tStep );
        if ( closestEdge.giveSize() != 2 ) {
            //no closest edge means no contact
            //return zeros
            answer.resize( 2 );
            return;
        }

		StructuralElement *elem = (StructuralElement *)this->giveDomain()->giveElement( closestEdge.at( 1 ) );
        int edgePos = closestEdge.at( 2 );
        IntArray edgeNodes;
        elem->giveBoundaryEdgeNodes( edgeNodes, edgePos );
        FloatArray nodalCoords;
        FloatArray nodalCoords2;

		elem->giveNode( edgeNodes( 0 ) )->giveUpdatedCoordinates( nodalCoords, tStep );
        elem->giveNode( edgeNodes( 1 ) )->giveUpdatedCoordinates( nodalCoords2, tStep );

		nodalCoords.append(nodalCoords2);
        FloatMatrix dNdXi;
        computedNdksi( dNdXi, node, tStep );
		answer.beProductOf( dNdXi, nodalCoords );
    }

    void Linear2dElementEdgeContactSegment::computeSegmentNMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {
        IntArray closestEdge;
        giveClosestEdge(closestEdge, node, tStep);
        if ( closestEdge.giveSize() != 2 ) {
            //no closest edge means no contact
            //return zeros
            answer.resize(2, 4);
            return;
        }

        StructuralElement* elem = (StructuralElement*)this->giveDomain()->giveElement(closestEdge.at(1));
        int edgePos = closestEdge.at(2);

        FloatMatrix N;
        FloatArray cPointLocal, nodeCoords, edgeNode1Coords, edgeNode2Coords;
        IntArray edgeNodes;
        elem->giveBoundaryEdgeNodes(edgeNodes, edgePos);

        node->giveUpdatedCoordinates(nodeCoords, tStep);
        elem->giveNode(edgeNodes(0))->giveUpdatedCoordinates(edgeNode1Coords, tStep);
		elem->giveNode(edgeNodes(1))->giveUpdatedCoordinates(edgeNode2Coords, tStep);
        bool inbetween = computeContactPoint(cPointLocal, node, elem, edgePos, tStep);
        //all the previous just to compute the contact point...
        elem->computeEdgeNMatrix(N, edgePos, cPointLocal);

        /*answer.resize(N.giveNumberOfRows(), N.giveNumberOfColumns() + 2);
        answer.zero();

        FloatMatrix extension(2, 2);
        extension.beUnitMatrix();

		N.times(-1);
        answer.setSubMatrix(extension, 1, 1);
        answer.setSubMatrix(N, 1, 3);*/

		//change to have the extensions added later in the contact condition class
		answer = N; 
	
    }

    void Linear2dElementEdgeContactSegment::computeSegmentBMatrix( FloatMatrix &answer, Node *node, TimeStep *tStep )
    {
		//for linear segments, this is always the same
		/*answer = {{0,0},{0,0}, {-0.5, 0},{0, -0.5},{0.5,0},{0,0.5}};
		answer.times(2);*/

		//change to have the extensions added later in the contact condition classs
		computedNdksi(answer, node, tStep);
	}

    void Linear2dElementEdgeContactSegment::computedNdksi( FloatMatrix &answer, Node *node, TimeStep *tStep )
    {
		//for linear segments, this is always the same
		answer = {{-0.5, 0},{0, -0.5},{0.5,0},{0,0.5}};
        answer.times( 2 );
    }
  

    double Linear2dElementEdgeContactSegment::computePenetration(Node * node, TimeStep * tStep)
    {
        IntArray closestEdge;
        giveClosestEdge(closestEdge, node, tStep);
        if ( closestEdge.giveSize() != 2 ) {
            //no closest edge means no contact
            return 0.0;
        }

		FloatArray cPointLocal, projection, nodeCoords, edgeNode1Coords;
		IntArray edgeNodes;

        node->giveUpdatedCoordinates(nodeCoords, tStep);


		StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(closestEdge.at(1));

		element->giveBoundaryEdgeNodes(edgeNodes, closestEdge.at(2));
		element->giveNode(edgeNodes(0))->giveUpdatedCoordinates(edgeNode1Coords, tStep);

        bool inbetween = computeContactPoint(cPointLocal, node, element, closestEdge.at(2), tStep);
		if(inbetween == false) {
			 return 0;
		}
        //projection.beDifferenceOf(nodeCoords, cPoint);
		projection.beDifferenceOf(nodeCoords, edgeNode1Coords);

        ////test whether initial and current vector are on different sides of line
		FloatArray normal;
		this->computeNormal(normal, node, tStep);
		return projection.dotProduct(normal);
    }

    void Linear2dElementEdgeContactSegment::computeMetricTensor( FloatMatrix &answer, Node *node, TimeStep *tStep )
    {
        answer.resize( 1, 1 );
        //the answer is length of segment, which is in fact the size of tangent, squared
		FloatArray tangent;
        computeTangent( tangent, node, tStep );
		double l = tangent.computeNorm();
        answer.at( 1, 1 ) = l * l;
    }

    void Linear2dElementEdgeContactSegment::giveLocationArray(const IntArray & dofIdArray, IntArray & answer, const UnknownNumberingScheme & c_s) const
    {
        if ( lastEdge.giveSize() == 2 ) {
            StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(lastEdge.at(1));
            int edgePos = lastEdge.at(2);

            IntArray boundaryNodes;
            element->giveBoundaryEdgeNodes(boundaryNodes, edgePos);
            element->giveBoundaryLocationArray(answer, boundaryNodes, c_s, nullptr);
        } else {
            //if no segment was worked with, returns zeros
            answer.resize(4);
        }

    }

    void Linear2dElementEdgeContactSegment::giveLocationArrays(const IntArray & dofIdArray, IntArray & answer, const UnknownNumberingScheme & c_s)
    {
        answer.resize(0);
        IntArray edgeloc, boundaryNodes;
        //iterate over all edges and add their locarrays
        for ( int pos = 0; pos < edges.giveSize() / 2; pos++ ) {
            StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(edges(pos * 2));
            int edgePos = pos * 2 + 1;

            element->giveBoundaryEdgeNodes(boundaryNodes, edges(edgePos));
            element->giveBoundaryLocationArray(edgeloc, boundaryNodes, c_s, nullptr);

            answer.followedBy(edgeloc);
        }
    }

	bool Linear2dElementEdgeContactSegment::computeContactPoint(FloatArray &ksi, Node* node, StructuralElement* element, int elemedge, TimeStep* tStep){
		//@todo: general function for the closest point projection should be introduced in a parent class
		// Really? This all seems pretty specific for 2D linear segments. Moved the definition, this is now an inherited protected virtual function

		//retrieve coordinates of points in deformed configuration
		FloatArray externalPoint, linePoint1, linePoint2;
		node->giveUpdatedCoordinates( externalPoint, tStep );

		IntArray edgeNodes;
		element->giveBoundaryEdgeNodes( edgeNodes, elemedge );

		element->giveNode(edgeNodes(0))->giveUpdatedCoordinates(linePoint1, tStep);
		element->giveNode(edgeNodes(1))->giveUpdatedCoordinates(linePoint2, tStep);

		//compute the contact point
		ksi.resize(1);
		FloatArray tangent;
		tangent.beDifferenceOf(linePoint2, linePoint1);
        double l = tangent.computeNorm();
		tangent.times(1./l);
		//compute vector x_s - x_1
		FloatArray v_s1;
		v_s1.beDifferenceOf(externalPoint, linePoint1);
		// compute parametric coord of the contact point by projection of x_s-x_1 onto tangent and divide by l (in interval -1 1)
		ksi.at(1) = 2. * tangent.dotProduct(v_s1) / l - 1.;
		return (ksi.at(1) <= 1 && ksi.at(1) >= -1);
    }



} //end namespace oofem


