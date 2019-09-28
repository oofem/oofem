#include "boundarycontactsegment.h"

namespace oofem {
	IRResultType oofem::BoundaryContactSegment::initializeFrom(InputRecord * ir)
	{
		IRResultType result;
		IR_GIVE_FIELD(ir, setnum, _IFT_BoundaryContactSegment_edgeSet);

		updateMode = UM_EachStep;
		int updateModeInt = 0;
		IR_GIVE_OPTIONAL_FIELD(ir, updateModeInt, _IFT_BoundaryContactSegment_pairUpdateMode);
		if (result == IRRT_OK) {
			updateMode = (UpdateMode)updateModeInt;
			if (updateModeInt < 0 || updateModeInt > 2) OOFEM_ERROR("Contact segment pairing update mode can be only 0, 1 or 2");
		}
		return ContactSegment::initializeFrom(ir);
	}

	bool BoundaryContactSegment::hasNonLinearGeometry(Node *node, TimeStep *tStep)
	{
		IntArray closestEdge;
		giveClosestEdge(closestEdge, node, tStep);
		if (closestEdge.giveSize() != 2) {
			//no closest edge means no contact
			return false;
		}
		NLStructuralElement *elem = dynamic_cast<NLStructuralElement *> (this->giveDomain()->giveElement(closestEdge.at(1)));
		return elem != nullptr && elem->giveGeometryMode() == 1;
	}

	void BoundaryContactSegment::postInitialize()
	{
		Set* set = this->giveDomain()->giveSet(this->setnum);
		if (set == nullptr) OOFEM_ERROR("Contact segment can not find set no. " + setnum);
		this->edges = set->giveBoundaryList();
		if (edges.giveSize() <= 0) OOFEM_WARNING("Contact segment's edge list is empty");
	}

	void BoundaryContactSegment::updateYourself(TimeStep *tStep)
		// Updates the receiver at end of step.
	{
		if (updateMode != UM_Never) {
			knownNodes.clear();
			knownClosestEdges.clear();
		}
		lastEdge.resize(0);
	}

	void BoundaryContactSegment::giveClosestEdge(IntArray & answer, Node * node, TimeStep * tStep)
	{
		int knownIndex = giveIndexOfKnownNode(node);
		if (updateMode != UM_EveryQuery && knownIndex != -1 && knownClosestEdges.at(knownIndex).giveSize() == 2) {
			answer = knownClosestEdges.at(knownIndex);
			lastEdge = answer;
			return;
		}
		//if previous if failed, it means that we don't know the closest edge to this node yet
		answer.resize(2);

		FloatArray cPoint, cPointLocal, projection, nodeCoords;
		double answerSize = -1.;

		node->giveUpdatedCoordinates(nodeCoords, tStep);

		//iterate over all edges to find the closest one
		for (int pos = 0; pos < edges.giveSize() / 2; pos++) {
			StructuralElement* element = (StructuralElement*)this->giveDomain()->giveElement(edges(pos * 2));
			int edgePos = pos * 2 + 1;

			bool inbetween = computeContactPoint(cPointLocal, node, element, edges(edgePos), tStep);
			FEInterpolation2d* interpolation = dynamic_cast<FEInterpolation2d*>(element->giveInterpolation());
			interpolation->boundaryEdgeLocal2Global(cPoint, edges(edgePos), cPointLocal, FEIElementDeformedGeometryWrapper(element, tStep));

			projection.beDifferenceOf(cPoint, nodeCoords);

			if (inbetween) {
				double normalSize = projection.computeNorm();
				//@todo: solve this issue
				if (answerSize == -1. || normalSize < answerSize) {
					//if ( answerSize == -1. || fabs(normalSize) < answerSize ) {
					//new minimum found, update answer
					answerSize = normalSize;
					answer(0) = edges(pos * 2);
					answer(1) = edges(edgePos);
				}
			}
		}
		//if answer size is still -1 it means no edge is properly in contact (no inbetween projections)
		if (answerSize == -1.) answer.resize(0);
		else {
			knownNodes.push_back(node);
			IntArray answerToStore(answer); //make copy to store to prevent it going out of scope - necessary??
			knownClosestEdges.push_back(answerToStore);
		}
		lastEdge = answer;
	}

}//end namespace oofem
