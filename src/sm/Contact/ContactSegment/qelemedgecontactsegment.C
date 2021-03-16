
#include "qelemedgecontactsegment.h"

namespace oofem {
    REGISTER_ContactSegment(QElementEdgeContactSegment);

    IRResultType QElementEdgeContactSegment::initializeFrom(InputRecord* ir) {
        IRResultType result;
        //IR_GIVE_FIELD(ir, this->elemSet, _IFT_ElementEdgeContactSegment_elemSet);
        IR_GIVE_FIELD(ir, setnum, _IFT_QElementEdgeContactSegment_edgeSet);
        return ContactSegment::initializeFrom(ir);
    }

    void QElementEdgeContactSegment::computeNormal(FloatArray & answer, Node * node, TimeStep * tstep)
    {

    }

    void QElementEdgeContactSegment::computeTangent( FloatArray &answer, Node *node, TimeStep *tstep )
    {
    }

    void QElementEdgeContactSegment::computeSegmentNMatrix(FloatMatrix & answer, Node * node, TimeStep * tStep)
    {

    }

    double QElementEdgeContactSegment::computePenetration(Node * node, TimeStep * tStep)
    {
        return 0.0;
    }

    void QElementEdgeContactSegment::giveLocationArray(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s) const 
    {

    }

    void QElementEdgeContactSegment::giveLocationArrays(const IntArray & dofIdArray, IntArray & s_loc, const UnknownNumberingScheme & c_s)
    {
    }

    void QElementEdgeContactSegment::updateYourself(TimeStep * tStep)
    {
        knownNodes.clear();
        knownClosestEdges.clear();
        lastEdge.resize(0);
    }

    void QElementEdgeContactSegment::postInitialize()
    {
        Set* set = this->giveDomain()->giveSet(this->setnum);
        if ( set == nullptr ) OOFEM_ERROR("Contact segment can not find set no. " + setnum);
        this->edges = set->giveBoundaryList();
        if ( edges.giveSize() <= 0 ) OOFEM_WARNING("Contact segment's edge list is empty");
    }

    void QElementEdgeContactSegment::giveClosestEdge(IntArray & answer, Node * node, TimeStep * tStep)
    {
    }

}
