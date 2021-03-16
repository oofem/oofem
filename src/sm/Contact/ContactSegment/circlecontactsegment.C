#include "circlecontactsegment.h"

namespace oofem {
    REGISTER_ContactSegment(CircleContactSegment);

    IRResultType CircleContactSegment::initializeFrom(InputRecord* ir) {
        //read circle parameters
        IRResultType result;
        IR_GIVE_FIELD(ir, centerPoint, _IFT_CircleContactSegment_centerpoint);
        IR_GIVE_FIELD(ir, radius, _IFT_CircleContactSegment_radius);

        return FunctionContactSegment::initializeFrom(ir);
    }

    void CircleContactSegment::computeContactPoint(FloatArray & answer, FloatArray& normal, const FloatArray & nodeCoords)
    {
        if ( nodeCoords.giveSize() != centerPoint.giveSize() ) {
            OOFEM_ERROR("Node coordinate dimension (%i) does not match circle/sphere center point dimension (%i)", nodeCoords.giveSize(), centerPoint.giveSize());
        }
        /*normal.beDifferenceOf(centerPoint, nodeCoords);
        double centerDistance = answer.computeNorm();
        normal.times((centerDistance - radius) / centerDistance);*/


        //foolproof (hopefully) approach
        //first compute projection from node to center of circle
        FloatArray projection;
        projection.beDifferenceOf(centerPoint, nodeCoords);

        //normalized projection is the normal vector
        normal = projection;
        normal.normalize();

        normal.times(-1.);

        //now compute contact point by adding to node coords the distance to circle
        double centerDistance = projection.computeNorm();
        projection.times((centerDistance - radius) / centerDistance);
        answer = nodeCoords;
        answer.add(projection);
    }
}//end namespace oofem