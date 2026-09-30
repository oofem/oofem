/*
 * OOFEM generalized closest-point projection data.
 */

#ifndef contactprojection_h
#define contactprojection_h

#include "floatarray.h"
#include <limits>

namespace oofem {

/**
 * Dimension of the master feature carrying a closest-point projection.
 *
 * Values intentionally equal the feature dimension. This makes restart data
 * stable and gives deterministic tie-breaking in favor of the smoother,
 * higher-dimensional feature at an exact common boundary.
 */
enum class ContactFeatureType : int
{
    Vertex = 0,
    Edge = 1,
    Surface = 2
};

/**
 * Geometry returned by the generalized closest-point procedure.
 *
 * The representative element/local coordinates retain ordinary FE
 * interpolation for force assembly. featureIndex identifies the local edge or
 * vertex inside that representative element; it is zero for a surface
 * projection.
 */
struct ContactProjection
{
    bool valid = false;
    ContactFeatureType featureType = ContactFeatureType::Surface;
    int elementId = -1;
    int featureIndex = 0;
    FloatArray localCoordinates;
    double gap = std::numeric_limits<double>::infinity();
    double distanceSquared = std::numeric_limits<double>::infinity();
    FloatArray normal;
    FloatArray tangent1;
    FloatArray tangent2;
};

} // namespace oofem

#endif // contactprojection_h
