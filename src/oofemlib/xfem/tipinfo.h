/*
 * tipinfo.h
 *
 *  Created on: Aug 9, 2013
 *      Author: svennine
 */

#ifndef TIPINFO_H_
#define TIPINFO_H_

#include "floatarray.h"

namespace oofem {
/**
 * TipInfo gathers useful information about a crack tip,
 * like its position and tangent direction. The purpose of the
 * struct is to keep the interface clean when evaluating
 * XFEM branch functions.
 *
 * @author Erik Svenning
 *  August 2013
 */
class TipInfo
{
public:
    TipInfo() { }
    ~TipInfo() { }

    FloatArray mGlobalCoord;
    double mArcPos;
    FloatArray mTangDir;
    FloatArray mNormalDir;
    int mTipIndex;
    int mEdgeIndex; /// Local number of which edge the crack enters the element (2d)
};

struct TipPropagation {
    TipPropagation() { }
    ~TipPropagation() { }

    int mTipIndex;
    FloatArray mPropagationDir;
    double mPropagationLength;
};
} // end namespace oofem


#endif /* TIPINFO_H_ */
