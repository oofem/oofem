/*
 * XFEMDebugTools.h
 *
 *  Created on: Jun 5, 2013
 *      Author: svennine
 */

#ifndef XFEMDEBUGTOOLS_H_
#define XFEMDEBUGTOOLS_H_

#include "oofemcfg.h"
#include "geometry.h"

#include <fstream>

namespace oofem {
class OOFEM_EXPORT XFEMDebugTools
{
public:
    XFEMDebugTools();
    virtual ~XFEMDebugTools();

    static void WriteTrianglesToVTK(const std :: string &iName, const std :: vector< Triangle > &iTriangles);
    static void WritePointsToVTK(const std :: string &iName, const std :: vector< FloatArray > &iPoints);

    static void WriteArrayToMatlab(const std :: string &iName, const std :: vector< double > &iX, const std :: vector< double > &iY);
    static void WriteArrayToGnuplot(const std :: string &iName, const std :: vector< double > &iX, const std :: vector< double > &iY);
};
} /* namespace oofem */
#endif /* XFEMDEBUGTOOLS_H_ */
