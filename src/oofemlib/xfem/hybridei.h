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

#ifndef HYBRIDEI_H_
#define HYBRIDEI_H_

#define _IFT_HybridEI_Name "hybridei"

#include "xfem/geometrybasedei.h"

namespace oofem {
class XfemManager;
class Domain;

/**
 * EnrichmentItem with hybrid geometry description in the following sense:
 * We have a BasicGeometry to describe the underlying geometry, and we use this
 * BasicGeometry to compute nodal level set fields. The enrichment is based on
 * interpolation of these nodal fields.
 * @author Erik Svenning
 * @date Sep 9, 2014
 */
class OOFEM_EXPORT HybridEI : public GeometryBasedEI
{
public:
    HybridEI(int n, XfemManager *xm, Domain *aDomain);
    virtual ~HybridEI();

    virtual const char *giveClassName() const { return "HybridEI"; }
    virtual const char *giveInputRecordName() const { return _IFT_HybridEI_Name; }

    virtual void evalLevelSetNormal(double &oLevelSet, const FloatArray &iGlobalCoord, const FloatArray &iN, const IntArray &iNodeInd) const;
    virtual void evalLevelSetTangential(double &oLevelSet, const FloatArray &iGlobalCoord, const FloatArray &iN, const IntArray &iNodeInd) const;
    virtual void evalGradLevelSetNormal(FloatArray &oGradLevelSet, const FloatArray &iGlobalCoord, const FloatMatrix &idNdX, const IntArray &iNodeInd) const;

    // By templating the function this way, we may choose if we want to pass iNodeInd as
    // an IntArray, a std::vector<int> or something else.
    // Any container that contains int and implements [] is legal.
    //    template< typename T >
    void interpLevelSet(double &oLevelSet, const FloatArray &iN, const IntArray &iNodeInd) const;

    void interpLevelSetTangential(double &oLevelSet, const FloatArray &iN, const IntArray &iNodeInd) const;

    void interpGradLevelSet(FloatArray &oGradLevelSet, const FloatMatrix &idNdX, const IntArray &iNodeInd) const;
};
} /* namespace oofem */

#endif /* HYBRIDEI_H_ */
