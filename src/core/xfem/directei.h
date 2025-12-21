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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef DIRECTEI_H_
#define DIRECTEI_H_

#define _IFT_DirectEI_Name "directei"

#include "xfem/geometrybasedei.h"

namespace oofem {
class XfemManager;
class Domain;

/**
 * EnrichmentItem with direct geometry description in the following sense:
 * We have a BasicGeometry to describe the underlying geometry, and we use this
 * BasicGeometry directly to compute level set fields.
 * @author Erik Svenning
 * @date Sep 10, 2014
 */
class DirectEI : public GeometryBasedEI
{
public:
    DirectEI(int n, XfemManager *xm, Domain *aDomain);
    virtual ~DirectEI();

    const char *giveClassName() const override { return "DirectEI"; }
    const char *giveInputRecordName() const override { return _IFT_DirectEI_Name; }

    void evalLevelSetNormal(double &oLevelSet, const FloatArray &iGlobalCoord, const FloatArray &iN, const IntArray &iNodeInd) const override;
    void evalLevelSetTangential(double &oLevelSet, const FloatArray &iGlobalCoord, const FloatArray &iN, const IntArray &iNodeInd) const override;
    void evalGradLevelSetNormal(FloatArray &oGradLevelSet, const FloatArray &iGlobalCoord, const FloatMatrix &idNdX, const IntArray &iNodeInd) const override;
};
} /* namespace oofem */

#endif /* DIRECTEI_H_ */
