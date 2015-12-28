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

#ifndef GEOMETRYBASEDEI_H_
#define GEOMETRYBASEDEI_H_

#define _IFT_GeometryBasedEI_Name "geometrybasedei"

#include "enrichmentitem.h"

#include <memory>

namespace oofem {
class XfemManager;
class Domain;

/**
 * EnrichmentItem with geometry described by BasicGeometry
 * @author Erik Svenning
 * @date Sep 9, 2014
 */
class OOFEM_EXPORT GeometryBasedEI : public EnrichmentItem
{
public:
    GeometryBasedEI(int n, XfemManager *xm, Domain *aDomain);
    virtual ~GeometryBasedEI();

    virtual int instanciateYourself(DataReader *dr);


    virtual void appendInputRecords(DynamicDataReader &oDR);

    virtual const char *giveClassName() const { return "GeometryBasedEI"; }
    virtual const char *giveInputRecordName() const { return _IFT_GeometryBasedEI_Name; }

    virtual void updateGeometry();
    virtual void updateNodeEnrMarker(XfemManager &ixFemMan);

    void updateLevelSets(XfemManager &ixFemMan);

    virtual void evaluateEnrFuncInNode(std :: vector< double > &oEnrFunc, const Node &iNode) const;

    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl) const;
    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl, const FloatArray &iN, const IntArray &iElNodes) const;

    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl) const;
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl, const FloatArray &iN, const FloatMatrix &idNdX, const IntArray &iElNodes) const;

    // TODO: Consider moving this function to a separate Cohesive Zone Interface /ES
    void evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, int iNodeInd, GaussPoint &iGP, bool iGPLivesOnCurrentCrack) const;


    virtual void computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, std :: vector< int > &oIntersectedEdgeInd, Element *element, std :: vector< double > &oMinDistArcPos) const;
    virtual void computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, std :: vector< int > &oIntersectedEdgeInd, Element *element, const Triangle &iTri, std :: vector< double > &oMinDistArcPos) const;

    virtual void writeVtkDebug() const;

    void giveSubPolygon(std :: vector< FloatArray > &oPoints, const double &iXiStart, const double &iXiEnd) const;

    virtual void propagateFronts(bool &oFrontsHavePropagated);
    virtual bool giveElementTipCoord(FloatArray &oCoord, double &oArcPos, Element &iEl, const FloatArray &iElCenter) const;

    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius);

    BasicGeometry *giveGeometry() { return mpBasicGeometry.get(); }

protected:
    std :: unique_ptr< BasicGeometry > mpBasicGeometry;
};
} /* namespace oofem */

#endif /* GEOMETRYBASEDEI_H_ */
