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

#ifndef iga_h
#define iga_h

#include "inputrecord.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "mathfem.h"
#include "feinterpol.h"

///@name Input fields for IGAElement
//@{
#define _IFT_IGAElement_KnotSpanParallelMode "knotspanparmode"
//@}

namespace oofem {
#ifdef __OOFEG
class StructuralElementEvaluator;
void drawIGAPatchDeformedGeometry(Element *elem, StructuralElementEvaluator *se, oofegGraphicContext &gc, UnknownType);
#endif

/**
 * Geometry wrapper for IGA elements.
 */
class OOFEM_EXPORT FEIIGAElementGeometryWrapper : public FEICellGeometry
{
public:
    const IntArray *knotSpan;
    ElementGeometry *elemGeometry;
public:
    FEIIGAElementGeometryWrapper(ElementGeometry * _elemGeometry, const IntArray * _knotSpan) : FEICellGeometry() {
        this->elemGeometry = _elemGeometry;
        this->knotSpan = _knotSpan;
    }
    FEIIGAElementGeometryWrapper(ElementGeometry * _elemGeometry) : FEICellGeometry() {
        this->elemGeometry = _elemGeometry;
        this->knotSpan = NULL;
    }

    int giveNumberOfVertices() const { return elemGeometry->giveNumberOfNodes(); }
    const FloatArray *giveVertexCoordinates(int i) const { return elemGeometry->giveNode(i)->giveCoordinates(); }
};


/**
 * IntegrationElement represent nonzero knot span, derived from Integration Rule.
 */
class OOFEM_EXPORT IGAIntegrationElement : public GaussIntegrationRule
{
protected:
    IntArray knotSpan;     // knot_span(nsd)
public:
    IGAIntegrationElement(int _n, ElementGeometry * _e, IntArray & _knotSpan) :
        GaussIntegrationRule(_n, _e, 0, 0, false),
        knotSpan(_knotSpan) { }
    const IntArray *giveKnotSpan() { return & this->knotSpan; }
    void setKnotSpan(IntArray &src) { this->knotSpan = src; }
};


/**
 * Implements base IGAElement, supposed to be a parent class of all elements with B-spline or NURBS based interpolation.
 */
class OOFEM_EXPORT IGAElement : public Element
{
protected:
    // FEInterpolation interpolation;
#ifdef __PARALLEL_MODE
    IntArray knotSpanParallelMode;
#endif
public:
    IGAElement(int n, Domain * aDomain) : Element(n, aDomain) { }
    IRResultType initializeFrom(InputRecord *ir);

#ifdef __PARALLEL_MODE
    elementParallelMode giveKnotSpanParallelMode(int) const;
#endif

#ifdef __OOFEG
    virtual void  drawRawGeometry(oofegGraphicContext &mode);
#endif

protected:
    virtual int giveNsd() = 0; // this info is available also from interpolation. Do we need it here ???
};


/**
 * IGATSplineElement setups integration rules differently from IGAElement.
 */
class OOFEM_EXPORT IGATSplineElement : public IGAElement
{
public:
    IGATSplineElement(int n, Domain * aDomain) : IGAElement(n, aDomain) { }
    IRResultType initializeFrom(InputRecord *ir);

protected:
    virtual int giveNsd() = 0;
};
} // end namespace oofem
#endif //iga_h
