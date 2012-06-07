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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef iga_h
#define iga_h

#include "inputrecord.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "mathfem.h"
#include "feinterpol.h"

namespace oofem {
#ifdef __OOFEG
class StructuralElementEvaluator;
void drawIGAPatchDeformedGeometry(Element * elem, StructuralElementEvaluator * se, oofegGraphicContext & gc, UnknownType);
#endif

/**
 * Geometry wrapper for IGA elements.
 */
class FEIIGAElementGeometryWrapper : public FEICellGeometry
{
public:
    const IntArray *knotSpan;
    Element *elem;
public:
    FEIIGAElementGeometryWrapper(Element *_elem, const IntArray *_knotSpan) : FEICellGeometry() {
        this->elem = _elem;
        this->knotSpan = _knotSpan;
    }
    FEIIGAElementGeometryWrapper(Element *_elem) : FEICellGeometry() {
        this->elem = _elem;
        this->knotSpan = NULL;
    }

    int giveNumberOfVertices() const { return elem->giveNumberOfNodes(); }
    const FloatArray *giveVertexCoordinates(int i) const { return elem->giveNode(i)->giveCoordinates(); }
};


/**
 * IntegrationElement represent nonzero knot span, derived from Integration Rule.
 */
class IGAIntegrationElement : public GaussIntegrationRule
{
protected:
    IntArray knotSpan;     // knot_span(nsd)
public:
    IGAIntegrationElement(int _n, Element *_e, IntArray &_knotSpan) :
        GaussIntegrationRule(_n, _e, 0, 0, false),
        knotSpan(_knotSpan) { }
    const IntArray *giveKnotSpan() { return & this->knotSpan; }
    void setKnotSpan(IntArray &src) { this->knotSpan = src; }
};


/**
 * Implements base IGAElement, supposed to be a parent class of all elements with B-spline or NURBS based interpolation.
 */
class IGAElement : public Element
{
protected:
    // FEInterpolation interpolation;
#ifdef __PARALLEL_MODE
    IntArray knotSpanParallelMode;
#endif
public:
    IGAElement(int n, Domain *aDomain) : Element(n, aDomain) { }
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
class IGATSplineElement : public IGAElement
{
public:
    IGATSplineElement(int n, Domain *aDomain) : IGAElement(n, aDomain) { }
    IRResultType initializeFrom(InputRecord *ir);

protected:
    virtual int giveNsd() = 0;
};
} // end namespace oofem
#endif //iga_h
