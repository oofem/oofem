/*
 * File:   subpatchintegrationrule.h
 * Author: chamrova
 *
 * Created on November 2, 2008, 1:27 PM
 */

#ifndef _PATCH_H
#define _PATCH_H

#include "delaunay.h"
#include "gausspnt.h"
#include "element.h"
#include "fei2dtrlin.h"
#include "node.h"

namespace oofem {

/** abstract representation of a part of element after subdivision  */
class Patch : public BasicGeometry
{
protected:
    /// parental element
    Element *parent;
    /// Material of the patch
    int material;
public:
    Patch(Element *parent, int material);
    Patch(Element *parent, AList< FloatArray > *vertices);
    virtual ~Patch() { }
    /// converts the GP into the parental system of an element
    virtual void convertGPIntoParental(GaussPoint *gp) = 0;
    int giveMaterial() { return this->material; }
    Element *giveParent() { return this->parent; }
#ifdef __OOFEG
    /// draw pure geometry
    void          draw(oofegGraphicContext &gc) { BasicGeometry :: draw(gc); }
    /// draw with vertex data
    virtual void  drawWD(oofegGraphicContext &gc, FloatArray &vd) { };
#endif
};

class TrianglePatch : public Patch
{
public:
    TrianglePatch(Element *parent, int material) : Patch(parent, material) { }
    TrianglePatch(Element *parent, AList< FloatArray > *vertices) : Patch(parent, vertices) { }
    ~TrianglePatch() { }
    // interpolation
    static FEI2dTrLin interpolation;
    void convertGPIntoParental(GaussPoint *gp);

#ifdef __OOFEG
    void          draw(oofegGraphicContext &gc);
    void  drawWD(oofegGraphicContext &gc, FloatArray &vd);
#endif
};

} // end namespace oofem
#endif
