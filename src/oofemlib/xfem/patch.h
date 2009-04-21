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

/** abstract representation of a part of element after subdivision  */
class Patch : public BasicGeometry
{
protected:
    /// parental element
    Element *parent;
    // here also Material will be added

public:
    Patch(Element *parent);
    /// converts the GP into the parental system of an element
    virtual void convertGPIntoParental(GaussPoint *gp) = 0;
};

class TrianglePatch : public Patch
{
public:
    TrianglePatch(Element *e) : Patch(e) {}
    // interpolation
    static FEI2dTrLin interpolation;
    void convertGPIntoParental(GaussPoint *gp);
};

#endif

