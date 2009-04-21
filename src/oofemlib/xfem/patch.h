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
#include "geotoolbox.h"

/** abstract representation of a part of element after subdivision  */
class Patch : public BasicGeometry
{
protected:
    /// parental element
    Element *parent;
    Material *mat;
    AList<GaussPoint> * gps;

public:
    Patch(Element *parent);
    Patch(Element *parent, AList<FloatArray> *vertices);
    ~Patch();
    /// converts the GP into the parental system of an element
    virtual void convertGPIntoParental(GaussPoint *gp) = 0;
    void computeMaterial();
    void addGps(GaussPoint *gp);
    GaussPoint *giveGaussPoint(int n) {return gps->at(n);}
    bool hasGaussPoint(GaussPoint *gp);
    Material * giveMaterial() {return this->mat; }
};

class TrianglePatch : public Patch
{
public:
    TrianglePatch(Element *e) : Patch(e) {}
    TrianglePatch(Element *e, AList<FloatArray> *vertices) : Patch(e, vertices) {}
    // interpolation
    static FEI2dTrLin interpolation;
    void convertGPIntoParental(GaussPoint *gp);
};


#endif

