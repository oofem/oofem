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

#ifndef patch_h
#define patch_h

#include "geometry.h"
#include "fei2dtrlin.h"

namespace oofem {
class Element;
class GaussPoint;
class Node;

/**
 * Abstract representation of a part of element after subdivision.
 * @author chamrova
 */
class Patch : public BasicGeometry
{
public:
    enum PatchType { PT_Unknown, PT_TrianglePatch };

protected:
    /// Parental element.
    Element *parent;
    /// Material of the patch.
    int material;

public:
    Patch(Element *parent);
    Patch(Element *parent, int material);
    Patch(Element *parent, AList< FloatArray > *vertices);
    virtual ~Patch() { }
    /// Converts the GP into the parental system of an element.
    virtual void convertGPIntoParental(GaussPoint *gp) = 0;
    /// Returns material id associated to receiver.
    int giveMaterial() { return this->material; }
    /// Returns reference to parent element.
    Element *giveParent() { return this->parent; }
    /// Returns patch type id of receiver.
    virtual PatchType givePatchType() { return PT_Unknown; }

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

#ifdef __OOFEG
    /// Draw pure geometry.
    void draw(oofegGraphicContext &gc) { BasicGeometry :: draw(gc); }
    /// Draw with vertex data.
    virtual void drawWD(oofegGraphicContext &gc, FloatArray &vd) { };
#endif
};

class TrianglePatch : public Patch
{
public:
    TrianglePatch(Element *parent) : Patch(parent) {}
    TrianglePatch(Element *parent, int material) : Patch(parent, material) { }
    TrianglePatch(Element *parent, AList< FloatArray > *vertices) : Patch(parent, vertices) { }
    virtual ~TrianglePatch() { }
    // interpolation
    static FEI2dTrLin interpolation;
    void convertGPIntoParental(GaussPoint *gp);
    /// Returns patch type id of receiver
    PatchType givePatchType() { return PT_TrianglePatch; }
#ifdef __OOFEG
    void draw(oofegGraphicContext &gc);
    void drawWD(oofegGraphicContext &gc, FloatArray &vd);
#endif
};
} // end namespace oofem
#endif // patch_h
