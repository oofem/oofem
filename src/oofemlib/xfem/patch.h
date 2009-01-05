/* 
 * File:   subpatchintegrationrule.h
 * Author: chamrova
 *
 * Created on November 2, 2008, 1:27 PM
 */

#ifndef _SUBPATCHINTEGRATIONRULE_H
#define	_SUBPATCHINTEGRATIONRULE_H

#include "delaunay.h"
#include "gausspnt.h"
#include "element.h"
#include "fei2dtrlin.h"

/* so far Patch done in order to implement the change of IntegrationRule,
 Patch inherits from Element, some functions of Element are redundant,
 the abstraction of the class is not solved yet, Patch is supposed to be a
 2d triangular patch with PlaneStress, this can be adjusted in the function
 setIntegrationRule */

class Patch : public Element{
    protected:
        // Delaunay Triangle
        Triangle *shape;
        // parental element
        Element *parent;
        // interpolation
        static FEI2dTrLin interpolation;
    public:
        /* constructor; destructor not needed since Triangle* is destructed in XfemManager,
         * since Patches share their vertices; Element is destructed in the Domain
        */
        Patch(int n, Domain *aDomain, Triangle *shape, Element* parent);
        // sets the Integration Rule
        void setIntegrationRule();
        // computes the volume around a GP
        double computeVolumeAround (GaussPoint* gp);
        // computes the total patch volume
        double computeVolume();
        // wrapper fro convertGPIntoParental
        void convertIntoParental();
        // converts the GP into the parental system of an element
        void convertGPIntoParental(GaussPoint *gp);
        // getter
        Triangle * getShape() {return shape;}
        
};

#endif

