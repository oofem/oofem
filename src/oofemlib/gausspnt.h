/* $Header: /home/cvs/bp/oofem/oofemlib/src/gausspnt.h,v 1.11.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */


//   *************************
//   *** CLASS GAUSS POINT ***
//   *************************

#ifndef gausspnt_h
#define gausspnt_h

#include "flotarry.h"
#include "integrationrule.h"
#include "matstatus.h"
#include "tdictionary.h"
#include "element.h"

class Material;
class LayeredCrossSection;
class MicroplaneMaterial;
class FiberedCrossSection;
class CrossSction;
class IntegrationRule;

/**
 * Class representing integration point in finite element program.
 * Integration point maintain its space position, integration
 * weight and corresponding material mode.
 * Link to related element which integration point belongs to
 * is also maintained.
 * Integration point generally can contain list of slave integration points
 * therefore is called as master point. Slaves are used for example to implement
 * layered or fibred cross sections by cross section class. Then in one
 * "macro" master gauss point, cross section creates few slaves (one per layer)
 * and puts them into master list. When cross sections completes requests for
 * particular master integration point, it performs integration over layers.
 * It therefore calls material class for each layer, sending corresponding
 * slave as parameter and integrates results.
 * Generally, every integretion point must hold
 * its own copy of history variables (which are related to corresponding
 * material model used). These material type dependent history variables
 * are stored in material type related material status, which can be
 * managed by integration point.
 * Each material model class should introduce related material status class
 * (derived from material status class or from its children), where necesary
 * history variables are kept and can be accesed by material.
 * Material class then creates unique copy of related status in all necessary
 * integration points. Because integration point is compulsory parameter of
 * all member functions of material class, particular material then can
 * easily access its associated status from integration point and therefore its
 * history variables for particular integration point.
 *
 * To provide support to integrate element contribution by parts from element subVolumes
 * (posibly with different material parameters etc), the integration point can maintain not
 * only its element natural coordinates, but also its subVolume local coordinates, that are neceessary
 * to compute its jacobian, for examle. These coordinates are stored in localCoordinates attribute.
 *
 */
class GaussPoint
{
    /*
     * This class implements a point for a gaussian quadrature. A Gauss point is
     * usually attribute of an element.
     * DESCRIPTION
     * A Gauss point is identified by its 'number' and by 'element' - the element
     * it belongs to. The array 'coordinates' defines its position in an axis sys-
     * tem defined by its element. 'weight' is the pondaration of the point in the
     * numerical integration ; 'weight' is naturally positive and smaller than 2.
     * stress and srain vectors  store the current values of strain and stresess ;
     * the size of these
     * arrays and the nature of their values is determined by the element .
     *
     * materialStatusDict stores corresponding MaterialStatuses (Base MaterialStatus)
     * as required by materials, yieldconditions,... They store their corresponding
     * instances of MaterialStatus Class with key = (int)obj->giveClassID(), where
     * obj is instance of material, yield crit which is associated with this gp.
     *
     * TASKS
     * - returning its coordinates, its weight and the strains and stresses ;
     * - reading/writing its attribute in a file.
     * REMARK
     * The Gauss point is a rather passive object : it does not compute its strains
     * and stresses - it just stores them. They are computed by the element. If ever
     * materially nonlinearity is implemented, the role of the Gauss point is
     * likely to grow.
     */

private:
    /// Number
    int number;
    /// Reference to parent integration rule.
    IntegrationRule *irule;
    /// Natural Element Coordinates of receiver.
    FloatArray *coordinates;
    /// Optional local subpatch (subpatches form element volume) coordinates of the receiver
    FloatArray *localCoordinates;
    /// Integraation weight.
    double weight;
    //    FloatArray*  strainVector ;  // full form
    //    FloatArray*  stressVector ;  // full form
    //    FloatArray*  stressIncrementVector;// increments are used mainly in nonlinear analysis
    //    FloatArray*  strainIncrementVector;// to find balanced state. for update call gp->update()
    //    FloatArray*  plasticStrainVector; // full form
    //    FloatArray*  plasticStrainIncrementVector;
    /// Material mode of receiver.
    MaterialMode materialMode;

protected:

    // layer and fibred material support
    /// Number of slaves.
    int numberOfGp;
    /// List of slave integration points.
    GaussPoint **gaussPointArray;
    /// Material status managed for material model.
    MaterialStatus *matStatus;

public:
    /**
     * Constructor.  Creates integration point belonging to given integration rule,
     * with given number, integration weight, coordinates and material mode.
     * @param ir integration rule to which integration point belongs to.
     * @param n integration point number
     * @param a coordinates
     * @param w integration weight
     * @param mode material mode
     */
    GaussPoint(IntegrationRule *ir, int n, FloatArray *a, double w, MaterialMode mode); // constructor
    /// Destructor
    virtual ~GaussPoint();                                  // destructor

    /// Returns i-th natural element coordinate of receiver
    double       giveCoordinate(int i)      { return coordinates->at(i); }
    /// Returns coordinate array of receiver.
    FloatArray *giveCoordinates()          { return coordinates; }
    void         setCoordinates(const FloatArray &c) { * coordinates = c; }

    /// Returns local subpatch coordinates of the receiver
    FloatArray *giveLocalCoordinates() { if ( localCoordinates ) { return localCoordinates; } else { return coordinates; } }
    void         setLocalCoordinates(const FloatArray &c)
    { if ( localCoordinates ) { * localCoordinates = c; } else { localCoordinates = c.GiveCopy(); } }


    //      FloatArray*  giveStrainVector ()         { return strainVector ;}
    //      FloatArray*  giveStressVector ()         { return stressVector ;}
    //    FloatArray*  givePlasticStrainVector()   { return plasticStrainVector;}
    //    FloatArray*  giveStrainIncrementVector() { return strainIncrementVector;}
    //    FloatArray*  giveStressIncrementVector() { return stressIncrementVector;}
    //    FloatArray*  givePlasticStrainIncrementVector() {return plasticStrainIncrementVector;}
    //    FloatArray*  giveHardeningParam ()       { return NULL; }

    /// Returns  integration weight odf receiver.
    virtual double       giveWeight()               { return weight; }
    void                 setWeight(double w) { weight = w; }
    /// Returns number of receiver.
    int          giveNumber()               { return number; }
    /// Returns corresponding integration rule to receiver
    IntegrationRule *giveIntegrationRule()   { return irule; }
    /// Returns corresponding element to receiver.
    Element *giveElement()               { return irule->giveElement(); }
    /// Sets element of gp
    //void         setElement (Element* e)     { element=e;}
    /// Returns corresponding material mode of receiver.
    MaterialMode giveMaterialMode()          { return this->materialMode; }
    /// Sets material mode of receiver.
    void         setMaterialMode(MaterialMode newMode) { this->materialMode = newMode; }
    /// Returns reference to material associated to related element of receiver.
    Material *giveMaterial()             { return giveElement()->giveMaterial(); }
    /// Returns reference to cross section associated to related element of receiver.
    CrossSection *giveCrossSection()        { return giveElement()->giveCrossSection(); }
    /// Returns reference to associated material status (NULL if not defined)
    MaterialStatus *giveMaterialStatus()    { return matStatus; }
    /**
     * Sets Material status managed by receiver.
     * Old status, if exist will be lost.
     * @param ptr poiter to new status of receiver.
     * @return pointer to new status.
     */
    MaterialStatus *setMaterialStatus(MaterialStatus *ptr)
    { delete matStatus;
      return ( matStatus = ptr ); }
    /**
     * Returns index-th slave gauss point of receiver.
     * @param index of retuned slave
     */
    GaussPoint *giveSlaveGaussPoint(int index);
    //      void         letStrainVectorBe (FloatArray* v)
    //   { delete strainVector ; strainVector = v ;}
    //      void         letStressVectorBe (FloatArray* v)
    //   { delete stressVector ; stressVector = v ;}
    //    void         letPlasticStrainVectorBe (FloatArray* v)
    //   { delete plasticStrainVector; plasticStrainVector = v;}
    //    void         letStressIncrementVectorBe (FloatArray* v)
    //   { delete stressIncrementVector ; stressIncrementVector = v ;}
    //    void         letStrainIncrementVectorBe (FloatArray* v)
    //   { delete strainIncrementVector ; strainIncrementVector = v ;}
    //    void         letPlasticStrainIncrementVectorBe (FloatArray* v)
    //   { delete plasticStrainIncrementVector; plasticStrainIncrementVector = v;}
    /**
     * Prints output of receiver to file. Corresponding printOutputAt  function for
     * associated status is called. The same fuction is also invoked for all available
     * slaves of receiver.
     */
    virtual void         printOutputAt(FILE *, TimeStep *);
    /**
     * Updates internal state of receiver after finishing time step.
     * Material::updateYourself (receiver, tStep) function is called to
     * update material status. Same fuction is also invoked for
     * all receiver's slaves.
     */
    virtual void         updateYourself(TimeStep *);
    //    void         resetYourself () ;

    // store & restore context functions
    /*
     * Stores receiver state to output stream, including associated material status.
     * Warning: Slaves are not saved, they must be saved by corresponding
     * cross section or material model, which creates these slaves.
     * This is because they may have special weights and coordinates,
     * and these are not saved into context file, because they remain the same
     * during whole solution (only variables which vary are stored into context).
     * Note: does not invoke FEMComponents saveContext, since this wrires only
     * class id header, but typically due to large number of IPs,
     * this is avoided.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    //contextIOResultType    saveContext (FILE* stream, void *obj = NULL);
    /*
     * Restores receiver state to output stream, including associated material status.
     * Warning: Slaves are not restored, they must be restored by corresponding
     * cross section or material model, which creates these slaves.
     * This is because they may have special weights and coordinates,
     * and these are not saved into context file, because they remain the same
     * during whole solution (only variables which vary are stored into context).
     * Note: does not invoke FEMComponents restoreContext, since this wrires only
     * class id header, but typically due to large number of IPs,
     * this is avoided.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    //contextIOResultType    restoreContext(FILE* stream, void *obj = NULL);
    /// Returns classType id of receiver.
    virtual classType giveClassID() const { return GaussPointClass; }
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "GaussPoint"; }
    ///Initializes receiver acording to object description stored in initString.
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

    friend class LayeredCrossSection;
    friend class MicroplaneMaterial;
    friend class FiberedCrossSection;
    /*   friend Material;
     * friend LayeredMaterial ; */
};


#endif // gausspnt_h





