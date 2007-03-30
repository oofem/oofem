/* $Header: $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

#ifndef levelsetpcs_h
#include "materialinterface.h"
#include "mathfem.h"
#include "geotoolbox.h"

class LevelSetPCS;

/**
   Element interface for LevelSetPCS class representing level-set like material interface.
   The elements should provide specific functionality in order to colaborate with LEPlic and this 
   required functionality is declared in this interface.
 */
class LevelSetPCSElementInterface : public Interface
{
 protected:
 public:
  LevelSetPCSElementInterface () {}

  /** Evaluetes F in level set equation of the form
      fi_t+F(grad(fi), x)*norm(grad(fi)) = 0
      where for interface position driven by flow with speed u: 
      F=dotProduct(u,grad(fi))/norm(grad(fi))
  */
  virtual double LS_PCS_computeF (LevelSetPCS*, TimeStep*) = 0;

  /** Returns gradient of shape functions (assumed constatnt <- linear approx)
   */
  virtual void LS_PCS_computedN (FloatMatrix& answer) = 0;
  /// Returns receiver's volume
  virtual double LS_PCS_computeVolume() = 0;

  /** Evaluetes S in level set equation of the form
      fi_t = S(fi)*(1-norm(grad(fi))) = 0
      where for interface position driven by flow with speed u: 
      S=fi/sqrt(fi^2+eps^2)
  */
  virtual double LS_PCS_computeS (LevelSetPCS*, TimeStep*) = 0;

  /**
     Returns VOF fractions for each material on element 
     according to nodal values of level set function (passed as parameter)
   */
  virtual void LS_PCS_computeVOFFractions (FloatArray& answer, FloatArray& fi) = 0;
};

/**
   Abstract base class representing Level Set representation of material interfaces. 
   The solution algorithm is based on positive coefficint scheme.
   This algorithm is limitted to d-dimensional simplexes with linear approximation.
   Its typical use to model moving interface (such as free surface)
   in a fixed-grid methods (as typically used in CFD).
   The basic tasks are representation of interface and its updating.
 */
class LevelSetPCS : public MaterialInterface
{
 protected:
  /// array used to store value of level set function for each node
  FloatArray levelSetValues, previousLevelSetValues;
  enum PCSEqType {PCS_levelSetUpdate, PCS_levelSetRedistance};
  Polygon initialRefMatVol;
  bool initialRefMatFlag;
  
 public:
  /** Constructor. Takes two two arguments. Creates 
      MaterialInterface instance with given number and belonging to given domain.
      @param n component number in particular domain. For instance, can represent  
      node number in particular domain.
      @param d domain to which component belongs to 
  */
  LevelSetPCS (int n,Domain* d) : MaterialInterface (n,d) { initialRefMatFlag = false;}

  /// initialize receiver
  virtual void initialize ();

  /**
     Updates the position of interface according to state reached in given solution step.
  */
  virtual void updatePosition (TimeStep* atTime);
  /**
     Updates element state after equlibrium in time step has been reached.
     All temporary history variables,
     which now describe equlibrium state should be copied into equilibrium ones.
     The existing internal state is used for update.
  */
  virtual void updateYourself (TimeStep* tStep) {previousLevelSetValues=levelSetValues;}
  double computeCriticalTimeStep (TimeStep*);
  virtual IRResultType initializeFrom (InputRecord* ir) ;

  virtual void redistance (TimeStep* atTime);


  /**
     Returns relative material contens at given point. Usually only one material is presented in given point,
     but some smoothing may be applied close to material interface to make transition smooth 
  */
  virtual void giveMaterialMixtureAt (FloatArray& answer, FloatArray& position);
  /**
     Returns volumetric (or other based measure) of relative material contens in given element.
  */
  virtual void giveElementMaterialMixture (FloatArray& answer, int ielem);
  /** Returns scalar value representation of material Interface at given point. For visualization */
  virtual double giveNodalScalarRepresentation (int i) {return levelSetValues.at(i);}

  /// Returns level set value in specific node
  double giveLevelSetDofManValue (int i) {return levelSetValues.at(i);}

  // identification
  const char* giveClassName () const { return "LevelSetPCS";}
  classType giveClassID ()      const { return LevelSetPCSClass;}


 protected:

  void pcs_stage1 (FloatArray& fs, FloatArray& w, TimeStep* atTime, PCSEqType t);
  double evalElemFContribution (PCSEqType t, int ie, TimeStep* atTime) ;
  double evalElemfContribution (PCSEqType t, int ie, TimeStep* atTime) ;
};

#define levelsetpcs_h
#endif
