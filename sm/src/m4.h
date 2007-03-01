/* $Header: /home/cvs/bp/oofem/sm/src/m4.h,v 1.4.4.1 2004/04/05 15:19:47 bp Exp $ */
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

  //*************************************************************************
  //*** ABSTRACT CLASS MICROPLANE MATERIAL ACCORDING TO BAZANTS APPROACH ***
  //*************************************************************************

#ifndef M4Material__h 

#include "microplanematerial_bazant.h"
#include "structuralms.h"
  /**
     Related material model status to BazantMBCMaterial class 
     for storing history  variables in particular integration point 
     (microplane). Unique copy for each microplane must exist.
  */

  class M4MaterialStatus : public StructuralMaterialStatus
{

 protected:
  // add history variables declaration here

      

 public:
  M4MaterialStatus (int n, Domain*d, GaussPoint* g) ;
  ~M4MaterialStatus ();

  // add declaration of access and update functions for history variables

  // definition
  /// Returns class name of receiver
  const char* giveClassName () const { return "M4MaterialStatus" ;}
  /// Returns class ID of receiver
  classType             giveClassID () const { return M4MaterialStatusClass; }

  /**
     Initializes temporary history variables of receiver to previous equilibrium values.
  */
  virtual void initTempStatus ();
  /**
     Updates history variables after new equilibrium state has been reached.
     (e.g. temporary variables are copied into equilibrium values
  */
  virtual void updateYourself(TimeStep*);

  /// Stores receiver's state to stream
  contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL);
  /// Restores receiver's state from stream
  contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL);
}; 


/**
   Implementation of microplane material model according to Bazant's boundary curve
   approach.
*/
class M4Material : public MicroplaneMaterial_Bazant
{

 protected:

  double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;
  double c13,c14,c15,c16,c17,c18,c19,c20; /*c... fixed empirical constants*/
  double k1,k2,k3,k4,k5,E,nu,mu;
  double EV,ED,ET;
  double talpha;

 public:
 
  /**
     Constructor. Creates  Bazant's Boundary Curve Microplane Material belonging 
     to domain d, with number n.
     @param n material number
     @param d domain to which newly created material belongs
  */
  M4Material (int n,Domain* d) ;
  /// Destructor.
  ~M4Material () {}


  /**
     Computes characteristic matrix of receiver.
  */
  virtual void  giveCharacteristicMatrix (FloatMatrix& answer,
					  MatResponseForm form,
					  MatResponseMode mode,
					  GaussPoint* gp,
					  TimeStep* atTime);
  /**
     Returns a vector of coefficients of thermal dilatation in direction
     of each material principal (local) axis.
     @param answer vector of thermal dilatation coefficients
     @param gp integration point
     @param tStep time step (most models are able to respond only when atTime is current time step)
  */
  void giveThermalDilatationVector (FloatArray &answer, GaussPoint*, TimeStep*);

  /**
     Computes real stress vector (the meaning of  values depends on particular implementaion,
     e.g, can contain volumetric, devatoric normal srtresses and shear streses on microplane)
     for given increment of microplane strains.
     @param answer computed result
     @param mplane pointer to microplane object, for which response is computed
     @param strain strain vector
     @param tStep time step
  */
  virtual void giveRealMicroplaneStressVector (FloatArray& answer, Microplane* mplane, const FloatArray& strain, TimeStep* tStep);

  /**
     Tests, if material supports material mode.
     @param mode required material mode
     @return nonzero if supported, zero otherwise
  */

  double macbra(double x);
  double FVplus (double ev,double k1,double c13,double c14,double c15,double Ev);
  double FVminus (double ev,double k1,double k3,double k4,double E);
  double FDminus(double ed,double k1,double c7,double c8,double c9,double E);
  double FDplus(double ed,double k1,double c5,double c6,double c7,double c20,double E);
  double FN(double en,double sv,double k1,double c1,double c2,double c3,double c4,
	    double E, double Ev);
  double FT (double sn,double ev,double k1,double k2,double c10,
	     double c11, double c12, double Et);

  void updateVolumetricStressTo (Microplane* mPlane, double sigv) ;

  virtual int  giveSizeOfReducedStressStrainVector (MaterialMode);
  virtual int hasMaterialModeCapability (MaterialMode mode) ;

  /// Instanciates receiver from input record.
  IRResultType initializeFrom (InputRecord* ir);
  /// Returns class name of the receiver.
  const char*    giveClassName () const { return "M4Material" ;}
  /// Returns classType id of receiver.
  classType giveClassID ()         const {return M4MaterialClass;}

 protected: 
 
  MaterialStatus* CreateMicroplaneStatus (GaussPoint *gp) {return new M4MaterialStatus (1,domain,gp);}

};

#define M4Material__h 
#endif
