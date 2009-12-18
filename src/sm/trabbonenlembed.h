/* $Header: /home/cvs/bp/oofem/sm/src/idmnl1.h,v 1.9 2003/04/06 14:08:30 bp Exp $ */
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

//   ********************************************
//   *** CLASS NONLOCAL TRABECULAR BONE MODEL ***
//   ********************************************

#ifndef trabbonenlembed_h 
 
#include "trabboneembed.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"

#include "sparsemtrx.h"
#include "dynalist.h"

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "conTable.h"
#endif

namespace oofem {

class GaussPoint ;


/////////////////////////////////////////////////////////////////
/////////TRABECULAR BONE NONLOCAL MATERIAL STATUS////////////////
/////////////////////////////////////////////////////////////////


class TrabBoneNLEmbedStatus : public TrabBoneEmbedStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
  protected:
 
  /////////////////////////////////////////////////////////////////
  // STATE VARIABLE DECLARATION
  //Equivalent strain for avaraging 

  double localCumPlastStrainForAverage;

  public:
 
  /////////////////////////////////////////////////////////////////
  // CONSTRUCTOR
  //

  TrabBoneNLEmbedStatus(int n, Domain*d, GaussPoint* g);


  /////////////////////////////////////////////////////////////////
  // DESTRUCTOR
  // 

  ~TrabBoneNLEmbedStatus ();


  /////////////////////////////////////////////////////////////////
  // OUTPUT PRINT 
  // Prints the receiver state to stream

  void   printOutputAt (FILE *file, TimeStep* tStep) ;


  /////////////////////////////////////////////////////////////////
  // STATE VARIABLE 
  // declare state variable access and modification methods

  double giveLocalCumPlastStrainForAverage ()     {return localCumPlastStrainForAverage;}
  void   setLocalCumPlastStrainForAverage (double ls) {localCumPlastStrainForAverage = ls;}

  
  ///////////////////////////////////////////////////////////////// 
  // DEFINITION
  // 

  const char* giveClassName () const { return "TrabBoneNLEmbedStatus" ;}
  classType   giveClassID () const { return TrabBoneEmbedStatusClass; }


  /////////////////////////////////////////////////////////////////
  // INITIALISATION OF TEMPORARY VARIABLES
  // Initializes the temporary internal variables, describing the current state according to
  // previously reached equilibrium internal variables.

  virtual void initTempStatus ();


  /////////////////////////////////////////////////////////////////
  // UPDATE VARIABLES 
  // Update equilibrium history variables according to temp-variables.
  // Invoked, after new equilibrium state has been reached.

  virtual void updateYourself(TimeStep*); // update after new equilibrium state reached


  /////////////////////////////////////////////////////////////////
  // SAVE CONTEXT - INTERUPT/RESTART
  // saves current context(state) into stream
  // Only non-temp internal history variables are stored.
  // @param stream stream where to write data
  // @param obj pointer to integration point, which invokes this method
  // @return contextIOResultType.
 
  contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL);


  /////////////////////////////////////////////////////////////////
  // RESTORE CONTEXT - INTERUPT/RESTART
  // Restores context of receiver from given stream. 
  // @param stream stream where to read data
  // @param obj pointer to integration point, which invokes this method
  // @return contextIOResultType.

  contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL);

  /////////////////////////////////////////////////////////////////
  // INTERFACE
  //Interface requesting service.
  //In the case of nonlocal constitutive models, 
  //the use of multiple inheritance is assumed. Typically, the class representing nonlocal 
  //constitutive model status is derived both from class representing local status and from class 
  //NonlocalMaterialStatusExtensionInterface or from one of its derived classes 
  //(which declare services and variables corresponding to specific analysis type).
  //@return In both cases, this function returns pointer to this object, obtained by 
  //returning adress of component or using pointer conversion from receiver to base class 
  //NonlocalMaterialStatusExtensionInterface. 

  virtual Interface* giveInterface (InterfaceType) ;
 
}; 


/////////////////////////////////////////////////////////////////
////////////TRABECULAR BONE NONLOCAL MATERIAL////////////////////
/////////////////////////////////////////////////////////////////


class TrabBoneNLEmbed : public TrabBoneEmbed, public StructuralNonlocalMaterialExtensionInterface
{
  protected :

  /////////////////////////////////////////////////////////////////
  // STATE VARIABLE DECLARATION
  // declare material properties here

  double R;
  double mParam;

  public :
 
  /////////////////////////////////////////////////////////////////
  // CONSTRUCTOR
  //

  TrabBoneNLEmbed(int n,Domain* d) ;

  /////////////////////////////////////////////////////////////////
  // DESTRUCTOR
  // 

  ~TrabBoneNLEmbed () ;

  /////////////////////////////////////////////////////////////////
  // INITIALIZATION OF FUNCTION/SUBROUTINE
  // 


  const char* giveClassName () const { return "TrabBoneNLEmbed" ;}
  classType   giveClassID ()   const { return TrabBoneEmbedClass;}
  const char*    giveInputRecordName() const {return "trabbonenlembed";}
 
  IRResultType initializeFrom (InputRecord* ir);

  virtual int giveInputRecordString(std::string &str, bool keyword = true);

  virtual Interface* giveInterface (InterfaceType);

  virtual void computeCumPlastStrain (double& alpha, GaussPoint* gp, TimeStep* atTime) ;

  void giveRealStressVector(FloatArray& answer,  MatResponseForm form, GaussPoint* gp , const FloatArray& strainVector, TimeStep* atTime);

  void computeLocalCumPlastStrain (double& alpha, const StrainVector& strain, GaussPoint* gp, TimeStep* atTime)
  {
    TrabBoneEmbed::computeCumPlastStrain (alpha, gp, atTime);
  }

  //@param CumPlastStrain equivalent strain vector in given integration point.
  //@param gp integration point to update.
  //@param atTime solution step indicating time of update.
  
  virtual void updateBeforeNonlocAverage (const FloatArray& strainVector, GaussPoint* gp, TimeStep* atTime) ;

  //Computes the value of nonlocal weight function in given point. 
  //@param src coordinates of source point.
  //@param coord coordinates of point, where nonlocal weight function is evaluated.
  //@return value of weight function.
  
  virtual double computeWeightFunction (const FloatArray& src, const FloatArray& coord) ;

  //Determines, whether receiver has bounded weighting function (limited support)
  //@return true if weighting function bounded, zero otherwise
  
  virtual int hasBoundedSupport () {return 1;}
 
  //Determines the width (radius) of limited support of weighting function
  
  virtual void giveSupportRadius (double& radius) {radius = this-> R;}

  //Computes the damage parameter from given equivalent strain in given integration point.
  
  protected: 

  // Creates the corresponding material status

  MaterialStatus* CreateStatus (GaussPoint *gp) const {return new TrabBoneNLEmbedStatus (1,TrabBoneEmbed::domain,gp);}

} ;

} // end namespace oofem
#define trabonenl3d_h
#endif
