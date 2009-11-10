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

#ifndef trabbonenl3d_h 
 
#include "trabbone3d.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"

#include "sparsemtrx.h"
#include "dynalist.h"

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "conTable.h"
#endif

class GaussPoint ;


/////////////////////////////////////////////////////////////////
/////////TRABECULAR BONE NONLOCAL MATERIAL STATUS////////////////
/////////////////////////////////////////////////////////////////


class TrabBoneNL3DStatus : public TrabBone3DStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
  protected:
  // STATE VARIABLE DECLARATION
  // Equivalent strain for avaraging 
  double localCumPlastStrainForAverage;

  public:
  // CONSTRUCTOR
  TrabBoneNL3DStatus(int n, Domain*d, GaussPoint* g);

  // DESTRUCTOR
  ~TrabBoneNL3DStatus ();

  // OUTPUT PRINT 
  // Prints the receiver state to stream
  void   printOutputAt (FILE *file, TimeStep* tStep) ;

  // STATE VARIABLE 
  // declare state variable access and modification methods
  double giveLocalCumPlastStrainForAverage ()     {return localCumPlastStrainForAverage;}
  const FloatArray *giveLTangentContrib();
  void setLocalCumPlastStrainForAverage (double ls) {localCumPlastStrainForAverage = ls;}

  // DEFINITION
  const char* giveClassName () const { return "TrabBoneNL3DStatus" ;}
  classType   giveClassID () const { return TrabBone3DStatusClass; }

  // INITIALISATION OF TEMPORARY VARIABLES
  // Initializes the temporary internal variables, describing the current state according to
  // previously reached equilibrium internal variables.
  virtual void initTempStatus ();

  // UPDATE VARIABLES 
  // Update equilibrium history variables according to temp-variables.
  // Invoked, after new equilibrium state has been reached.
  virtual void updateYourself(TimeStep*); // update after new equilibrium state reached

  // INTERFACE
  // Interface requesting service.
  // In the case of nonlocal constitutive models, 
  // the use of multiple inheritance is assumed. Typically, the class representing nonlocal 
  // constitutive model status is derived both from class representing local status and from class 
  // NonlocalMaterialStatusExtensionInterface or from one of its derived classes 
  // (which declare services and variables corresponding to specific analysis type).
  // @return In both cases, this function returns pointer to this object, obtained by 
  // returning adress of component or using pointer conversion from receiver to base class 
  // NonlocalMaterialStatusExtensionInterface. 
  virtual Interface* giveInterface (InterfaceType) ;
 
}; 


/////////////////////////////////////////////////////////////////
////////////TRABECULAR BONE NONLOCAL MATERIAL////////////////////
/////////////////////////////////////////////////////////////////


class TrabBoneNL3D : public TrabBone3D, public StructuralNonlocalMaterialExtensionInterface, 
public NonlocalMaterialStiffnessInterface
{
  protected :
  // STATE VARIABLE DECLARATION
  // declare material properties here
  double R;
  double mParam;

  public :
  // CONSTRUCTOR
  TrabBoneNL3D(int n,Domain* d) ;

  // DESTRUCTOR
  ~TrabBoneNL3D () ;

  // INITIALIZATION OF FUNCTION/SUBROUTINE
  const char* giveClassName () const { return "TrabBoneNL3D" ;}
  classType   giveClassID ()   const { return TrabBone3DClass;}
  const char*    giveInputRecordName() const {return "trabbonenl3d";}
 
  // Initializes the receiver from given record
  IRResultType initializeFrom (InputRecord* ir);

  // Returns input record name of the receiver.
  virtual int giveInputRecordString(std::string &str, bool keyword = true);

  // Interface requesting service
  virtual Interface* giveInterface (InterfaceType);

  // Computes the nonlocal cumulated plastic strain from its local form.
  // @param kappa return param, containing the corresponding cumulated plastic strain
  // @param gp integration point
  // @param atTime time step
  virtual void computeCumPlastStrain (double& kappa, GaussPoint* gp, TimeStep* atTime) ;

  // Computes the local cumulated plastic strain from given strain vector (full form).
  // @param kappa return param, comtaining the corresponding cumulated plastic strain
  // @param strain total strain vector in full form
  // @param gp integration point
  // @param atTime time step
  void computeLocalCumPlastStrain (double& kappa, const StrainVector& strain, GaussPoint* gp, TimeStep* atTime) 
  {
    TrabBone3D::computeCumPlastStrain (kappa, gp, atTime);
  }

  /////////////////////////////////////////////////////////////////
  // NONLOCAL STIFFNESS FUNCTION/SUBROUTINE

  virtual void give3dMaterialStiffnessMatrix (FloatMatrix& answer,  MatResponseForm, MatResponseMode,GaussPoint* gp,  TimeStep* atTime);

  #ifdef __OOFEG
  // Plots the sparse structure of stiffness contribution.
    virtual void NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(GaussPoint* gp, oofegGraphicContext& gc, TimeStep* atTime);
  #endif

  // @name Services required by NonlocalMaterialStiffnessInterface and related ones to support Nonlocal Stiffness*/
  // @{compute ans add IP contributions to destination matrix
  virtual void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx& dest, GaussPoint* gp, TimeStep* atTime);

  // Returns integration list of receiver. Contains localIntegrationRecord structures, containing 
  // references to integration points and their weights that influence to nonlocal average in 
  // receiver's associated integration point.
  virtual dynaList<localIntegrationRecord>* NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint* gp);

  // Computes the "local" part of nonlocal stiffness contribution assembled for given integration point.
  // @param gp source integration point
  // @param loc local code numbers
  // @param lcontrib "local" contribution
  // @return nonzero if local point contributes (loading) or zero if not (unloading in elastic range, elastic)
  int giveLocalNonlocalStiffnessContribution (GaussPoint* gp, IntArray& loc, FloatArray& lcontrib, TimeStep* atTime);

  // Computes the "remote" part of nonlocal stiffness contribution assembled for given integration point.
  // @param gp remote integration point
  // @param loc remote element code numbers
  // @param rcontrib "remote" contribution
  void giveRemoteNonlocalStiffnessContribution (GaussPoint* gp, IntArray& rloc, FloatArray& rcontrib, TimeStep* atTime);

  // Computes elastic stiffness for normal stress components
  // @param answer result of size (3,3)
  // @param mode determines the MatResponseMode
  // @param gp integration point
  // @param atTime time step
//  void giveNormalElasticStiffnessMatrix (FloatMatrix& answer, MatResponseMode rMode, GaussPoint*gp, TimeStep* atTime) ;

  // NONLOCAL STIFFNESS FUNCTION/SUBROUTINE
  /////////////////////////////////////////////////////////////////

  void giveRealStressVector(FloatArray& answer,  MatResponseForm form, GaussPoint* gp , const FloatArray& strainVector, TimeStep* atTime);

  // Implements the service updating local variables in given integration points, 
  // which take part in nonlocal average process. Actually, no update is necessary,
  // because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
  // computation. It is therefore necessary only to store local strain in corresponding status.
  // @param CumPlastStrain equivalent strain vector in given integration point.
  // @param gp integration point to update.
  // @param atTime solution step indicating time of update.
  virtual void updateBeforeNonlocAverage (const FloatArray& strainVector, GaussPoint* gp, TimeStep* atTime) ;

  // Computes the value of nonlocal weight function in given point. 
  // @param src coordinates of source point.
  // @param coord coordinates of point, where nonlocal weight function is evaluated.
  // @return value of weight function.
  virtual double computeWeightFunction (const FloatArray& src, const FloatArray& coord) ;

  // Determines, whether receiver has bounded weighting function (limited support)
  // @return true if weighting function bounded, zero otherwise
  virtual int hasBoundedSupport () {return 1;}
 
  // Determines the width (radius) of limited support of weighting function
  virtual void giveSupportRadius (double& radius) {radius = this-> R;}


  #ifdef __PARALLEL_MODE
  // Updates domain before nonloc average (using updateDomainBeforeNonlocAverage service)
  // to ensure, that the localStrainVectorForAverage variable is correctly updated and
  // pack this localStrainVectorForAverage into given buffer.
  // @see Material::packUnknowns for description.
  // @param buff communication buffer
  // @param stepN solution step
  // @param ip integration point
  int packUnknowns (CommunicationBuffer& buff, TimeStep* stepN, GaussPoint* ip);
  
  // Unpack localStrainVectorForAverage value from given buffer.
  // @see Material::unpackAndUpdateUnknowns service.
  // @param buff communication buffer
  // @param stepN solution step.
  // @param ip integration point
  int unpackAndUpdateUnknowns (CommunicationBuffer& buff, TimeStep* stepN, GaussPoint* ip);

  // Estimates the necessary pack size to hold all packed data of receiver.
  int estimatePackSize (CommunicationBuffer& buff, GaussPoint* ip);
  #endif


  protected: 
  // Creates the corresponding material status
  MaterialStatus* CreateStatus (GaussPoint *gp) const {return new TrabBoneNL3DStatus (1,TrabBone3D::domain,gp);}

};

#define trabonenl3d_h
#endif




