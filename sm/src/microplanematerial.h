/* $Header: /home/cvs/bp/oofem/sm/src/microplanematerial.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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

//   *****************************************
//   *** CLASS GENERAL MICROPLANE MATERIAL ***
//   *****************************************

#ifndef microplanematerial_h 

#include "structuralmaterial.h"
#include "gausspnt.h"
#include "matconst.h"

class Microplane;

#define MAX_NUMBER_OF_MICROPLANES 61

/**
 Abstract base class for all microplane models. 

 This class provides methods for setting up the integration weights and normals 
 for particular microplanes. Also projection tensors are computed in advance 
 (see method initializeData) and stored. 
 Methods for computing macro strain projections onto particular microplanes and
 for homogenization of stress vector are provided.
*/
class MicroplaneMaterial : public StructuralMaterial
{
protected:
 /// Number of microplanes;
 int numberOfMicroplanes;

 /// Integration weights of microplanes
 double microplaneWeights [MAX_NUMBER_OF_MICROPLANES];
 /// Normals of microplanes.
 double microplaneNormals [MAX_NUMBER_OF_MICROPLANES][3];

/** kronecker's delta */
  double Kronecker [6];

  /** 
  Normal projection tensors for all microplanes. 
  Due to symmetry, compressed form is stored.
  */
  double N [MAX_NUMBER_OF_MICROPLANES][6];
  /**
  Shear projection tensors (m direction) for all microplanes.
  Due to symmetry, compressed form is stored.
  */
  double M [MAX_NUMBER_OF_MICROPLANES][6];
  /**
  Shear projection tensors (l direction) for all microplanes.
  Due to symmetry, compressed form is stored.
  */
  double L [MAX_NUMBER_OF_MICROPLANES][6];

public: 
 
 /**
  Constructor. Creates Microplane Material belonging to domain d, with number n.
  @param n material number
  @param d domain to which newly created material belongs
  */
  MicroplaneMaterial (int n,Domain* d) : StructuralMaterial (n,d) {numberOfMicroplanes=0;}
  /// Destructor.
  ~MicroplaneMaterial ()                {}
 
  /**
  Computes real stress vector on given microplane
  (the meaning of  values depends on particular implementaion,
  e.g, can contain volumetric, devatoric normal srtresses and shear streses on microplane)
  for given increment of microplane strains.
  @param answer computed result
  @param mplane pointer to microplane object, for which response is computed
  @param strain strain vector 
  @param tStep time step
  */
  virtual void giveRealMicroplaneStressVector (FloatArray& answer, Microplane* mplane, 
                        const FloatArray& strain, TimeStep* tStep) = 0;

/**
 Computes the length of normal strain vector on given microplane.
 */
double computeNormalStrainComponent (Microplane* mplane, const FloatArray& macroStrain);
/**
 Computes the normal volumetric component of macro strain on given microplane.
*/
double computeNormalVolumetricStrainComponent (Microplane* mplane, const FloatArray& macroStrain);
/**
 Computes the normal deviatoric component of macro strain on given microplane.
 */
double computeNormalDeviatoricStrainComponent (Microplane* mplane, const FloatArray& macroStrain);
/**
 Computes the shear component (in m direction) of macro strain on given microplane.
 */
double computeShearMStrainComponent (Microplane* mplane, const FloatArray& macroStrain);
/**
 Computes the shear component (in l direction) of macro strain on given microplane.
 */
double computeShearLStrainComponent (Microplane* mplane, const FloatArray& macroStrain);
/**
 Computes the vector of all micro stress components (Ev, En, Em, El) of macro strain 
 vector on given microplane. 
 */
void computeStrainVectorComponents (FloatArray& answer, Microplane* mplane, 
                  const FloatArray& macroStrain);


  /**
  Computes normal of given microplane.
  @param answer normal of given microplane
  @mplane microplane, which normal will be computed
  */
  virtual void giveMicroplaneNormal (FloatArray& answer, Microplane* mplane);
  /**
  Returns microplane integration weight.
  @param mplane microplane
  @return integration weight of given microplane
  */
  virtual double giveMicroplaneIntegrationWeight (Microplane* mplane);

 /// Returns i-th microplane belonging to master-macro-integration point. )-based indexing.
  Microplane* giveMicroplane (int i, GaussPoint* masterGp);

  /** 
  Initializes internal data (integration weights, 
  microplane normals and computes projection tensors)
  @param numberOfMicroplanes number of required microplanes
  */
  virtual void initializeData (int numberOfMicroplanes);

/**
 Returns corresponding material mode for microplane according to macro integration mode
 */
  virtual MaterialMode giveCorrespondingSlaveMaterialMode (MaterialMode masterMode) ;

  // store & restore context functions
 /**
  Stores context of receiver into given stream. This method is called from 
  integration point saveContext function, to store material related staus in 
  integration point. Integration point passes itself as obj parameter, when
  invokes this method. This implementation loops over all slaves ip (microplanes) 
  and saves their satuses as wel as status of master is saved. 
  @param stream stream where to write data
  @param obj pointer to integration point, which invokes this method
  @return contextIOResultType.
  */
  contextIOResultType    saveContext (FILE* stream, void *obj = NULL);
 /**
  Restores context of receiver from given stream. This method is called from 
  integration point restoeContext function, to restore material related staus in 
  integration point. Integration point passes itself as obj parameter, when
  invokes this method. This implementation loops over all slaves ip (microplanes) 
  and loads their satuses as wel as status of master is loaded. 
  @param stream stream where to read data
  @param obj pointer to integration point, which invokes this method
  @return contextIOResultType.
  */
  contextIOResultType    restoreContext(FILE* stream, void *obj = NULL);

  /// Instanciates receiver from input record.
 IRResultType initializeFrom (InputRecord* ir);

 // identification and auxiliary functions
/// Returns class name of the receiver.
 const char*    giveClassName () const { return "MicroplaneMaterial" ;}
/// Returns classType id of receiver.
 classType giveClassID ()         const {return MicroplaneMaterialClass;}

    virtual MaterialStatus* giveMicroplaneStatus (GaussPoint* gp);
protected:
    virtual MaterialStatus* CreateMicroplaneStatus (GaussPoint* gp) = 0;
    virtual void initTempStatus (GaussPoint* gp);
} ;


#define microplanematerial_h
#endif
