/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.h,v 1.3 2003/04/23 14:22:15 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2002   Borek Patzak                                       



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

#ifndef tr1_2d_supg_axi_h
#define tr1_2d_supg_axi_h


#include "tr1_2d_supg.h"

/**
  Class representing 2d linear axisymmetric triangular element
  for solving incompressible fluid with SUPG solver

  This class is similar to TR1_2D_SUPG2_AXI, but diference is in handling 
  multiple fluids. This class uses rule of mixture which interpolates the
  properties using VOF value, requiring the use of twofluidmaterial class
  as material model for this situation.
 */
class TR1_2D_SUPG_AXI : public TR1_2D_SUPG
{
protected:
  /// radius at element center
  double rc;
public:
  // constructor
  TR1_2D_SUPG_AXI (int,Domain*) ;
  ~TR1_2D_SUPG_AXI () ;                        // destructor
  
  /**
     Computes acceleration terms (generalized mass matrix with stabilization terms ) for momentum balance equations(s)
   */
  void computeAccelerationTerm_MB (FloatMatrix& answer, TimeStep* atTime);
  /**
     Computes nonlinear advection terms for momentum balance equations(s)
  */
  void computeAdvectionTerm_MB (FloatArray& answer, TimeStep* atTime); 
  /**
     Computes the derivative of advection terms for momentum balance equations(s)
     with respect to nodal velocities
  */
  void computeAdvectionDerivativeTerm_MB (FloatMatrix &answer, TimeStep* atTime) ;
  /** 
      Computes diffusion terms for momentum balance equations(s)
  */
  void computeDiffusionTerm_MB (FloatArray& answer, TimeStep* atTime);
  /** Computes the derivative of diffusion terms for momentum balance equations(s)
      with respect to nodal velocities
  */
  void computeDiffusionDerivativeTerm_MB (FloatMatrix& answer, MatResponseMode mode, TimeStep* atTime);
  /** Computes pressure terms for momentum balance equations(s) */
  void computePressureTerm_MB (FloatMatrix& answer, TimeStep* atTime);
  /** Computes SLIC stabilization term for momentum balance equation(s) */
  void computeLSICStabilizationTerm_MB (FloatMatrix& answer, TimeStep* atTime);
  /** Computes the linear advection term for mass conservation equation */
  void computeLinearAdvectionTerm_MC (FloatMatrix& answer, TimeStep* atTime);
  /**
     Computes advection terms for mass conservation equation 
  */
  void computeAdvectionTerm_MC (FloatArray& answer, TimeStep* atTime);
  /** Computes the derivative of advection terms for mass conservation equation 
      with respect to nodal velocities
  */
  void computeAdvectionDerivativeTerm_MC (FloatMatrix& answer, TimeStep* atTime);
  /**
     Computes diffusion terms for mass conservation equation 
  */
  void computeDiffusionDerivativeTerm_MC (FloatMatrix& answer, TimeStep* atTime);
  void computeDiffusionTerm_MC (FloatArray& answer, TimeStep* atTime) ;
  /**
     Computes acceleration terms for mass conservation equation 
  */
  void  computeAccelerationTerm_MC (FloatMatrix& answer, TimeStep* atTime);
  /**
     Computes pressure terms for mass conservation equation 
  */
  void computePressureTerm_MC (FloatMatrix& answer, TimeStep* atTime);

  void     updateStabilizationCoeffs (TimeStep* );
  // calculates critical time step
  // virtual double        computeCriticalTimeStep (TimeStep* tStep);
  /**
     Computes Rhs terms due to boundary conditions
  */
  void  computeBCRhsTerm_MB (FloatArray& answer, TimeStep* atTime);
  /**
     Computes Rhs terms due to boundary conditions
  */
  void  computeBCRhsTerm_MC (FloatArray& answer, TimeStep* atTime);

  double             computeVolumeAround (GaussPoint*) ;

  // definition
  const char* giveClassName () const { return "TR1_2D_SUPG_AXI" ;}
  classType                giveClassID () const { return SUPGElementClass; }
  
 protected:
 void                  computeGaussPoints ();
 virtual void computeDeviatoricStress (FloatArray& answer, GaussPoint* gp, TimeStep* );
 virtual void initGeometry            ();
 double computeRadiusAt (GaussPoint*);
 void   computeBMtrx (FloatMatrix& answer, GaussPoint* gp);
 void   computeNVector (FloatArray& answer, GaussPoint* gp);

} ;

#endif // tr1_2d_supg_axi_h







