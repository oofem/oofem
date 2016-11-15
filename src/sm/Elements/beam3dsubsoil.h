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
 *               Copyright (C) 1993 - 2014   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef beam3dsubsoil_h
#define beam3dsubsoil_h

#include "../sm/Elements/structuralelement.h"
#include "interface.h"


#define _IFT_Beam3dSubsoil_Name "beam3dsubsoil"
#define  _IFT_Beam3dSubsoil_springConstants "springconstants"
#define _IFT_Beam3dSubsoil_localFormulation "localformulation"
#define _IFT_Beam3dSubsoil_master "master"

namespace oofem {
class FEI3dLineLin;
  /**
   * This class implements an subsoil element for any beam3d (qubic, linear, etc)
   * element implementing Beam3dSubsoilElementIterface.  
   * The soil is assumed to act in z-direction.
   *
   * Modes supported: Winkler (spring support)
   *                  Winkler-Pasternak (not yet)
   *
   * From principle of Wirtual work the contribution of Winkler model is K=\int N^t K N dx,
   * where N is the displacement interpolation and K is diagonal matrix of spring coefficients.
   *
   * Reference:
   * 
   * @author Borek Patzak
   */
  class OOFEM_EXPORT Beam3dSubsoil : public StructuralElement {
    protected:
      /// Geometry interpolator only.
      static FEI3dLineLin interp;
      /// Enum to define type of subsoil model 
      enum modeEnum {_Winkler, _WinklerPasternak};
      /// Type of subsoil model
      modeEnum mode;
      /** array of spring (subsoil) coefficients
	  for Winkler model for each DOF one elastic constant
      **/
      FloatArray springConstants;
      /// Flag indicating whether subsoil model defined in global or element local c.s.
      bool local_formulation;
      /// Master element (from which interpolation is obtained)
      int master;
      
    public:
      Beam3dSubsoil(int n, Domain * d);
      virtual ~Beam3dSubsoil() { }
      
      virtual FEInterpolation *giveInterpolation() const; 
      virtual MaterialMode giveMaterialMode()  { return _Unknown; }
      virtual int testElementExtension(ElementExtension ext) { return 0;}
      
      // definition & identification
      virtual const char *giveInputRecordName() const { return _IFT_Beam3dSubsoil_Name; }
      virtual const char *giveClassName() const { return "Beam3dSubsoil"; }
      virtual IRResultType initializeFrom(InputRecord *ir);
      
      virtual int computeNumberOfDofs() { return 12; }
      virtual void giveDofManDofIDMask(int inode, IntArray &) const;
      
      virtual double computeVolumeAround(GaussPoint *gp);
      
      virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
      virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
      { computeLumpedMassMatrix(answer, tStep); }
      
      virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
      void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
      
    protected:
      virtual void computeGaussPoints();
      virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
      virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
      virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
      
  };
  
  /**
     Interface defining required functionality of master element. Mostly, an interpolation 
     matrix and transformation matrices are to be provided in order to abstract from particular
     element type and interpolation.
   */
  class OOFEM_EXPORT Beam3dSubsoilElementInterface : public Interface
  {
  public:
    /// Constructor
    Beam3dSubsoilElementInterface() ;

    /** Evaluate element interpolation matrix N, where $u^l = Nr^l$.
	The interpolation matrix should contain 6 rows corresponding to 6 expected DOFs 
	(3 displacements and three rotations)
    */
    virtual void B3SSI_getNMatrix (FloatMatrix &answer, GaussPoint *gp) = 0;
    /// Evaluate transformation matrix for reciver DOFs
    virtual void B3SSI_getGtoLRotationMatrix (FloatMatrix &answer) = 0;
    /// Evalueates transfromation matrix for element unknowns
    virtual void B3SSI_getNodalGtoLRotationMatrix(FloatMatrix& answer) = 0;
    /// Computes Jacobian at given integration point
    virtual double B3SSI_computeVolumeAround (GaussPoint* gp) = 0;
  };
  
  
} // end namespace oofem
#endif // beam3dsubsoil_h
