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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#ifndef contactpoint_h
#define contactpoint_h

#include "gausspoint.h"
#include "fecontactsurface.h"
#include "feinterpol.h"

namespace oofem {
class IntArray;
class FloatArray;

/**
 * Abstract base class for all contact finite elements. Derived classes should be  base
 * classes for specific analysis type (for example base class for structural analysis,
 * thermal analysis or magnetostatics one). These derived classes then declare
 * analysis-specific part of interface and they provide default implementation
 * for these methods.
 * This abstract class declares (and possibly implements) general data and methods
 * common to all element types. General methods for obtaining characteristic vectors,
 * matrices and values are introduced and should be used instead of calling directly
 * specific member functions (these must be overloaded by derived analysis-specific
 * classes in order to invoke proper method according to type of component requested).
 */

class ContactPoint
{
protected:
  int surface_dimension;
public:
  ContactPoint() {}
  ~ContactPoint(){;}
  virtual void computeNmatrix(FloatMatrix &answer) = 0;
  virtual void compute_dNdxi_matrix(FloatMatrix &Bs, int i) = 0;
  //  
  virtual void computeCurvature(FloatMatrix &G, const FloatArray &normal, TimeStep *tStep) = 0;
  //
  virtual const FloatArray &giveLocalCoordinates() = 0;
  virtual FloatArray giveGlobalCoordinates() = 0;
  //
  virtual bool giveLocationArray(IntArray &locationArray,const IntArray &dofIDArry, 
                           const UnknownNumberingScheme &s) const = 0;
  virtual void giveUnknownVector(FloatArray &answer, const IntArray &dofMask, ValueModeType mode, TimeStep *tStep, bool padding = false) = 0;
  virtual FloatArray giveNormalVector() = 0;
  virtual void updateYourself(TimeStep *tStep){;}
  virtual void init(){;}
  virtual bool inContact() = 0;
  virtual void computeVectorOf(ValueModeType u, TimeStep *tStep, FloatArray &answer) = 0;
  virtual void giveUpdatedCoordinates(FloatArray &coords, TimeStep* tStep) = 0;
  int giveSurfaceDimension(){return surface_dimension;}

};


class FEContactPoint : public ContactPoint
{
protected:
  int contactElementId;
  std::unique_ptr<FEContactSurface> contactSurface;

public:
  FEContactPoint(FEContactSurface *cs, int ceId, int sd) : ContactPoint(), contactElementId(ceId), contactSurface(std::move(cs)){this->surface_dimension  = sd;}
  ~FEContactPoint(){;}
  //
  void computeNmatrix(FloatMatrix &answer) override;
  //  
  void compute_dNdxi_matrix(FloatMatrix &Bs, int i) override;
  //
  FloatArray giveNormalVector() override;
  //
  void computeVectorOf(ValueModeType mode, TimeStep *tStep, FloatArray &answer) override;
  //  
  void computeCurvature(FloatMatrix &G, const FloatArray &normal, TimeStep *tStep) override;
  //
  bool giveLocationArray(IntArray &locationArray,const IntArray &dofIDArry, 
			 const UnknownNumberingScheme &s) const override;
  //
  void giveUnknownVector(FloatArray &answer, const IntArray &dofMask, ValueModeType mode, TimeStep *tStep, bool padding = false) override;
  void giveUpdatedCoordinates(FloatArray &coords, TimeStep* tStep) override;   
  ///////////////////////////////////////////////////////////////////////////////////////
  bool inContact() override {return(contactElementId < 0 ? false : true);}
  //
  const FloatArray &giveLocalCoordinates() override = 0;
  FloatArray giveGlobalCoordinates() override = 0;
  //
   FEInterpolation* giveInterpolation();
  int giveContactElementId(){return contactElementId;}
  void setContactElementId(int ceId){contactElementId = ceId;}
  /////////////////////////////////////////////////////////////////////////////////////////
};

  
class FEContactPoint_Slave : public FEContactPoint
{
protected:
  std::unique_ptr<GaussPoint> slave_point;
  
public:
  FEContactPoint_Slave(FEContactSurface *cs, int ceId, int sd,GaussPoint *gp)  : FEContactPoint(cs, ceId, sd),slave_point(std::move(gp)){;}
  ~FEContactPoint_Slave(){;}
  //
  const FloatArray &giveLocalCoordinates() override {return slave_point->giveNaturalCoordinates();}
  FloatArray giveGlobalCoordinates() override {return slave_point->giveGlobalCoordinates();}
};

  
  
class FEContactPoint_Master : public FEContactPoint
{
protected:
  FloatArray localCoordinates;  
  
public:
  FEContactPoint_Master(FEContactSurface *cs, int ceId, int sd, FloatArray lc)  : FEContactPoint(cs, ceId, sd),localCoordinates(lc){;}
  ~FEContactPoint_Master(){;}
  //
  const FloatArray &giveLocalCoordinates() override {return this->localCoordinates;}
  FloatArray giveGlobalCoordinates() override;
};

  
} // end namespace oofem
#endif //contactpoint_h
