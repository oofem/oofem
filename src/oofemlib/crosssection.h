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

#ifndef crosssection_h
#define crosssection_h

#include "femcmpnn.h"
#include "gausspoint.h"
#include "materialmode.h"
#include "matresponsemode.h"
#include "material.h"
#include "internalstatetype.h"
#include "internalstatevaluetype.h"
#include "crosssectextension.h"

namespace oofem {

/// List of properties possibly stored in a cross section.
enum CrossSectionProperty {
    CS_Thickness=400,  ///< Thickness
    CS_Width,          ///< Width
    CS_BeamShearCoeff, ///< Shear coefficient of beam
    CS_Area,           ///< Area
    CS_InertiaMomentY, ///< Moment of inertia around y-axis
    CS_InertiaMomentZ, ///< Moment of inertia around z-axis
    CS_TorsionMomentX, ///< Moment of inertia around x-axis
    CS_SHEAR_AREA_Y,   ///< Shear area in y direction
    CS_SHEAR_AREA_Z,   ///< Shear area in z direction
    CS_TopZCoord,      ///< Top z coordinate
    CS_BottomZCoord,   ///< Bottom z coordinate
    CS_NumLayers,      ///< Number of layers that makes up the cross section
};

/**
 * Base abstract class representing cross section in finite element mesh.
 *
 * The main idea, why cross section has been introduced, is to hide all details
 * of cross section description from particular element. Generally elements
 * do not communicate directly with material but communicate through cross section interface,
 * which therefore can perform necessary integration (for example over layers of fibers).
 *
 * The cross section returns properties like thickness and area.
 *
 * The derived classes are supposed to be base cross section classes for particular
 * type of analysis. They should declare general interface methods necessary.
 *
 * In particular cross section implementation, where is necessary to perform integration
 * over cross section volume (over layers, fibers, ...) and therefore generally one must keep
 * complete load history in these integration points, the concept of master-slave integration
 * points should be used.
 *
 * Integration point generally can contain list of slave integration points
 * therefore is called as master point. Slaves are used for example to implement
 * layered or fibered cross sections by cross section class. Then in one
 * "macro" master Gauss point, cross section creates few slaves (one per layer)
 * and puts them into master list. When cross sections completes requests for
 * particular master integration point, it performs integration over layers.
 * It therefore calls material class for each layer, sending corresponding
 * slave as parameter and integrates results.
 * @see GaussPoint class for more detail.
 */
class CrossSection : public FEMComponent
{
protected:
    /**
     * Dictionary for storing cross section parameters (like dimensions).
     * More preferably, (due to slow access into dictionary values) one should use
     * corresponding variables declared inside class
     */
    Dictionary *propertyDictionary;

public:
    /**
     * Constructor. Creates cross section with number n belonging to domain d.
     * @param n Cross section number.
     * @param d Domain.
     */
    CrossSection(int n, Domain *d) : FEMComponent(n, d)
    { propertyDictionary = new Dictionary(); }
    /// Destructor.
    virtual ~CrossSection() { delete propertyDictionary; }

    /**
     * Returns the value of cross section property.
     * @param a Id of requested property.
     * @return Property value.
     */
    virtual double give(CrossSectionProperty a);
    /**
     * Check for symmetry of stiffness matrix.
     * Default implementation returns true.
     * It can be moved to base Cross section class in the future.
     * @param rMode Response mode of material.
     * @param mat Material index.
     * @return True if stiffness matrix of receiver is symmetric.
     */
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode, int mat);

    virtual void printYourself();

    /**
     * Returns nonzero, if receiver implements required extension.
     * @param ext Required extension.
     * @return Nonzero, if supported, zero otherwise.
     */
    virtual int testCrossSectionExtension(CrossSectExtension ext) { return 0; }

    /**
     * Returns the integration point corresponding value in Reduced form.
     * @param answer contain corresponding ip value, zero sized if not available
     * @param ip Integration point.
     * @param type Determines the type of internal variable.
     * @param atTime Time step.
     * @return Nonzero if o.k, zero otherwise.
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *atTime)
    { return ip->giveMaterial()->giveIPValue(answer, ip, type, atTime); }
    /**
     * Returns the corresponding integration point value size in Reduced form.
     * @param type Determines the type of internal variable.
     * @param ip Integration point to check.
     * @return Nonzero if o.k, zero otherwise.
     */
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *ip)
    { return ip->giveMaterial()->giveIPValueSize(type, ip); }
    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer Mask of full vector size, with components being the indexes to reduced form vectors.
     * @param type Determines the internal variable requested (physical meaning).
     * @param mmode Material mode.
     * @param mat Material to check.
     * @return Nonzero if o.k or error is generated for unknown mat mode.
     */
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode, Material *mat)
    { return mat->giveIntVarCompFullIndx(answer, type, mmode); }
    /**
     * Returns the type of internal variable (scalar, vector, tensor,...).
     * @param type Determines the type of internal variable.
     * @param mat Material to check.
     * @return Type of internal variable.
     */
    virtual InternalStateValueType giveIPValueType(InternalStateType type, Material *mat)
    { return mat->giveIPValueType(type); }

#ifdef __PARALLEL_MODE
    /**
     * Pack all necessary data of integration point (according to element parallel_mode)
     * into given communication buffer. The corresponding material model service for particular integration point
     * is invoked. The nature of packed data is material model dependent.
     * Typically, for material of "local" response (response depends only on integration point local state)
     * no data are exchanged. For "nonlocal" constitutive models the send/receive of local values which
     * undergo averaging is performed between local and corresponding remote elements.
     * @param buff Communication buffer.
     * @param stepN Solution step.
     * @param ip Integration point.
     * @return Nonzero if successful.
     */
    virtual int packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
    { return ip->giveMaterial()->packUnknowns(buff, stepN, ip); }
    /**
     * Unpack and updates all necessary data of given integration point (according to element parallel_mode)
     * into given communication buffer.
     * @see packUnknowns service.
     * @param buff Communication buffer.
     * @param stepN Solution step.
     * @param ip Integration point.
     * @return Nonzero if successful.
     */
    virtual int unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
    { return ip->giveMaterial()->unpackAndUpdateUnknowns(buff, stepN, ip); }
    /**
     * Estimates the necessary pack size to hold all packed data of receiver.
     * The corresponding material model  service is invoked. The
     * nature of packed data is typically material model dependent.
     * @param buff Communication buffer.
     * @param ip Integration point.
     * @return Estimate of pack size.
     */
    virtual int estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip)
    { return ip->giveMaterial()->estimatePackSize(buff, ip); }
    /**
     * Returns the weight representing relative computational cost of receiver
     * The reference cross section is integral model in plane stress.
     * Its weight is equal to 1.0
     * Default implementation computes average computational cost of material model
     * and multiplies it by cross section type weight (obtained by giveRelativeSelfComputationalCost())
     * The other cross section models should compare to this reference.
     * @param ip Integration point.
     * @return Prediction of the computational cost.
     */
    virtual double predictRelativeComputationalCost(GaussPoint *ip);
    /**
     * Returns the weight representing relative computational cost of receiver
     * The reference element is integral model in plane stress.
     * Its weight is equal to 1.0.
     * The other cross section models should compare to this reference.
     * @return Relative computational cost of self.
     */
    virtual double giveRelativeSelfComputationalCost() { return 1.0; }
    /**
     * @return Relative redistribution cost of the receiver.
     */
    virtual double predictRelativeRedistributionCost(GaussPoint *gp) { return 1.0; }
#endif

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "CrossSection"; }
    virtual classType giveClassID() const { return CrossSectionClass; }

    /**
     * Stores integration point state to output stream.
     * @param stream Output stream.
     * @param mode Determines amount of info required in stream (state, definition, ...).
     * @param gp integration point.
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);
    /**
     * Reads integration point state to output stream.
     * @param stream Output stream.
     * @param mode Determines amount of info required in stream (state, definition, ...).
     * @param gp integration point.
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);

    friend class Material;
};
} // end namespace oofem
#endif // crosssection_h

