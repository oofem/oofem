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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "dynamicinputrecord.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "tensor/tensor2.h"
#include "tensor/tensor4.h"

#pragma once


///@name Input fields for VolumetricEnergyInterface
//@{
#define _IFT_BaseHyperElasticMaterial_k "k"
#define _IFT_BaseHyperElasticMaterial "type"
//@}


namespace oofem {
class InputRecord;
/**
 * Abstract base class for hyperelastic materials
 * It provides calculation of first and second derivative of principal invariants of the left Cauchy green tensor with respect to deformation gradient
 * It provides calculation of volumetric energy(currently only a logarithmic form is supported)
 * All the calculations are based on the core/tensor/tensoralgebra.h file which uses FTensor library
 * More to be added gradually
 * @author Martin Horak
 **/

class OOFEM_EXPORT BaseHyperElasticMaterial
{
protected:
    /**
     *  Type characterizing the volumetric energy of the hyperelastic material
     * Type 0: Logarithmic \f$W(J) = 0.5 * K * (ln J)^2 \f$
     * More types can be added ...
     *
     * @author Martin Horak
     **/
    enum VolumetricEnergyType {
        VET_Logarithmic = 0
    };
    /// Bulk modulus
    double K = 0;
    /// volumetric energy type
    VolumetricEnergyType VET_Type = VET_Logarithmic;


public:

    BaseHyperElasticMaterial() {; }
    ~BaseHyperElasticMaterial() {; }

    /// initialization for the input file
    void initializeFrom(InputRecord &ir);

    /**
     * Compute the first invariant of Cauchy-Green deformation tensor C
     * @param Deformation gradient (second-order tensor in 3d)
     * @return Trace of Cauchy-Green deformation tensor C
     **/
    inline double compute_I1_C_from_F(const Tensor2_3d &F) const
    {
        return ( F(k_3, 0) * F(k_3, 0) + F(k_3, 1) * F(k_3, 1) + F(k_3, 2) * F(k_3, 2) );
    }

    /**
     * Compute the second invariant of Cauchy-Green deformation tensor C
     * @param Deformation gradient (second-order tensor in 3d)
     * @return The second invariant of Cauchy-Green deformation tensor C
     **/

    inline double compute_I2_C_from_F(const Tensor2_3d &F) const
    {
        return ( 0.5 * ( ( F(k_3, 0) * F(k_3, 0) + F(k_3, 1) * F(k_3, 1) + F(k_3, 2) * F(k_3, 2) ) * ( F(k_3, 0) * F(k_3, 0) + F(k_3, 1) * F(k_3, 1) + F(k_3, 2) * F(k_3, 2) ) - F(k_3, 0) * F(k_3, l_3) * F(m_3, l_3) * F(m_3, 0) - F(k_3, 1) * F(k_3, l_3) * F(m_3, l_3) * F(m_3, 1) - F(k_3, 2) * F(k_3, l_3) * F(m_3, l_3) * F(m_3, 2) ) );
    }


    /**
     * Compute the third invariant of Cauchy-Green deformation tensor C
     * @param Deformation gradient (second-order tensor in 3d)
     * @return Determinant of Cauchy-Green deformation tensor C
     **/
    inline double compute_I3_C_from_F(const Tensor2_3d &F) const
    {
        double J = F.compute_determinant();
        return ( J * J );
    }

    /**
     * Compute the first invariant of Cauchy-Green deformation tensor C
     * @param Cauchy-Green deformation tensor C(second-order tensor in 3d)
     * @return Trace of Cauchy-Green deformation tensor C
     **/
    inline double compute_I1_C_from_C(const Tensor2_3d &C) const
    {
        return C(i_3, i_3);
    }

    /**
     * Compute the second invariant of Cauchy-Green deformation tensor C
     * @param Cauchy-Green deformation tensor C(second-order tensor in 3d)
     * @return The second invariant of Cauchy-Green deformation tensor C
     **/
    inline double compute_I2_C_from_C(const Tensor2_3d &C) const
    {
        double trC = C(i_3, i_3);
        return ( 0.5 * ( trC * trC - ( C(k_3, 0) * C(k_3, 0) + C(k_3, 1) * C(k_3, 1) + C(k_3, 2) * C(k_3, 2) ) ) );
    }

    /**
     * Compute the third of Cauchy-Green deformation tensor C
     * @param Cauchy-Green deformation tensor C(second-order tensor in 3d)
     * @return Determinant of Cauchy-Green deformation tensor C
     **/
    inline double compute_I3_C_from_C(const Tensor2_3d &C) const
    {
        double J2 = C.compute_determinant();
        return J2;
    }


    /**
     * Compute the first derivative of the first invariant of the Cauchy-Green deformation tensor C wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return first derivative of the Cauchy-Green deformation tensor C wrt deformation gradient F(second-order tensor in 3d)
     **/
    inline Tensor2_3d compute_dI1_C_dF(const Tensor2_3d &F) const
    {
        Tensor2_3d dI1_C_dF;
        dI1_C_dF(i_3, j_3) = 2. * F(i_3, j_3);
        return dI1_C_dF;
    }


    /**
     * Compute the first derivative of the second invariant of the Cauchy-Green deformation tensor C wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return first derivative of the second invariant of the Cauchy-Green deformation tensor C wrt deformation gradient F(second-order tensor in 3d)
     **/
    inline Tensor2_3d compute_dI2_C_dF(const Tensor2_3d &F) const
    {
        Tensor2_3d dI2_C_dF;
        dI2_C_dF(i_3, j_3) =  2. * ( F(k_3, 0) * F(k_3, 0) + F(k_3, 1) * F(k_3, 1) + F(k_3, 2) * F(k_3, 2) ) * F(i_3, j_3) - 2. * F(i_3, m_3) * F(k_3, m_3) * F(k_3, j_3);
        return dI2_C_dF;
    }



    /**
     * Compute the first derivative of the third invariant of the Cauchy-Green deformation tensor C wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return first derivative of the third invariant of the Cauchy-Green deformation tensor C wrt deformation gradient F(second-order tensor in 3d)
     **/
    inline Tensor2_3d compute_dI3_C_dF(const Tensor2_3d &F) const
    {
        auto [ detF, cofF ] = F.compute_determinant_and_cofactor();
        Tensor2_3d dI3_C_dF;
        dI3_C_dF(i_3, j_3) = 2. * detF * cofF(i_3, j_3);
        return dI3_C_dF;
    }


    /**
     * Compute the first derivative of the Jacobian(determinant of  F) wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return the first derivative of the Jacobian(determinant of  F) wrt deformation gradient F (second-order tensor in 3d)
     **/
    inline Tensor2_3d compute_dJ_dF(const Tensor2_3d &F) const
    {
        return F.compute_cofactor();
    }


    /**
     * Compute the second derivative of the first invariant of the Cauchy-Green deformation tensor wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return the second derivative of the first invariant of the Cauchy-Green deformation tensor wrt deformation gradient F (fourth-order tensor in 3d)
     **/
    inline Tensor4_3d compute_d2I1_C_dF2(const Tensor2_3d &F) const
    {
        Tensor2_3d I(1.0, 0., 0., 0., 1.0, 0., 0., 0., 1.0);
        Tensor4_3d d2I1_C_dF2;
        d2I1_C_dF2(i_3, j_3, k_3, l_3) = 2. *  I(i_3, k_3) * I(j_3, l_3);
        return d2I1_C_dF2;
    }

    /**
     * Compute the second derivative of the second invariant of the Cauchy-Green deformation tensor wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return the second derivative of the second invariant of the Cauchy-Green deformation tensor wrt deformation gradient F (fourth-order tensor in 3d)
     **/
    inline Tensor4_3d compute_d2I2_C_dF2(const Tensor2_3d &F) const
    {
        Tensor2_3d I(1.0, 0., 0., 0., 1.0, 0., 0., 0., 1.0);
        Tensor4_3d d2I2_C_dF2;
        d2I2_C_dF2(i_3, j_3, k_3, l_3) = 4. * F(i_3, j_3) * F(k_3, l_3) - 2. * F(i_3, l_3) * F(k_3, j_3) +  2. * ( F(m_3, 0) * F(m_3, 0) + F(m_3, 1) * F(m_3, 1) + F(m_3, 2) * F(m_3, 2) ) * I(i_3, k_3) * I(j_3, l_3) - 2. * F(m_3, l_3) * F(m_3, j_3) * I(i_3, k_3) - 2. * F(i_3, m_3) * F(k_3, m_3) * I(j_3, l_3);
        return d2I2_C_dF2;
    }



    /**
     * Compute the second derivative of the third invariant of the Cauchy-Green deformation tensor wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return the second derivative of the third invariant of the Cauchy-Green deformation tensor wrt deformation gradient F (fourth-order tensor in 3d)
     **/
    inline Tensor4_3d compute_d2I3_C_dF2(const Tensor2_3d &F) const
    {
        auto [ detF, cofF ] = F.compute_determinant_and_cofactor();
        Tensor4_3d d2I3_C_dF2;
        d2I3_C_dF2(i_3, j_3, k_3, l_3) = 2. * cofF(i_3, j_3) * cofF(k_3, l_3) + 2 * detF * F.compute_tensor_cross_product()(i_3, j_3, k_3, l_3);
        return d2I3_C_dF2;
    }

    /**
     * Compute the first invariant of the deviatoric Cauchy-Green deformation tensor C
     * @param Deformation gradient (second-order tensor in 3d)
     * @return Trace of the deviatoric Cauchy-Green deformation tensor C
     **/
    inline double compute_I1_Cdev_from_F(const Tensor2_3d &F) const
    {
        return pow(F.compute_determinant(), -2. / 3.) * compute_I1_C_from_F(F);
    }

    /**
     * Compute the second invariant of the deviatoric Cauchy-Green deformation tensor C
     * @param Deformation gradient (second-order tensor in 3d)
     * @return Second invariant of the deviatoric Cauchy-Green deformation tensor C
     **/
    inline double compute_I2_Cdev_from_F(const Tensor2_3d &F) const
    {
        return pow(F.compute_determinant(), -4. / 3.) * compute_I2_C_from_F(F);
    }


    /**
     * Compute the first derivatiove of the first invariant of the deviatoric Cauchy-Green deformation tensor C wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return first derivatiove of the first invariant of the deviatoric Cauchy-Green deformation tensor C (second-order tensor)
     **/
    inline Tensor2_3d compute_dI1_Cdev_dF(const Tensor2_3d &F) const
    {
        Tensor2_3d dI1_Cdev_dF;
        auto [ J, cofF ] = F.compute_determinant_and_cofactor();
        dI1_Cdev_dF(i_3, j_3) = 2. * pow(J, -2. / 3.) * ( F(i_3, j_3) - compute_I1_C_from_F(F) / ( 3. * J ) * cofF(i_3, j_3) );
        return dI1_Cdev_dF;
    }

    /**
     * Compute the first derivatiove of the second invariant of the deviatoric Cauchy-Green deformation tensor C wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return first derivatiove of the second invariant of the deviatoric Cauchy-Green deformation tensor C (second-order tensor)
     **/
    inline Tensor2_3d compute_dI2_Cdev_dF(const Tensor2_3d &F) const
    {
        Tensor2_3d iF, dI2_Cdev_dF, C;
        auto [ J, cofF ] = F.compute_determinant_and_cofactor();
        C(i_3, j_3) = F(k_3, i_3) * F(k_3, j_3);
        auto I2 = compute_I2_C_from_C(C);

        dI2_Cdev_dF(i_3, j_3) = 2. * pow(F.compute_determinant(), -4. / 3.) * ( C(k_3, k_3) * F(i_3, j_3) - F(i_3, k_3) * C(k_3, j_3) - 2. * I2 / ( 3. * J ) * cofF(i_3, j_3) );
        return dI2_Cdev_dF;
    }

    /**
     * Compute the second derivatiove of the first invariant of the deviatoric Cauchy-Green deformation tensor C wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return second derivatiove of the first invariant of the deviatoric Cauchy-Green deformation tensor C (fourth-order tensor)
     **/
    inline Tensor4_3d compute_d2I1_Cdev_dF2(const Tensor2_3d &F) const
    {
        Tensor2_3d I(1.0, 0., 0., 0., 1.0, 0., 0., 0., 1.0);
        auto I1 = compute_I1_C_from_F(F);
        auto iF = F.compute_inverse();
        Tensor4_3d d2I1_Cdev_dF2;
        d2I1_Cdev_dF2(i_3, j_3, k_3, l_3) = 2. / 3. * pow(F.compute_determinant(), -2. / 3.) * ( 3. * I(i_3, k_3) * I(j_3, l_3) + I1 * iF(j_3, k_3) * iF(l_3, i_3) + 2. / 3. * I1 * iF(j_3, i_3) * iF(l_3, k_3) - 2. * iF(l_3, k_3) * F(i_3, j_3) - 2. * iF(j_3, i_3) * F(k_3, l_3) );
        return d2I1_Cdev_dF2;
    }

    /**
     * Compute the second derivatiove of the second invariant of the deviatoric Cauchy-Green deformation tensor C wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return Second derivative of the second invariant of the deviatoric Cauchy-Green deformation tensor C (second-order tensor)
     **/
    inline Tensor4_3d compute_d2I2_Cdev_dF2(const Tensor2_3d &F) const
    {
        Tensor2_3d I(1.0, 0., 0., 0., 1.0, 0., 0., 0., 1.0);
        auto [ J, cofF ] = F.compute_determinant_and_cofactor();
        auto J_43 = pow(J, -4. / 3.);
        Tensor2_3d iF, C, B, FC;
        iF(i_3, j_3) = cofF(j_3, i_3) / J;
        C(i_3, j_3) = F(k_3, i_3) * F(k_3, j_3);
        B(i_3, j_3) = F(i_3, k_3) * F(j_3, k_3);
        FC(i_3, j_3) = F(i_3, k_3) * C(k_3, j_3);
        auto I1 = compute_I1_C_from_C(C);
        auto I2 = compute_I2_C_from_C(C);
        Tensor4_3d d2I2_Cdev_dF2;
        d2I2_Cdev_dF2(i_3, j_3, k_3, l_3) = 2. * J_43 * ( I1 * I(i_3, k_3) * I(j_3, l_3) + 2. * F(i_3, j_3) * F(k_3, l_3) - 4. / 3. * I1 * F(i_3, j_3)  * iF(l_3, k_3) + 8. / 9. * I2 * iF(j_3, i_3)  * iF(l_3, k_3) - 4. / 3. * I1 * iF(j_3, i_3)  * F(k_3, l_3) + 4. / 3. * FC(k_3, l_3) *  iF(j_3, i_3) + 2. / 3. * I2 * iF(j_3, k_3)  * iF(l_3, i_3) + 4. / 3. * FC(i_3, j_3) * iF(l_3, k_3) - C(l_3, j_3) * I(i_3, k_3) - F(i_3, l_3) * F(k_3, j_3) - B(i_3, k_3) * I(l_3, j_3) );
        return d2I2_Cdev_dF2;
    }


    /**
     * Compute the second derivative of the jacobian(determinant of the deformation gradient) wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return second derivative of the jacobian(determinant of the deformation gradient (fourth-order tensor)
     **/
    inline Tensor4_3d compute_d2J_dF2(const Tensor2_3d &F) const
    {
        return F.compute_tensor_cross_product();
    }

    /**
     * Compute the first derivative of the volumetric energy wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return the first derivative of the volumetric energy (second-order tensor)
     **/
    Tensor2_3d compute_dVolumetricEnergy_dF(const Tensor2_3d &F) const;

    /**
     * Compute the second derivative of the volumetric energy wrt deformation gradient F
     * @param Deformation gradient (second-order tensor in 3d)
     * @return the second derivative of the volumetric energy (fourth-order tensor)
     **/
    Tensor4_3d compute_d2VolumetricEnergy_dF2(const Tensor2_3d &F) const;
};
} // end namespace oofem
