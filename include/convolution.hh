// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn
// 
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <array>

#include <general.hh>
#include <grid_fft.hh>

/// @brief base class for convolutions of two or three fields
/// @tparam data_t 
/// @tparam derived_t 
template <typename data_t, typename derived_t>
class BaseConvolver
{
protected:
    std::array<size_t, 3> np_;
    std::array<real_t, 3> length_;
    size_t n_ffts_{0};

public:

    /// @brief Construct a new Base Convolver object
    /// @param N linear grid size
    /// @param L physical box size
    BaseConvolver(const std::array<size_t, 3> &N, const std::array<real_t, 3> &L)
        : np_(N), length_(L), n_ffts_(0) {}


    /// @brief Construct a new Base Convolver object [deleted copy constructor]  
    BaseConvolver( const BaseConvolver& ) = delete;
    
    /// @brief destructor (virtual)
    virtual ~BaseConvolver() {}

    /// @brief return number of convolver FFTs executed so far
    size_t get_fft_count() const noexcept { return n_ffts_; }

    /// @brief reset convolver FFT counter
    void reset_fft_count() noexcept { n_ffts_ = 0; }

    /// @brief implements auto-convolution (squaring) of one Fourier-space field
    /// @tparam kfunc field
    /// @tparam opp output operator
    template <typename kfunc, typename opp>
    void convolve_square(kfunc kf, opp op) {}

    /// @brief implements convolution of two Fourier-space fields
    /// @tparam kfunc1 field 1
    /// @tparam kfunc2 field 2
    /// @tparam opp output operator
    template <typename kfunc1, typename kfunc2, typename opp>
    void convolve2(kfunc1 kf1, kfunc2 kf2, opp op) {}

    /// @brief implements convolution of three Fourier-space fields
    /// @tparam kfunc1 field 1
    /// @tparam kfunc2 field 2
    /// @tparam kfunc3 field 3
    /// @tparam opp output operator
    template <typename kfunc1, typename kfunc2, typename kfunc3, typename opp>
    void convolve3(kfunc1 kf1, kfunc2 kf2, kfunc3 kf3, opp op) {}

    /// @brief implements triple convolution where the first two fields are identical (kf12^2 * kf3)
    /// @tparam kfunc12 repeated inner field
    /// @tparam kfunc3 outer field
    /// @tparam opp output operator
    template <typename kfunc12, typename kfunc3, typename opp>
    void convolve3_square(kfunc12 kf12, kfunc3 kf3, opp op) {}

    /// @brief begin accumulating real-space convolution products
    void begin_accumulate() {}

    /// @brief accumulate weighted square of one field in real space
    template <typename kfunc>
    void accumulate_square(kfunc kf, data_t weight = data_t(1.0)) {}

    /// @brief accumulate weighted product of two fields in real space
    template <typename kfunc1, typename kfunc2>
    void accumulate2(kfunc1 kf1, kfunc2 kf2, data_t weight = data_t(1.0)) {}

    /// @brief accumulate weighted triple convolution in real space
    template <typename kfunc1, typename kfunc2, typename kfunc3>
    void accumulate3(kfunc1 kf1, kfunc2 kf2, kfunc3 kf3, data_t weight = data_t(1.0)) {}

    /// @brief accumulate weighted triple convolution with squared inner field in real space
    template <typename kfunc12, typename kfunc3>
    void accumulate3_square(kfunc12 kf12, kfunc3 kf3, data_t weight = data_t(1.0)) {}

    /// @brief transform accumulated real-space sum forward and write out
    template <typename opp>
    void finish_accumulate(opp op) {}

public:

    /// @brief convolve two gradient fields in Fourier space a_{,i} * b_{,j}
    /// @tparam opp output operator type
    /// @param inl left input field a
    /// @param d1l direction of first gradient (,i)
    /// @param inr right input field b
    /// @param d1r direction of second gradient (,j)
    /// @param output_op output operator
    template <typename opp>
    void convolve_Gradients(Grid_FFT<data_t> &inl, const std::array<int, 1> &d1l,
                            Grid_FFT<data_t> &inr, const std::array<int, 1> &d1r,
                            opp output_op)
    {
        // transform to FS in case fields are not
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        // perform convolution of two gradients
        static_cast<derived_t &>(*this).convolve2(
            // first gradient
            [&inl,&d1l](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inl.gradient(d1l[0],{i,j,k});
                return grad1*inl.kelem(i, j, k);
            },
            // second gradient
            [&inr,&d1r](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inr.gradient(d1r[0],{i,j,k});
                return grad1*inr.kelem(i, j, k);
            },
            // -> output operator
            output_op);
    }

    /// @brief convolve a gradient and a Hessian field in Fourier space a_{,i} * b_{,jk}
    /// @tparam opp output operator type
    /// @param inl left input field a
    /// @param d1l direction of gradient (,i)
    /// @param inr right input field b
    /// @param d2r directions of Hessian (,jk)
    /// @param output_op output operator
    template <typename opp>
    void convolve_Gradient_and_Hessian(Grid_FFT<data_t> &inl, const std::array<int, 1> &d1l,
                                       Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r,
                                       opp output_op)
    {
        // transform to FS in case fields are not
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        // perform convolution of gradient and Hessian
        static_cast<derived_t &>(*this).convolve2(
            // gradient
            [&](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto kk = inl.template get_k<real_t>(i, j, k);
                return ccomplex_t(0.0, -kk[d1l[0]]) * inl.kelem(i, j, k);
            },
            // Hessian
            [&](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto kk = inr.template get_k<real_t>(i, j, k);
                return -kk[d2r[0]] * kk[d2r[1]] * inr.kelem(i, j, k);
            },
            // -> output operator
            output_op);
    }

    /// @brief auto-convolve (square) a Hessian field in Fourier space (a_{,ij})^2
    /// @tparam opp output operator type
    /// @param in input field a
    /// @param d2 directions of Hessian (,ij)
    /// @param output_op output operator
    template <typename opp>
    void convolve_HessianSquare(Grid_FFT<data_t> &in, const std::array<int, 2> &d2, opp output_op)
    {
        // transform to FS in case field is not
        in.FourierTransformForward();
        // perform auto-convolution of Hessian
        static_cast<derived_t &>(*this).convolve_square(
            [&in, &d2](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = in.gradient(d2[0],{i,j,k});
                auto grad2 = in.gradient(d2[1],{i,j,k});
                return grad1*grad2*in.kelem(i, j, k);
            },
            output_op);
    }

    /// @brief auto-convolve (square) the sum of three Hessian fields in Fourier space (a_{,ij} + a_{,kl} + a_{,mn})^2
    /// @tparam opp output operator type
    /// @param in input field a
    /// @param d2a directions of first Hessian (,ij)
    /// @param d2b directions of second Hessian (,kl)
    /// @param d2c directions of third Hessian (,mn)
    /// @param output_op output operator
    template <typename opp>
    void convolve_SumOfHessiansSquare(Grid_FFT<data_t> &in,
                                      const std::array<int, 2> &d2a,
                                      const std::array<int, 2> &d2b,
                                      const std::array<int, 2> &d2c,
                                      opp output_op)
    {
        // transform to FS in case field is not
        in.FourierTransformForward();
        // perform auto-convolution of sum of three Hessians
        static_cast<derived_t &>(*this).convolve_square(
            [&in, &d2a, &d2b, &d2c](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad11 = in.gradient(d2a[0],{i,j,k});
                auto grad12 = in.gradient(d2a[1],{i,j,k});
                auto grad21 = in.gradient(d2b[0],{i,j,k});
                auto grad22 = in.gradient(d2b[1],{i,j,k});
                auto grad31 = in.gradient(d2c[0],{i,j,k});
                auto grad32 = in.gradient(d2c[1],{i,j,k});
                return (grad11*grad12 + grad21*grad22 + grad31*grad32)*in.kelem(i, j, k);
            },
            output_op);
    }

    /// @brief convolve two Hessian fields in Fourier space a_{,ij} * b_{,kl}
    /// @tparam opp output operator type
    /// @param inl left input field a
    /// @param d2l directions of first Hessian (,ij)
    /// @param inr right input field b
    /// @param d2r directions of second Hessian (,kl)
    /// @param output_op output operator
    template <typename opp>
    void convolve_Hessians(Grid_FFT<data_t> &inl, const std::array<int, 2> &d2l,
                           Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r,
                           opp output_op)
    {
        if (&inl == &inr && ((d2l[0] == d2r[0] && d2l[1] == d2r[1]) || (d2l[0] == d2r[1] && d2l[1] == d2r[0])))
        {
            convolve_HessianSquare(inl, d2l, output_op);
            return;
        }
        // transform to FS in case fields are not
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        // perform convolution of Hessians
        static_cast<derived_t &>(*this).convolve2(
            // first Hessian
            [&inl,&d2l](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inl.gradient(d2l[0],{i,j,k});
                auto grad2 = inl.gradient(d2l[1],{i,j,k});
                return grad1*grad2*inl.kelem(i, j, k);
            },
            // second Hessian
            [&inr,&d2r](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inr.gradient(d2r[0],{i,j,k});
                auto grad2 = inr.gradient(d2r[1],{i,j,k});
                return grad1*grad2*inr.kelem(i, j, k);
            },
            // -> output operator
            output_op);
    }

    /// @brief convolve a Fourier-space field with a Hessian field a * b_{,kl}
    /// @tparam opp output operator type
    /// @param inl left input field a
    /// @param inr right input field b
    /// @param d2r directions of Hessian (,kl)
    /// @param output_op output operator
    template <typename opp>
    void convolve_Field_and_Hessian(Grid_FFT<data_t> &inl,
                                    Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r,
                                    opp output_op)
    {
        // transform to FS in case fields are not
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        // perform convolution of field and Hessian
        static_cast<derived_t &>(*this).convolve2(
            // first field
            [&inl](size_t i, size_t j, size_t k) -> ccomplex_t {
                return inl.kelem(i, j, k);
            },
            // second Hessian
            [&inr, &d2r](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inr.gradient(d2r[0],{i,j,k});
                auto grad2 = inr.gradient(d2r[1],{i,j,k});
                return grad1*grad2*inr.kelem(i, j, k);
            },
            // -> output operator
            output_op);
    }

    /// @brief convolve squared Hessian field with another Hessian field in Fourier space (a_{,ij})^2 * c_{,mn}
    /// @tparam opp output operator
    /// @param inl first input field a (squared)
    /// @param d2l directions of squared Hessian (,ij)
    /// @param inr outer input field c
    /// @param d2r directions of outer Hessian (,mn)
    /// @param output_op output operator
    template <typename opp>
    void convolve_HessianSquare_and_Hessian(Grid_FFT<data_t> &inl, const std::array<int, 2> &d2l,
                                            Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r,
                                            opp output_op)
    {
        // transform to FS in case fields are not
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        // perform triple convolution with squared inner Hessian
        static_cast<derived_t &>(*this).convolve3_square(
            // squared Hessian
            [&inl, &d2l](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inl.gradient(d2l[0],{i,j,k});
                auto grad2 = inl.gradient(d2l[1],{i,j,k});
                return grad1*grad2*inl.kelem(i, j, k);
            },
            // outer Hessian
            [&inr, &d2r](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inr.gradient(d2r[0],{i,j,k});
                auto grad2 = inr.gradient(d2r[1],{i,j,k});
                return grad1*grad2*inr.kelem(i, j, k);
            },
            // -> output operator
            output_op);
    }

    /// @brief convolve three Hessian fields in Fourier space a_{,ij} * b_{,kl} * c_{,mn}
    /// @tparam opp output operator
    /// @param inl first input field a
    /// @param d2l directions of first Hessian (,ij)
    /// @param inm second input field b
    /// @param d2m directions of second Hessian (,kl)
    /// @param inr third input field c
    /// @param d2r directions of third Hessian (,mn)
    /// @param output_op output operator
    template <typename opp>
    void convolve_Hessians(Grid_FFT<data_t> &inl, const std::array<int, 2> &d2l,
                           Grid_FFT<data_t> &inm, const std::array<int, 2> &d2m,
                           Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r,
                           opp output_op)
    {
        if (&inl == &inm && ((d2l[0] == d2m[0] && d2l[1] == d2m[1]) || (d2l[0] == d2m[1] && d2l[1] == d2m[0])))
        {
            convolve_HessianSquare_and_Hessian(inl, d2l, inr, d2r, output_op);
            return;
        }
        // transform to FS in case fields are not
        inl.FourierTransformForward();
        inm.FourierTransformForward();
        inr.FourierTransformForward();
        // perform convolution of Hessians
        static_cast<derived_t &>(*this).convolve3(
            // first Hessian
            [&inl, &d2l](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inl.gradient(d2l[0],{i,j,k});
                auto grad2 = inl.gradient(d2l[1],{i,j,k});
                return grad1*grad2*inl.kelem(i, j, k);
            },
            // second Hessian
            [&inm, &d2m](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inm.gradient(d2m[0],{i,j,k});
                auto grad2 = inm.gradient(d2m[1],{i,j,k});
                return grad1*grad2*inm.kelem(i, j, k);
            },
            // third Hessian
            [&inr, &d2r](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inr.gradient(d2r[0],{i,j,k});
                auto grad2 = inr.gradient(d2r[1],{i,j,k});
                return grad1*grad2*inr.kelem(i, j, k);
            },
            // -> output operator
            output_op);
    }

    /// @brief convolve Hessian field with sum of two Hessian fields in Fourier space a_{,ij} * (b_{,kl} + c_{,mn})
    /// @tparam opp output operator type
    /// @param inl left input field a
    /// @param d2l directions of first Hessian (,ij)
    /// @param inr right input field b
    /// @param d2r1 directions of second Hessian (,kl)
    /// @param d2r2 directions of third Hessian (,mn)
    /// @param output_op output operator
    template <typename opp>
    void convolve_SumOfHessians(Grid_FFT<data_t> &inl, const std::array<int, 2> &d2l,
                                Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r1, const std::array<int, 2> &d2r2,
                                opp output_op)
    {
        // transform to FS in case fields are not
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        // perform convolution of Hessians
        static_cast<derived_t &>(*this).convolve2(
            // first Hessian
            [&inl, &d2l](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inl.gradient(d2l[0],{i,j,k});
                auto grad2 = inl.gradient(d2l[1],{i,j,k});
                return grad1*grad2*inl.kelem(i, j, k);
            },
            // second two Hessian and sum
            [&inr, &d2r1, &d2r2](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad11 = inr.gradient(d2r1[0],{i,j,k});
                auto grad12 = inr.gradient(d2r1[1],{i,j,k});
                auto grad21 = inr.gradient(d2r2[0],{i,j,k});
                auto grad22 = inr.gradient(d2r2[1],{i,j,k});
                return (grad11*grad12+grad21*grad22)*inr.kelem(i, j, k);
            },
            // -> output operator
            output_op);
    }

    /// @brief convolve Hessian field with difference of two Hessian fields in Fourier space a_{,ij} * (b_{,kl} - c_{,mn})
    /// @tparam opp output operator type
    /// @param inl left input field a
    /// @param d2l directions of first Hessian (,ij)
    /// @param inr right input field b
    /// @param d2r1 directions of second Hessian (,kl)
    /// @param d2r2 directions of third Hessian (,mn)
    /// @param output_op output operator
    template <typename opp>
    void convolve_DifferenceOfHessians(Grid_FFT<data_t> &inl, const std::array<int, 2> &d2l,
                                       Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r1, const std::array<int, 2> &d2r2,
                                       opp output_op)
    {
        // transform to FS in case fields are not
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        // perform convolution of Hessians
        static_cast<derived_t &>(*this).convolve2(
            // first Hessian
            [&inl, &d2l](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inl.gradient(d2l[0],{i,j,k});
                auto grad2 = inl.gradient(d2l[1],{i,j,k});
                return grad1*grad2*inl.kelem(i, j, k);
            },
            // second two Hessian and difference
            [&inr, &d2r1, &d2r2](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad11 = inr.gradient(d2r1[0],{i,j,k});
                auto grad12 = inr.gradient(d2r1[1],{i,j,k});
                auto grad21 = inr.gradient(d2r2[0],{i,j,k});
                auto grad22 = inr.gradient(d2r2[1],{i,j,k});
                return (grad11*grad12-grad21*grad22)*inr.kelem(i, j, k);
            },
            output_op);
    }
    
    /// @brief performs the multiplication of two fields by convolving them in Fourier space
    /// @tparam opp output operator type
    /// @param inl left input field a)
    /// @param inr right input field b
    /// @param output_op output operator
    template <typename opp> // TOMA
    void multiply_field(Grid_FFT<data_t> &inl, Grid_FFT<data_t> &inr, opp output_op)
    {
        if (&inl == &inr)
        {
            inl.FourierTransformForward();
            static_cast<derived_t &>(*this).convolve_square(
                [&inl](size_t i, size_t j, size_t k) -> ccomplex_t { return inl.kelem(i, j, k); },
                output_op);
            return;
        }
        // transform to FS in case fields are not
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        // perform convolution of Hessians
        static_cast<derived_t &>(*this).convolve2(
            [&inl](size_t i, size_t j, size_t k) -> ccomplex_t {return  inl.kelem(i, j, k); },
            [&inr](size_t i, size_t j, size_t k) -> ccomplex_t {return  inr.kelem(i, j, k); },
            output_op);
    }

    /// @brief accumulate weighted auto-convolution (square) of a Hessian field (a_{,ij})^2
    void accumulate_HessianSquare(Grid_FFT<data_t> &in, const std::array<int, 2> &d2, data_t weight = data_t(1.0))
    {
        in.FourierTransformForward();
        static_cast<derived_t &>(*this).accumulate_square(
            [&in, &d2](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = in.gradient(d2[0],{i,j,k});
                auto grad2 = in.gradient(d2[1],{i,j,k});
                return grad1*grad2*in.kelem(i, j, k);
            },
            weight);
    }

    /// @brief accumulate weighted square of the sum of three Hessian fields (a_{,ij} + a_{,kl} + a_{,mn})^2
    void accumulate_SumOfHessiansSquare(Grid_FFT<data_t> &in,
                                        const std::array<int, 2> &d2a,
                                        const std::array<int, 2> &d2b,
                                        const std::array<int, 2> &d2c,
                                        data_t weight = data_t(1.0))
    {
        in.FourierTransformForward();
        static_cast<derived_t &>(*this).accumulate_square(
            [&in, &d2a, &d2b, &d2c](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad11 = in.gradient(d2a[0],{i,j,k});
                auto grad12 = in.gradient(d2a[1],{i,j,k});
                auto grad21 = in.gradient(d2b[0],{i,j,k});
                auto grad22 = in.gradient(d2b[1],{i,j,k});
                auto grad31 = in.gradient(d2c[0],{i,j,k});
                auto grad32 = in.gradient(d2c[1],{i,j,k});
                return (grad11*grad12 + grad21*grad22 + grad31*grad32)*in.kelem(i, j, k);
            },
            weight);
    }

    /// @brief accumulate weighted product of two Hessian fields a_{,ij} * b_{,kl}
    void accumulate_Hessians(Grid_FFT<data_t> &inl, const std::array<int, 2> &d2l,
                             Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r,
                             data_t weight = data_t(1.0))
    {
        if (&inl == &inr && ((d2l[0] == d2r[0] && d2l[1] == d2r[1]) || (d2l[0] == d2r[1] && d2l[1] == d2r[0])))
        {
            accumulate_HessianSquare(inl, d2l, weight);
            return;
        }
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        static_cast<derived_t &>(*this).accumulate2(
            [&inl,&d2l](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inl.gradient(d2l[0],{i,j,k});
                auto grad2 = inl.gradient(d2l[1],{i,j,k});
                return grad1*grad2*inl.kelem(i, j, k);
            },
            [&inr,&d2r](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inr.gradient(d2r[0],{i,j,k});
                auto grad2 = inr.gradient(d2r[1],{i,j,k});
                return grad1*grad2*inr.kelem(i, j, k);
            },
            weight);
    }

    /// @brief accumulate weighted product of a Fourier-space field with a Hessian field a * b_{,kl}
    void accumulate_Field_and_Hessian(Grid_FFT<data_t> &inl,
                                      Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r,
                                      data_t weight = data_t(1.0))
    {
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        static_cast<derived_t &>(*this).accumulate2(
            [&inl](size_t i, size_t j, size_t k) -> ccomplex_t {
                return inl.kelem(i, j, k);
            },
            [&inr, &d2r](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inr.gradient(d2r[0],{i,j,k});
                auto grad2 = inr.gradient(d2r[1],{i,j,k});
                return grad1*grad2*inr.kelem(i, j, k);
            },
            weight);
    }

    /// @brief accumulate weighted triple convolution (a_{,ij})^2 * c_{,mn}
    void accumulate_HessianSquare_and_Hessian(Grid_FFT<data_t> &inl, const std::array<int, 2> &d2l,
                                              Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r,
                                              data_t weight = data_t(1.0))
    {
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        static_cast<derived_t &>(*this).accumulate3_square(
            [&inl, &d2l](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inl.gradient(d2l[0],{i,j,k});
                auto grad2 = inl.gradient(d2l[1],{i,j,k});
                return grad1*grad2*inl.kelem(i, j, k);
            },
            [&inr, &d2r](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inr.gradient(d2r[0],{i,j,k});
                auto grad2 = inr.gradient(d2r[1],{i,j,k});
                return grad1*grad2*inr.kelem(i, j, k);
            },
            weight);
    }

    /// @brief accumulate weighted triple convolution of three Hessian fields a_{,ij} * b_{,kl} * c_{,mn}
    void accumulate_Hessians(Grid_FFT<data_t> &inl, const std::array<int, 2> &d2l,
                             Grid_FFT<data_t> &inm, const std::array<int, 2> &d2m,
                             Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r,
                             data_t weight = data_t(1.0))
    {
        if (&inl == &inm && ((d2l[0] == d2m[0] && d2l[1] == d2m[1]) || (d2l[0] == d2m[1] && d2l[1] == d2m[0])))
        {
            accumulate_HessianSquare_and_Hessian(inl, d2l, inr, d2r, weight);
            return;
        }
        inl.FourierTransformForward();
        inm.FourierTransformForward();
        inr.FourierTransformForward();
        static_cast<derived_t &>(*this).accumulate3(
            [&inl, &d2l](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inl.gradient(d2l[0],{i,j,k});
                auto grad2 = inl.gradient(d2l[1],{i,j,k});
                return grad1*grad2*inl.kelem(i, j, k);
            },
            [&inm, &d2m](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inm.gradient(d2m[0],{i,j,k});
                auto grad2 = inm.gradient(d2m[1],{i,j,k});
                return grad1*grad2*inm.kelem(i, j, k);
            },
            [&inr, &d2r](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inr.gradient(d2r[0],{i,j,k});
                auto grad2 = inr.gradient(d2r[1],{i,j,k});
                return grad1*grad2*inr.kelem(i, j, k);
            },
            weight);
    }

    /// @brief accumulate weighted product a_{,ij} * (b_{,kl} + b_{,mn})
    void accumulate_SumOfHessians(Grid_FFT<data_t> &inl, const std::array<int, 2> &d2l,
                                  Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r1, const std::array<int, 2> &d2r2,
                                  data_t weight = data_t(1.0))
    {
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        static_cast<derived_t &>(*this).accumulate2(
            [&inl, &d2l](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inl.gradient(d2l[0],{i,j,k});
                auto grad2 = inl.gradient(d2l[1],{i,j,k});
                return grad1*grad2*inl.kelem(i, j, k);
            },
            [&inr, &d2r1, &d2r2](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad11 = inr.gradient(d2r1[0],{i,j,k});
                auto grad12 = inr.gradient(d2r1[1],{i,j,k});
                auto grad21 = inr.gradient(d2r2[0],{i,j,k});
                auto grad22 = inr.gradient(d2r2[1],{i,j,k});
                return (grad11*grad12+grad21*grad22)*inr.kelem(i, j, k);
            },
            weight);
    }

    /// @brief accumulate weighted product a_{,ij} * (b_{,kl} - b_{,mn})
    void accumulate_DifferenceOfHessians(Grid_FFT<data_t> &inl, const std::array<int, 2> &d2l,
                                         Grid_FFT<data_t> &inr, const std::array<int, 2> &d2r1, const std::array<int, 2> &d2r2,
                                         data_t weight = data_t(1.0))
    {
        inl.FourierTransformForward();
        inr.FourierTransformForward();
        static_cast<derived_t &>(*this).accumulate2(
            [&inl, &d2l](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad1 = inl.gradient(d2l[0],{i,j,k});
                auto grad2 = inl.gradient(d2l[1],{i,j,k});
                return grad1*grad2*inl.kelem(i, j, k);
            },
            [&inr, &d2r1, &d2r2](size_t i, size_t j, size_t k) -> ccomplex_t {
                auto grad11 = inr.gradient(d2r1[0],{i,j,k});
                auto grad12 = inr.gradient(d2r1[1],{i,j,k});
                auto grad21 = inr.gradient(d2r2[0],{i,j,k});
                auto grad22 = inr.gradient(d2r2[1],{i,j,k});
                return (grad11*grad12-grad21*grad22)*inr.kelem(i, j, k);
            },
            weight);
    }
};

//! low-level implementation of convolutions -- naive convolution class, ignoring aliasing (no padding)
template <typename data_t>
class NaiveConvolver : public BaseConvolver<data_t, NaiveConvolver<data_t>>
{
protected:
    /// @brief buffer for Fourier transformed fields
    Grid_FFT<data_t> *fbuf1_, *fbuf2_, *facc_;
    bool acc_active_{false}, acc_has_terms_{false};

    /// @brief number of points in each direction
    using BaseConvolver<data_t, NaiveConvolver<data_t>>::np_;

    /// @brief length of each direction
    using BaseConvolver<data_t, NaiveConvolver<data_t>>::length_;

    /// @brief FFT counter
    using BaseConvolver<data_t, NaiveConvolver<data_t>>::n_ffts_;

public:
    /// @brief constructor
    /// @param N number of points in each direction
    /// @param L length of each direction
    NaiveConvolver(const std::array<size_t, 3> &N, const std::array<real_t, 3> &L)
        : BaseConvolver<data_t, NaiveConvolver<data_t>>(N, L)
    {
        fbuf1_ = new Grid_FFT<data_t>(N, length_, true, kspace_id);
        fbuf2_ = nullptr;
        facc_ = nullptr;
    }

    /// @brief destructor
    ~NaiveConvolver()
    {
        delete fbuf1_;
        delete fbuf2_;
        delete facc_;
    }

    /// @brief auto-convolution (squaring) of one field
    template <typename kfunc, typename opp>
    void convolve_square(kfunc kf, opp output_op)
    {
        //... prepare data
        fbuf1_->FourierTransformForward(false);
        this->copy_in(kf, *fbuf1_);

        //... square in real space
        fbuf1_->FourierTransformBackward();
        ++n_ffts_;

#pragma omp parallel for
        for (size_t i = 0; i < fbuf1_->ntot_; ++i)
        {
            (*fbuf1_).relem(i) *= (*fbuf1_).relem(i);
        }
        fbuf1_->FourierTransformForward();
        ++n_ffts_;
//... copy data back
#pragma omp parallel for
        for (size_t i = 0; i < fbuf1_->ntot_; ++i)
        {
            output_op(i, (*fbuf1_)[i]);
        }
    }

    /// @brief convolution of two fields
    template <typename kfunc1, typename kfunc2, typename opp>
    void convolve2(kfunc1 kf1, kfunc2 kf2, opp output_op)
    {
        if (fbuf2_ == nullptr)
        {
            fbuf2_ = new Grid_FFT<data_t>(np_, length_, true, kspace_id);
        }
        //... prepare data 1
        fbuf1_->FourierTransformForward(false);
        this->copy_in(kf1, *fbuf1_);

        //... prepare data 2
        fbuf2_->FourierTransformForward(false);
        this->copy_in(kf2, *fbuf2_);

        //... convolve
        fbuf1_->FourierTransformBackward();
        fbuf2_->FourierTransformBackward();
        n_ffts_ += 2;

#pragma omp parallel for
        for (size_t i = 0; i < fbuf1_->ntot_; ++i)
        {
            (*fbuf2_).relem(i) *= (*fbuf1_).relem(i);
        }
        fbuf2_->FourierTransformForward();
        ++n_ffts_;
        // fbuf2_->dealias();
//... copy data back
#pragma omp parallel for
        for (size_t i = 0; i < fbuf2_->ntot_; ++i)
        {
            output_op(i, (*fbuf2_)[i]);
        }


    }

    /// @brief convolution of three fields
    template <typename kfunc1, typename kfunc2, typename kfunc3, typename opp>
    void convolve3(kfunc1 kf1, kfunc2 kf2, kfunc3 kf3, opp output_op)
    {
        if (fbuf2_ == nullptr)
        {
            fbuf2_ = new Grid_FFT<data_t>(np_, length_, true, kspace_id);
        }
        //... prepare data 1
        fbuf1_->FourierTransformForward(false);
        this->copy_in(kf1, *fbuf1_);

        //... prepare data 2
        fbuf2_->FourierTransformForward(false);
        this->copy_in(kf2, *fbuf2_);

        //... convolve
        fbuf1_->FourierTransformBackward();
        fbuf2_->FourierTransformBackward();
        n_ffts_ += 2;

#pragma omp parallel for
        for (size_t i = 0; i < fbuf1_->ntot_; ++i)
        {
            (*fbuf2_).relem(i) *= (*fbuf1_).relem(i);
        }

        //... prepare data 2
        fbuf1_->FourierTransformForward(false);
        this->copy_in(kf3, *fbuf1_);

        //... convolve
        fbuf1_->FourierTransformBackward();
        ++n_ffts_;

#pragma omp parallel for
        for (size_t i = 0; i < fbuf1_->ntot_; ++i)
        {
            (*fbuf2_).relem(i) *= (*fbuf1_).relem(i);
        }

        fbuf2_->FourierTransformForward();
        ++n_ffts_;
//... copy data back
#pragma omp parallel for
        for (size_t i = 0; i < fbuf2_->ntot_; ++i)
        {
            output_op(i, (*fbuf2_)[i]);
        }
    }

    /// @brief triple convolution where the first two fields are identical
    template <typename kfunc12, typename kfunc3, typename opp>
    void convolve3_square(kfunc12 kf12, kfunc3 kf3, opp output_op)
    {
        if (fbuf2_ == nullptr)
        {
            fbuf2_ = new Grid_FFT<data_t>(np_, length_, true, kspace_id);
        }
        //... prepare and square inner field
        fbuf2_->FourierTransformForward(false);
        this->copy_in(kf12, *fbuf2_);
        fbuf2_->FourierTransformBackward();
        ++n_ffts_;

#pragma omp parallel for
        for (size_t i = 0; i < fbuf2_->ntot_; ++i)
        {
            (*fbuf2_).relem(i) *= (*fbuf2_).relem(i);
        }

        //... prepare outer field
        fbuf1_->FourierTransformForward(false);
        this->copy_in(kf3, *fbuf1_);
        fbuf1_->FourierTransformBackward();
        ++n_ffts_;

#pragma omp parallel for
        for (size_t i = 0; i < fbuf1_->ntot_; ++i)
        {
            (*fbuf2_).relem(i) *= (*fbuf1_).relem(i);
        }

        fbuf2_->FourierTransformForward();
        ++n_ffts_;
//... copy data back
#pragma omp parallel for
        for (size_t i = 0; i < fbuf2_->ntot_; ++i)
        {
            output_op(i, (*fbuf2_)[i]);
        }
    }

    /// @brief begin accumulating real-space products
    void begin_accumulate()
    {
        if (facc_ == nullptr)
        {
            facc_ = new Grid_FFT<data_t>(np_, length_, true, rspace_id);
        }
        facc_->FourierTransformBackward(false);
        acc_active_ = true;
        acc_has_terms_ = false;
    }

    /// @brief accumulate weighted square of one field in real space
    template <typename kfunc>
    void accumulate_square(kfunc kf, data_t weight = data_t(1.0))
    {
        assert(acc_active_);
        fbuf1_->FourierTransformForward(false);
        this->copy_in(kf, *fbuf1_);
        fbuf1_->FourierTransformBackward();
        ++n_ffts_;

        if (!acc_has_terms_)
        {
#pragma omp parallel for
            for (size_t i = 0; i < fbuf1_->ntot_; ++i)
            {
                (*facc_).relem(i) = weight * (*fbuf1_).relem(i) * (*fbuf1_).relem(i);
            }
            acc_has_terms_ = true;
        }
        else
        {
#pragma omp parallel for
            for (size_t i = 0; i < fbuf1_->ntot_; ++i)
            {
                (*facc_).relem(i) += weight * (*fbuf1_).relem(i) * (*fbuf1_).relem(i);
            }
        }
    }

    /// @brief accumulate weighted product of two fields in real space
    template <typename kfunc1, typename kfunc2>
    void accumulate2(kfunc1 kf1, kfunc2 kf2, data_t weight = data_t(1.0))
    {
        assert(acc_active_);
        if (fbuf2_ == nullptr)
        {
            fbuf2_ = new Grid_FFT<data_t>(np_, length_, true, kspace_id);
        }
        fbuf1_->FourierTransformForward(false);
        this->copy_in(kf1, *fbuf1_);
        fbuf2_->FourierTransformForward(false);
        this->copy_in(kf2, *fbuf2_);

        fbuf1_->FourierTransformBackward();
        fbuf2_->FourierTransformBackward();
        n_ffts_ += 2;

        if (!acc_has_terms_)
        {
#pragma omp parallel for
            for (size_t i = 0; i < fbuf1_->ntot_; ++i)
            {
                (*facc_).relem(i) = weight * (*fbuf1_).relem(i) * (*fbuf2_).relem(i);
            }
            acc_has_terms_ = true;
        }
        else
        {
#pragma omp parallel for
            for (size_t i = 0; i < fbuf1_->ntot_; ++i)
            {
                (*facc_).relem(i) += weight * (*fbuf1_).relem(i) * (*fbuf2_).relem(i);
            }
        }
    }

    /// @brief accumulate weighted triple product in real space
    template <typename kfunc1, typename kfunc2, typename kfunc3>
    void accumulate3(kfunc1 kf1, kfunc2 kf2, kfunc3 kf3, data_t weight = data_t(1.0))
    {
        assert(acc_active_);
        if (fbuf2_ == nullptr)
        {
            fbuf2_ = new Grid_FFT<data_t>(np_, length_, true, kspace_id);
        }
        fbuf1_->FourierTransformForward(false);
        this->copy_in(kf1, *fbuf1_);
        fbuf2_->FourierTransformForward(false);
        this->copy_in(kf2, *fbuf2_);

        fbuf1_->FourierTransformBackward();
        fbuf2_->FourierTransformBackward();
        n_ffts_ += 2;

#pragma omp parallel for
        for (size_t i = 0; i < fbuf1_->ntot_; ++i)
        {
            (*fbuf2_).relem(i) *= (*fbuf1_).relem(i);
        }

        fbuf1_->FourierTransformForward(false);
        this->copy_in(kf3, *fbuf1_);
        fbuf1_->FourierTransformBackward();
        ++n_ffts_;

        if (!acc_has_terms_)
        {
#pragma omp parallel for
            for (size_t i = 0; i < fbuf1_->ntot_; ++i)
            {
                (*facc_).relem(i) = weight * (*fbuf2_).relem(i) * (*fbuf1_).relem(i);
            }
            acc_has_terms_ = true;
        }
        else
        {
#pragma omp parallel for
            for (size_t i = 0; i < fbuf1_->ntot_; ++i)
            {
                (*facc_).relem(i) += weight * (*fbuf2_).relem(i) * (*fbuf1_).relem(i);
            }
        }
    }

    /// @brief accumulate weighted triple product with squared inner field in real space
    template <typename kfunc12, typename kfunc3>
    void accumulate3_square(kfunc12 kf12, kfunc3 kf3, data_t weight = data_t(1.0))
    {
        assert(acc_active_);
        if (fbuf2_ == nullptr)
        {
            fbuf2_ = new Grid_FFT<data_t>(np_, length_, true, kspace_id);
        }
        fbuf2_->FourierTransformForward(false);
        this->copy_in(kf12, *fbuf2_);
        fbuf2_->FourierTransformBackward();
        ++n_ffts_;

#pragma omp parallel for
        for (size_t i = 0; i < fbuf2_->ntot_; ++i)
        {
            (*fbuf2_).relem(i) *= (*fbuf2_).relem(i);
        }

        fbuf1_->FourierTransformForward(false);
        this->copy_in(kf3, *fbuf1_);
        fbuf1_->FourierTransformBackward();
        ++n_ffts_;

        if (!acc_has_terms_)
        {
#pragma omp parallel for
            for (size_t i = 0; i < fbuf1_->ntot_; ++i)
            {
                (*facc_).relem(i) = weight * (*fbuf2_).relem(i) * (*fbuf1_).relem(i);
            }
            acc_has_terms_ = true;
        }
        else
        {
#pragma omp parallel for
            for (size_t i = 0; i < fbuf1_->ntot_; ++i)
            {
                (*facc_).relem(i) += weight * (*fbuf2_).relem(i) * (*fbuf1_).relem(i);
            }
        }
    }

    /// @brief transform accumulated real-space sum forward and write out
    template <typename opp>
    void finish_accumulate(opp output_op)
    {
        assert(acc_active_);
        if (!acc_has_terms_)
        {
            facc_->zero();
        }
        facc_->FourierTransformForward();
        ++n_ffts_;
#pragma omp parallel for
        for (size_t i = 0; i < facc_->ntot_; ++i)
        {
            output_op(i, (*facc_)[i]);
        }
        acc_active_ = false;
        acc_has_terms_ = false;
    }

//--------------------------------------------------------------------------------------------------------
private:

    /// @brief copy data into a grid
    /// @tparam kfunc abstract function type generating data
    /// @param kf abstract function generating data
    /// @param g grid to copy data into
    template <typename kfunc>
    void copy_in(kfunc kf, Grid_FFT<data_t> &g)
    {
#pragma omp parallel for
        for (size_t i = 0; i < g.size(0); ++i)
        {
            for (size_t j = 0; j < g.size(1); ++j)
            {
                for (size_t k = 0; k < g.size(2); ++k)
                {
                    g.kelem(i, j, k) = kf(i, j, k);
                }
            }
        }
    }
};

//! convolution class, respecting Orszag's 3/2 rule (padding in Fourier space to avoid aliasing)
template <typename data_t>
class OrszagConvolver : public BaseConvolver<data_t, OrszagConvolver<data_t>>
{
private:
    /// @brief buffer for Fourier transformed fields
    Grid_FFT<data_t> *f1p_, *f2p_, *faccp_, *fbuf_;
    bool acc_active_{false}, acc_has_terms_{false};

    using BaseConvolver<data_t, OrszagConvolver<data_t>>::np_;
    using BaseConvolver<data_t, OrszagConvolver<data_t>>::length_;
    using BaseConvolver<data_t, OrszagConvolver<data_t>>::n_ffts_;

    ccomplex_t *crecvbuf_; //!< receive buffer for MPI (complex)
    real_t *recvbuf_;     //!< receive buffer for MPI (real)
    size_t maxslicesz_;  //!< maximum size of a slice
    std::vector<ptrdiff_t> offsets_, offsetsp_; //!< offsets for MPI
    std::vector<size_t> sizes_, sizesp_;       //!< sizes for MPI

    /// @brief get task index for a given index
    /// @param index index
    /// @param offsets offsets
    /// @param sizes sizes
    /// @param ntasks number of tasks
    int get_task(ptrdiff_t index, const std::vector<ptrdiff_t> &offsets, const std::vector<size_t> &sizes, const int ntasks)
    {
        int itask = 0;
        while (itask < ntasks - 1 && offsets[itask + 1] <= index)
            ++itask;
        return itask;
    }

public:

    /// @brief constructor
    /// @param N grid size
    /// @param L grid length
    OrszagConvolver(const std::array<size_t, 3> &N, const std::array<real_t, 3> &L)
        : BaseConvolver<data_t, OrszagConvolver<data_t>>({3 * N[0] / 2, 3 * N[1] / 2, 3 * N[2] / 2}, L)
    {
        //... create temporaries
        f1p_ = new Grid_FFT<data_t>(np_, length_, true, kspace_id);
        f2p_ = nullptr; // allocated lazily only if a two-field product is requested
        faccp_ = nullptr; // allocated lazily only if real-space accumulation is requested
        fbuf_ = new Grid_FFT<data_t>(N, length_, true, kspace_id); // needed for MPI, or for triple conv.

#if defined(USE_MPI)
        maxslicesz_ = f1p_->sizes_[1] * f1p_->sizes_[3] * 2;

        crecvbuf_ = new ccomplex_t[maxslicesz_ / 2];
        recvbuf_ = reinterpret_cast<real_t *>(&crecvbuf_[0]);

        int ntasks(MPI::get_size());

        offsets_.assign(ntasks, 0);
        offsetsp_.assign(ntasks, 0);
        sizes_.assign(ntasks, 0);
        sizesp_.assign(ntasks, 0);

        size_t tsize = N[0], tsizep = f1p_->size(0);

        MPI_Allgather(&fbuf_->local_1_start_, 1, MPI_LONG_LONG, &offsets_[0], 1,
                      MPI_LONG_LONG, MPI_COMM_WORLD);
        MPI_Allgather(&f1p_->local_1_start_, 1, MPI_LONG_LONG, &offsetsp_[0], 1,
                      MPI_LONG_LONG, MPI_COMM_WORLD);
        MPI_Allgather(&tsize, 1, MPI_LONG_LONG, &sizes_[0], 1, MPI_LONG_LONG,
                      MPI_COMM_WORLD);
        MPI_Allgather(&tsizep, 1, MPI_LONG_LONG, &sizesp_[0], 1, MPI_LONG_LONG,
                      MPI_COMM_WORLD);
#endif
    }

    /// @brief destructor
    ~OrszagConvolver()
    {
        delete f1p_;
        delete f2p_;
        delete faccp_;
        delete fbuf_;
#if defined(USE_MPI)
        delete[] crecvbuf_;
#endif
    }

    /// @brief auto-convolve (square) one field using a single padded buffer
    /// @tparam kfunc abstract function type generating data for the field
    /// @tparam opp abstract function type for the output operation
    template <typename kfunc, typename opp>
    void convolve_square(kfunc kf, opp output_op)
    {
        //... prepare data
        f1p_->FourierTransformForward(false);
        this->pad_insert(kf, *f1p_);

        //... square in real space
        f1p_->FourierTransformBackward();
        ++n_ffts_;

#pragma omp parallel for
        for (size_t i = 0; i < f1p_->ntot_; ++i)
        {
            (*f1p_).relem(i) *= (*f1p_).relem(i);
        }
        f1p_->FourierTransformForward();
        ++n_ffts_;
        //... copy data back
        unpad(*f1p_, output_op);
    }

    /// @brief convolve two fields
    /// @tparam kfunc1 abstract function type generating data for the first field
    /// @tparam kfunc2 abstract function type generating data for the second field
    /// @tparam opp abstract function type for the output operation
    template <typename kfunc1, typename kfunc2, typename opp>
    void convolve2(kfunc1 kf1, kfunc2 kf2, opp output_op)
    {
        if (f2p_ == nullptr)
        {
            f2p_ = new Grid_FFT<data_t>(np_, length_, true, kspace_id);
        }

        //... prepare data 1
        f1p_->FourierTransformForward(false);
        this->pad_insert(kf1, *f1p_);

        //... prepare data 2
        f2p_->FourierTransformForward(false);
        this->pad_insert(kf2, *f2p_);

        //... convolve
        f1p_->FourierTransformBackward();
        f2p_->FourierTransformBackward();
        n_ffts_ += 2;

#pragma omp parallel for
        for (size_t i = 0; i < f1p_->ntot_; ++i)
        {
            (*f2p_).relem(i) *= (*f1p_).relem(i);
        }
        f2p_->FourierTransformForward();
        ++n_ffts_;
        //... copy data back
        unpad(*f2p_, output_op);
    }

    /// @brief convolve three fields
    /// @tparam kfunc1 abstract function type generating data for the first field
    /// @tparam kfunc2 abstract function type generating data for the second field
    /// @tparam kfunc3 abstract function type generating data for the third field
    /// @tparam opp abstract function type for the output operation
    template <typename kfunc1, typename kfunc2, typename kfunc3, typename opp>
    void convolve3(kfunc1 kf1, kfunc2 kf2, kfunc3 kf3, opp output_op)
    {
        auto assign_to = [](auto &g) { return [&](auto i, auto v) { g[i] = v; }; };
        fbuf_->FourierTransformForward(false);
        convolve2(kf1, kf2, assign_to(*fbuf_));
        convolve2([&](size_t i, size_t j, size_t k) -> ccomplex_t { return fbuf_->kelem(i, j, k); }, kf3, output_op);
    }

    /// @brief convolve squared field with outer field (kf12^2 * kf3)
    /// @tparam kfunc12 abstract function type generating data for the squared inner field
    /// @tparam kfunc3 abstract function type generating data for the outer field
    /// @tparam opp abstract function type for the output operation
    template <typename kfunc12, typename kfunc3, typename opp>
    void convolve3_square(kfunc12 kf12, kfunc3 kf3, opp output_op)
    {
        auto assign_to = [](auto &g) { return [&](auto i, auto v) { g[i] = v; }; };
        fbuf_->FourierTransformForward(false);
        convolve_square(kf12, assign_to(*fbuf_));
        convolve2([&](size_t i, size_t j, size_t k) -> ccomplex_t { return fbuf_->kelem(i, j, k); }, kf3, output_op);
    }

    /// @brief begin accumulating real-space products on the padded grid
    void begin_accumulate()
    {
        if (faccp_ == nullptr)
        {
            faccp_ = new Grid_FFT<data_t>(np_, length_, true, rspace_id);
        }
        faccp_->FourierTransformBackward(false);
        acc_active_ = true;
        acc_has_terms_ = false;
    }

    /// @brief accumulate weighted square of one field on the padded real-space grid
    template <typename kfunc>
    void accumulate_square(kfunc kf, data_t weight = data_t(1.0))
    {
        assert(acc_active_);
        f1p_->FourierTransformForward(false);
        this->pad_insert(kf, *f1p_);
        f1p_->FourierTransformBackward();
        ++n_ffts_;

        if (!acc_has_terms_)
        {
#pragma omp parallel for
            for (size_t i = 0; i < f1p_->ntot_; ++i)
            {
                (*faccp_).relem(i) = weight * (*f1p_).relem(i) * (*f1p_).relem(i);
            }
            acc_has_terms_ = true;
        }
        else
        {
#pragma omp parallel for
            for (size_t i = 0; i < f1p_->ntot_; ++i)
            {
                (*faccp_).relem(i) += weight * (*f1p_).relem(i) * (*f1p_).relem(i);
            }
        }
    }

    /// @brief accumulate weighted product of two fields on the padded real-space grid
    template <typename kfunc1, typename kfunc2>
    void accumulate2(kfunc1 kf1, kfunc2 kf2, data_t weight = data_t(1.0))
    {
        assert(acc_active_);
        if (f2p_ == nullptr)
        {
            f2p_ = new Grid_FFT<data_t>(np_, length_, true, kspace_id);
        }

        f1p_->FourierTransformForward(false);
        this->pad_insert(kf1, *f1p_);

        f2p_->FourierTransformForward(false);
        this->pad_insert(kf2, *f2p_);

        f1p_->FourierTransformBackward();
        f2p_->FourierTransformBackward();
        n_ffts_ += 2;

        if (!acc_has_terms_)
        {
#pragma omp parallel for
            for (size_t i = 0; i < f1p_->ntot_; ++i)
            {
                (*faccp_).relem(i) = weight * (*f1p_).relem(i) * (*f2p_).relem(i);
            }
            acc_has_terms_ = true;
        }
        else
        {
#pragma omp parallel for
            for (size_t i = 0; i < f1p_->ntot_; ++i)
            {
                (*faccp_).relem(i) += weight * (*f1p_).relem(i) * (*f2p_).relem(i);
            }
        }
    }

    /// @brief accumulate weighted dealiased triple product on the padded real-space grid
    template <typename kfunc1, typename kfunc2, typename kfunc3>
    void accumulate3(kfunc1 kf1, kfunc2 kf2, kfunc3 kf3, data_t weight = data_t(1.0))
    {
        assert(acc_active_);
        auto assign_to = [](auto &g) { return [&](auto i, auto v) { g[i] = v; }; };
        fbuf_->FourierTransformForward(false);
        convolve2(kf1, kf2, assign_to(*fbuf_));
        accumulate2([&](size_t i, size_t j, size_t k) -> ccomplex_t { return fbuf_->kelem(i, j, k); }, kf3, weight);
    }

    /// @brief accumulate weighted dealiased triple product with squared inner field on the padded real-space grid
    template <typename kfunc12, typename kfunc3>
    void accumulate3_square(kfunc12 kf12, kfunc3 kf3, data_t weight = data_t(1.0))
    {
        assert(acc_active_);
        auto assign_to = [](auto &g) { return [&](auto i, auto v) { g[i] = v; }; };
        fbuf_->FourierTransformForward(false);
        convolve_square(kf12, assign_to(*fbuf_));
        accumulate2([&](size_t i, size_t j, size_t k) -> ccomplex_t { return fbuf_->kelem(i, j, k); }, kf3, weight);
    }

    /// @brief transform accumulated padded real-space sum forward and unpad to output
    template <typename opp>
    void finish_accumulate(opp output_op)
    {
        assert(acc_active_);
        if (!acc_has_terms_)
        {
            faccp_->zero();
        }
        faccp_->FourierTransformForward();
        ++n_ffts_;
        unpad(*faccp_, output_op);
        acc_active_ = false;
        acc_has_terms_ = false;
    }

    // template< typename opp >
    // void test_pad_unpad( Grid_FFT<data_t> & in, Grid_FFT<data_t> & res, opp op )
    // {
    //     //... prepare data 1
    //     f1p_->FourierTransformForward(false);
    //     this->pad_insert( [&in]( size_t i, size_t j, size_t k ){return in.kelem(i,j,k);}, *f1p_ );
    //     f1p_->FourierTransformBackward();
    //     f1p_->FourierTransformForward();
    //     res.FourierTransformForward();
    //     unpad(*f1p_, res, op);
    // }

private:

    /// @brief unpad the result of a convolution and copy it to a grid
    /// @tparam kdep_functor abstract function type generating data for the result
    /// @param kfunc abstract function generating data for the result
    /// @param fp grid to copy the result to
    template <typename kdep_functor>
    void pad_insert( kdep_functor kfunc, Grid_FFT<data_t> &fp)
    {
        const real_t rfac = std::pow(1.5, 1.5);

#if !defined(USE_MPI)
        const size_t nhalf[3] = {fp.n_[0] / 3, fp.n_[1] / 3, fp.n_[2] / 3};

        fp.zero();

        #pragma omp parallel for
        for (size_t i = 0; i < 2 * fp.size(0) / 3; ++i)
        {
            size_t ip = (i > nhalf[0]) ? i + nhalf[0] : i;
            for (size_t j = 0; j < 2 * fp.size(1) / 3; ++j)
            {
                size_t jp = (j > nhalf[1]) ? j + nhalf[1] : j;
                for (size_t k = 0; k < nhalf[2]+1; ++k)
                {
                    size_t kp = (k > nhalf[2]) ? k + nhalf[2] : k;
                    fp.kelem(ip, jp, kp) = kfunc(i, j, k) * rfac;
                }
            }
        }
#else
        fbuf_->FourierTransformForward(false);
        
        #pragma omp parallel for
        for (size_t i = 0; i < fbuf_->size(0); ++i)
        {
            for (size_t j = 0; j < fbuf_->size(1); ++j)
            {
                for (size_t k = 0; k < fbuf_->size(2); ++k)
                {
                    fbuf_->kelem(i, j, k) = kfunc(i, j, k) * rfac;
                }
            }
        }

        fbuf_->FourierInterpolateCopyTo( fp );
        
#endif //defined(USE_MPI)
    }

    /// @brief unpad the result of a convolution and write it to an output operator
    /// @tparam operator_t abstract function type for the output operation
    /// @param fp grid to copy the result from
    /// @param output_op abstract function to write the result to
    template <typename operator_t>
    void unpad( Grid_FFT<data_t> &fp, operator_t output_op)
    {
        const real_t rfac = std::sqrt(fp.n_[0] * fp.n_[1] * fp.n_[2]) / std::sqrt(fbuf_->n_[0] * fbuf_->n_[1] * fbuf_->n_[2]);

        // make sure we're in Fourier space...
        assert(fp.space_ == kspace_id);

#if !defined(USE_MPI) ////////////////////////////////////////////////////////////////////////////////////
        fbuf_->FourierTransformForward(false);
        size_t nhalf[3] = {fbuf_->n_[0] / 2, fbuf_->n_[1] / 2, fbuf_->n_[2] / 2};

        #pragma omp parallel for
        for (size_t i = 0; i < fbuf_->size(0); ++i)
        {
            size_t ip = (i > nhalf[0]) ? i + nhalf[0] : i;
            for (size_t j = 0; j < fbuf_->size(1); ++j)
            {
                size_t jp = (j > nhalf[1]) ? j + nhalf[1] : j;
                for (size_t k = 0; k < fbuf_->size(2); ++k)
                {
                    size_t kp = (k > nhalf[2]) ? k + nhalf[2] : k;
                    fbuf_->kelem(i, j, k) = fp.kelem(ip, jp, kp) / rfac;
                    // zero Nyquist modes since they are not unique after convolution
                    if( i==nhalf[0]||j==nhalf[1]||k==nhalf[2]){
                        fbuf_->kelem(i, j, k) = 0.0; 
                    }
                }
            }
        }

        //... copy data back
        #pragma omp parallel for
        for (size_t i = 0; i < fbuf_->ntot_; ++i)
        {
            output_op(i, (*fbuf_)[i]);
        }

#else /// then USE_MPI is defined //////////////////////////////////////////////////////////////
    
        fp.FourierInterpolateCopyTo( *fbuf_ );

        //... copy data back
        #pragma omp parallel for
        for (size_t i = 0; i < fbuf_->ntot_; ++i)
        {

            output_op(i, (*fbuf_)[i] / rfac);
        }

#endif //defined(USE_MPI)
    }
};
