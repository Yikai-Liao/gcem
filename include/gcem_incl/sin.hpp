/*################################################################################
  ##
  ##   Copyright (C) 2016-2024 Keith O'Hara
  ##
  ##   This file is part of the GCE-Math C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

/*
 * compile-time sine function using tan(x/2)
 * 
 * see eq. 5.4.8 in Numerical Recipes
 */

#ifndef _gcem_sin_HPP
#define _gcem_sin_HPP

namespace internal
{

template<typename T>
constexpr
T
sin_compute(const T x)
noexcept
{
    return T(2)*x/(T(1) + x*x);
}

template<typename T>
constexpr
T
sin_compute_chebyshev(const T x)
noexcept {
    constexpr T two_pie = 2 * T(GCEM_PI);
    // limit x to range(-pi, pi)
    const T x1 = x - two_pie * floor(x / two_pie) - T(GCEM_PI);
    const T x2 = x1 * x1;
    // Chebyshev approximation of sin(x) with degree 29
    // weights are from: https://publik-void.github.io/sin-cos-approximations/#_sin_rel_error_minimized_degree_29
    return -x1*(1. + x2*(-0.166666666666666666666666666666666629 + x2*(0.00833333333333333333333333333333222318 + x2*(-0.000198412698412698412698412698399723487 + x2*(2.75573192239858906525573184281026561e-6 + x2*(-2.50521083854417187750518138734779333e-8 + x2*(1.60590438368216145993211483359730933e-10 + x2*(-7.64716373181981646407639862566185217e-13 + x2*(2.81145725434551937513944561966331531e-15 + x2*(-8.22063524662315930664562316486611394e-18 + x2*(1.9572941062679701610307767721364446e-20 + x2*(-3.8681701397118386284141935821892017e-23 + x2*(6.44694091890808795934274904084719604e-26 + x2*(-9.18181103282828321298337267365314117e-29 + 1.10854165613066936680763885008011392e-31*x2))))))))))))));
}

template<typename T>
constexpr
T
sin_check(const T x)
noexcept
{
    return( // NaN check
            is_nan(x) ? \
                GCLIM<T>::quiet_NaN() :
            // indistinguishable from zero
            GCLIM<T>::min() > abs(x) ? \
                T(0) :
                sin_compute_chebyshev(x));
}

}

/**
 * Compile-time sine function
 *
 * @param x a real-valued input.
 * @return the sine function using \f[ \sin(x) = \frac{2\tan(x/2)}{1+\tan^2(x/2)} \f]
 */

template<typename T>
constexpr
return_t<T>
sin(const T x)
noexcept
{
    return internal::sin_check( static_cast<return_t<T>>(x) );
}

#endif
