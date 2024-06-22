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
 * compile-time cosine function using 30th degree Chebyshev approximation
 */

#ifndef _gcem_cos_HPP
#define _gcem_cos_HPP

namespace internal
{

template<typename T>
constexpr
T
cos_compute(const T x)
noexcept
{
    return( T(1) - x*x)/(T(1) + x*x );
}

template<typename T>
constexpr
T
cos_compute_chebyshev(const T x)
noexcept {
    // limit x to range(0, 2*pi)
    const T x1 = x - T(GCEM_2PI) * floor(x / T(GCEM_2PI));
    const T x2 = x1 * x1;
    // Chebyshev approximation of sin(x) with degree 30
    // weights are from: https://publik-void.github.io/sin-cos-approximations/#_cos_rel_error_minimized_degree_30
    return 1. + x2*(-0.499999999999999999999999999999999997 + x2*(0.0416666666666666666666666666666665864 + x2*(-0.00138888888888888888888888888888792216 + x2*(0.0000248015873015873015873015872954247692 + x2*(-2.75573192239858906525573168323374212e-7 + x2*(2.08767569878680989792094774070130099e-9 + x2*(-1.147074559772972471374259582051505e-11 + x2*(4.77947733238738528351115798998742201e-14 + x2*(-1.56192069685862134697357839719554288e-16 + x2*(4.11031762331127151219297331433042247e-19 + x2*(-8.8967913919984472368244650875181655e-22 + x2*(1.61173755446230144849643763129805917e-24 + x2*(-2.47959193649582299714152529144837331e-27 + x2*(3.27913485904815466449434894827885198e-30 - 3.69081529290172836989264646507263976e-33*x2))))))))))))));
}

template<typename T>
constexpr
T
cos_check(const T x)
noexcept
{
    return( // NaN check
            is_nan(x) ? \
                GCLIM<T>::quiet_NaN() :
            // indistinguishable from 0
            GCLIM<T>::min() > abs(x) ? 
                T(1) :
                cos_compute_chebyshev(x));
}

}

/**
 * Compile-time cosine function
 *
 * @param x a real-valued input.
 * @return the cosine function using \f[ \cos(x) = \frac{1-\tan^2(x/2)}{1+\tan^2(x/2)} \f]
 */

template<typename T>
constexpr
return_t<T>
cos(const T x)
noexcept
{
    return internal::cos_check( static_cast<return_t<T>>(x) );
}

#endif
