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

#define TEST_PRINT_PRECISION_1 6
#define TEST_PRINT_PRECISION_2 18

#include "gcem_tests.hpp"

int main()
{
    print_begin("fabsf");

    //

    GCEM_TEST_COMPARE_VALS(gcem::fabsf,std::fabsf, 0.0f);
    GCEM_TEST_COMPARE_VALS(gcem::fabsf,std::fabsf,-0.0f);
    GCEM_TEST_COMPARE_VALS(gcem::fabsf,std::fabsf, 1.0f);
    GCEM_TEST_COMPARE_VALS(gcem::fabsf,std::fabsf,-1.0f);

    GCEM_TEST_COMPARE_VALS(gcem::fabsf,std::fabsf, std::numeric_limits<float>::lowest());
    GCEM_TEST_COMPARE_VALS(gcem::fabsf,std::fabsf, std::numeric_limits<float>::min());
    GCEM_TEST_COMPARE_VALS(gcem::fabsf,std::fabsf, std::numeric_limits<float>::max());

    GCEM_TEST_COMPARE_VALS(gcem::fabsf,std::fabsf, std::numeric_limits<float>::quiet_NaN());
    GCEM_TEST_COMPARE_VALS(gcem::fabsf,std::fabsf, std::numeric_limits<float>::infinity());
    GCEM_TEST_COMPARE_VALS(gcem::fabsf,std::fabsf, -std::numeric_limits<float>::infinity());

    //

    print_final("fabsf");

    return 0;
}
