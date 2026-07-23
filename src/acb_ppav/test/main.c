/*
    Copyright (C) 2026 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "test_helpers.h"
#include "acb_ppav.h"

#include "t-agm_half_plane.c"
#include "t-agm_max_diff.c"
#include "t-siegel_coset_nb.c"
#include "t-theta_naive_ref.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_ppav_agm_half_plane),
    TEST_FUNCTION(acb_ppav_agm_max_diff),
    TEST_FUNCTION(acb_ppav_g2_siegel_coset_nb),
    TEST_FUNCTION(acb_ppav_theta_naive_ref),
};

/* main function *************************************************************/

TEST_MAIN(tests)

