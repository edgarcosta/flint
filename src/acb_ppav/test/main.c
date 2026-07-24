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

#include "t-agm.c"
#include "t-agm_half_plane.c"
#include "t-agm_max_diff.c"
#include "t-periods_from_theta2.c"
#include "t-periods_overlap_sp4.c"
#include "t-rosenhain.c"
#include "t-siegel_coset_nb.c"
#include "t-theta2.c"
#include "t-theta4.c"
#include "t-theta_naive_ref.c"
#include "t-thomae_chain.c"
#include "t-weierstrass.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_ppav_agm),
    TEST_FUNCTION(acb_ppav_agm_half_plane),
    TEST_FUNCTION(acb_ppav_agm_max_diff),
    TEST_FUNCTION(acb_ppav_periods_from_theta2),
    TEST_FUNCTION(acb_ppav_periods_from_theta2_regression),
    TEST_FUNCTION(acb_ppav_periods_from_theta2_scalar),
    TEST_FUNCTION(acb_ppav_periods_from_theta2_invalid),
    TEST_FUNCTION(acb_ppav_periods_from_theta2_boundary),
    TEST_FUNCTION(acb_ppav_periods_from_theta2_escalation),
    TEST_FUNCTION(acb_ppav_periods_overlap_sp4),
    TEST_FUNCTION(acb_ppav_g2_rosenhain),
    TEST_FUNCTION(acb_ppav_g2_siegel_coset_nb),
    TEST_FUNCTION(acb_ppav_g2_theta2),
    TEST_FUNCTION(acb_ppav_g2_theta4),
    TEST_FUNCTION(acb_ppav_theta_naive_ref),
    TEST_FUNCTION(acb_ppav_thomae_chain),
    TEST_FUNCTION(acb_ppav_weierstrass),
};

/* main function *************************************************************/

TEST_MAIN(tests)

