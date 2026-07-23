/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_ppav.h"

TEST_FUNCTION_START(acb_ppav_agm_half_plane, state)
{
    arf_t eps, beta;
    acb_ptr a;
    slong j;

    arf_init(eps);
    arf_init(beta);
    arf_set_d(eps, 1e-10);
    a = _acb_vec_init(4);

    /* (A) all positive reals: a common half-plane certainly exists. */
    for (j = 0; j < 4; j++)
        acb_set_ui(&a[j], j + 1);
    if (acb_ppav_agm_half_plane(beta, a, eps, 2, 128) != 1)
        TEST_FUNCTION_FAIL("A: expected 1\n");

    /* (B) 1, -1, i, -i point in four balanced directions: no common
       half-plane can exist. */
    acb_set_ui(&a[0], 1);
    acb_set_si(&a[1], -1);
    acb_onei(&a[2]);
    acb_neg(&a[3], &a[2]);
    if (acb_ppav_agm_half_plane(beta, a, eps, 2, 128) != 0)
        TEST_FUNCTION_FAIL("B: expected 0\n");

    /* (C) one entry is a fat ball straddling the origin at tiny precision:
       undecidable. */
    for (j = 0; j < 4; j++)
        acb_set_ui(&a[j], 1);
    arb_add_error_2exp_si(acb_realref(&a[3]), 0);
    arb_add_error_2exp_si(acb_imagref(&a[3]), 0);
    if (acb_ppav_agm_half_plane(beta, a, eps, 2, 32) != 2)
        TEST_FUNCTION_FAIL("C: expected 2\n");

    _acb_vec_clear(a, 4);
    arf_clear(eps);
    arf_clear(beta);

    TEST_FUNCTION_END(state);
}
