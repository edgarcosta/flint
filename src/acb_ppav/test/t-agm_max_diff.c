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

TEST_FUNCTION_START(acb_ppav_agm_max_diff, state)
{
    arf_t delta;
    acb_ptr a, c;

    arf_init(delta);

    /* (A) a = (2, 3): ratio 1/2. Expect 0.5 <= delta <= 0.5*(1+2^-40), finite. */
    a = _acb_vec_init(2);
    acb_set_ui(&a[0], 2);
    acb_set_ui(&a[1], 3);
    acb_ppav_agm_max_diff(delta, a, 1, 128);
    if (!arf_is_finite(delta) || arf_cmp_d(delta, 0.5) < 0
        || arf_cmp_d(delta, 0.5000001) > 0)
        TEST_FUNCTION_FAIL("A: delta = %{arf}\n", delta);
    _acb_vec_clear(a, 2);

    /* (B) a0 contains zero -> delta must be +inf (not garbage). */
    c = _acb_vec_init(2);
    acb_zero(&c[0]);
    arb_add_error_2exp_si(acb_realref(&c[0]), 0);
    acb_set_ui(&c[1], 1);
    acb_ppav_agm_max_diff(delta, c, 1, 64);
    if (!arf_is_pos_inf(delta))
        TEST_FUNCTION_FAIL("B: expected +inf, got %{arf}\n", delta);
    _acb_vec_clear(c, 2);

    arf_clear(delta);

    TEST_FUNCTION_END(state);
}
