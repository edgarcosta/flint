/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb.h"
#include "acb_ppav.h"

TEST_FUNCTION_START(acb_ppav_agm, state)
{
    arf_t eps;
    slong prec = 256;

    arf_init(eps);
    arf_set_d(eps, 1e-10);

    /* g = 1: the classical AGM(1, sqrt 2). Reference midpoint verified to 85
       digits by the review panel; here we pin the first 49 digits. */
    {
        acb_ptr a = _acb_vec_init(2), r = _acb_vec_init(1);
        arb_t want;

        acb_one(&a[0]);
        acb_set_ui(&a[1], 2);
        acb_sqrt(&a[1], &a[1], 400);

        if (acb_ppav_agm(r, a, eps, 1, prec) != 1)
            TEST_FUNCTION_FAIL("g1 code\n");

        arb_init(want);
        arb_set_str(want,
            "1.1981402347355922074399224922803238782272126632156 +/- 1e-48",
            prec);
        if (!arb_contains(want, acb_realref(&r[0]))
            || !arb_contains_zero(acb_imagref(&r[0])))
            TEST_FUNCTION_FAIL("g1 value %{acb}\n", &r[0]);

        /* accuracy must scale with prec: a RELATIVE bound, not a fat ball. */
        if (acb_rel_accuracy_bits(&r[0]) < prec - 30)
            TEST_FUNCTION_FAIL("g1 accuracy %wd < %wd: %{acb}\n",
                acb_rel_accuracy_bits(&r[0]), prec - 30, &r[0]);

        arb_clear(want);
        _acb_vec_clear(a, 2);
        _acb_vec_clear(r, 1);
    }

    /* g = 2: a = (1,1,1,1) is a fixed point of the Borchardt step, so the mean
       is exactly 1. This independently guards the 2^g normalization: dropping
       the /4 sends the iteration to infinity. */
    {
        acb_ptr a2 = _acb_vec_init(4), r2 = _acb_vec_init(1);
        acb_t one;
        slong j;

        for (j = 0; j < 4; j++)
            acb_one(&a2[j]);

        if (acb_ppav_agm(r2, a2, eps, 2, prec) != 1)
            TEST_FUNCTION_FAIL("g2 code\n");

        acb_init(one);
        acb_one(one);
        if (!acb_contains(&r2[0], one))
            TEST_FUNCTION_FAIL("g2 value %{acb}\n", &r2[0]);

        /* the all-ones input is exact, so the enclosure must be tight. */
        if (acb_rel_accuracy_bits(&r2[0]) < prec - 30)
            TEST_FUNCTION_FAIL("g2 accuracy %wd < %wd: %{acb}\n",
                acb_rel_accuracy_bits(&r2[0]), prec - 30, &r2[0]);

        acb_clear(one);
        _acb_vec_clear(a2, 4);
        _acb_vec_clear(r2, 1);
    }

    arf_clear(eps);

    TEST_FUNCTION_END(state);
}
