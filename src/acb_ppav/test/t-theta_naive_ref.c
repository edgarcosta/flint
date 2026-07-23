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
#include "acb_mat.h"
#include "acb_theta.h"
#include "acb_ppav.h"
#include "theta_naive_ref.h"

/* Pin the acb_theta characteristic ordering against the convention-explicit
   naive oracle: for a random reduced tau, every one of the 16 entries of
   acb_theta_all(..., sqr = 1) must enclose the same theta_{a,b}(0,tau)^2 the
   oracle computes from the definition. The index -> (a,b) map is the acb_theta
   packing (acb_theta.rst, char_set_slong_vec.c): the 4-bit index is (a,b) with
   most significant bit first, i.e. index = 8 a0 + 4 a1 + 2 b0 + b1. */

TEST_FUNCTION_START(acb_ppav_theta_naive_ref, state)
{
    slong prec = 200, k;
    acb_mat_t tau;
    acb_ptr th, z;

    acb_mat_init(tau, 2, 2);
    acb_siegel_randtest_reduced(tau, state, prec, 4);      /* acb_theta.h */

    th = _acb_vec_init(16);
    z = _acb_vec_init(2);                                  /* z = 0 */
    acb_theta_all(th, z, tau, 1, prec);                    /* squares, acb_theta.h */

    for (k = 0; k < 16; k++)
    {
        int a[2] = { (k >> 3) & 1, (k >> 2) & 1 };
        int b[2] = { (k >> 1) & 1, k & 1 };
        acb_t ref;

        acb_init(ref);
        _acb_ppav_naive_theta2_ab(ref, a, b, tau, prec);

        if (!acb_overlaps(th + k, ref))
        {
            flint_printf("FAIL: char %wd disagrees with naive theta^2\n", k);
            flint_printf("a = (%d,%d), b = (%d,%d)\n", a[0], a[1], b[0], b[1]);
            flint_printf("acb_theta_all: "); acb_printd(th + k, 15); flint_printf("\n");
            flint_printf("naive oracle : "); acb_printd(ref, 15); flint_printf("\n");
            acb_clear(ref);
            TEST_FUNCTION_FAIL("char %wd: acb_theta_all disagrees with naive theta^2\n", k);
        }
        acb_clear(ref);
    }

    _acb_vec_clear(th, 16);
    _acb_vec_clear(z, 2);
    acb_mat_clear(tau);

    TEST_FUNCTION_END(state);
}
