/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_mulmod_shoup, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong a, b, d, r1, r2, q, p1, p2, w_pr;

        d = n_randtest_not_zero(state) / 2 + 1;
        a = n_randtest(state) % d;
        b = n_randtest(state) % d;

        w_pr = n_mulmod_precomp_shoup(a, d);

        r1 = n_mulmod_shoup(a, b, w_pr, d);

        umul_ppmm(p1, p2, a, b);
        p1 %= d;
        udiv_qrnnd(q, r2, p1, p2, d);

        result = (r1 == r2);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "a = %wu, b = %wu, d = %wu, w_pr = %wu\n"
                    "q = %wu, r1 = %wu, r2 = %wu\n",
                    a, b, d, w_pr, q, r1, r2);
    }

    TEST_FUNCTION_END(state);
}
