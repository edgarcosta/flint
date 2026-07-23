/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_ppav.h"

TEST_FUNCTION_START(acb_ppav_g2_siegel_coset_nb, state)
{
    slong ells[] = {2, 3, 5, 7, 101};
    slong nb = sizeof(ells) / sizeof(*ells);
    slong ix;

    for (ix = 0; ix < nb; ix++)
    {
        slong ell = ells[ix];
        slong expected = ell * ell * ell + ell * ell + ell + 1;
        slong got = acb_ppav_g2_siegel_coset_nb(ell);
        slong got2 = acb_ppav_g2_siegel_2step_coset_nb(ell);

        if (got != expected)
            TEST_FUNCTION_FAIL("ell = %wd: siegel_coset_nb = %wd, expected %wd\n",
                ell, got, expected);

        if (got2 != ell * got)
            TEST_FUNCTION_FAIL("ell = %wd: siegel_2step_coset_nb = %wd, expected %wd\n",
                ell, got2, ell * got);
    }

    TEST_FUNCTION_END(state);
}
