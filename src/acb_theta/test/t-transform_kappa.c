/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("transform_kappa....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: kappa^2 on [u, 0; 0, u^-t] is det(u) */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        fmpz_mat_t U, mat;
        fmpz_t det;
        slong kappa;
        slong bits = 1 + n_randint(state, 10);

        fmpz_mat_init(U, g, g);
        fmpz_mat_init(mat, 2 * g, 2 * g);
        fmpz_init(det);

        fmpz_mat_one(U);
        if (iter % 2 == 0)
        {
            fmpz_set_si(fmpz_mat_entry(U, 0, 0), -1);
        }

        fmpz_mat_randops(U, state, 2 * bits);
        sp2gz_block_diag(mat, U);
        fmpz_mat_det(det, U);
        kappa = acb_theta_transform_kappa(mat);

        /* det is 1 or -1; k2 is 0 or 2 mod 4 */
        if ((kappa % 4) != 1 - fmpz_get_si(det))
        {
            flint_printf("FAIL\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("\n");
            flint_printf("k2: %wd\n", kappa);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(U);
        fmpz_mat_clear(mat);
        fmpz_clear(det);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}