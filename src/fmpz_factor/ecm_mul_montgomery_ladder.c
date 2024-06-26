/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "ulong_extras.h"
#include "fmpz_factor.h"

/* P (x0 : z0) <- kP using Montgomery ladder algorithm */

/* tstbit uses 0 based indexing */

void
fmpz_factor_ecm_mul_montgomery_ladder(nn_ptr x, nn_ptr z, nn_ptr x0, nn_ptr z0,
                                      ulong k, nn_ptr n, ecm_t ecm_inf)
{
    nn_ptr x1, z1, x2, z2;      /* Q (x1 : z1), P (x2 : z2) */
    ulong len;

    TMP_INIT;

    if (k == 0)
    {
        mpn_zero(x, ecm_inf->n_size);
        mpn_zero(z, ecm_inf->n_size);
        return;
    }

    if (k == 1)
    {
        flint_mpn_copyi(x, x0, ecm_inf->n_size);
        flint_mpn_copyi(z, z0, ecm_inf->n_size);
        return;
    }

    TMP_START;
    x1 = TMP_ALLOC(ecm_inf->n_size * sizeof(ulong));
    z1 = TMP_ALLOC(ecm_inf->n_size * sizeof(ulong));
    x2 = TMP_ALLOC(ecm_inf->n_size * sizeof(ulong));
    z2 = TMP_ALLOC(ecm_inf->n_size * sizeof(ulong));


    flint_mpn_copyi(x1, x0, ecm_inf->n_size);    /* Q <- P0 */
    flint_mpn_copyi(z1, z0, ecm_inf->n_size);
    mpn_zero(x2, ecm_inf->n_size);
    mpn_zero(z2, ecm_inf->n_size);

    fmpz_factor_ecm_double(x2, z2, x0, z0, n, ecm_inf);   /* P <- 2P0 */

    len = n_sizeinbase(k, 2) - 2;

    while (1)
    {
        if (((UWORD(1) << len) & k) != 0)       /* ith bit is 1 */
        {
            /* Q <- P + Q */
            fmpz_factor_ecm_add(x1, z1, x1, z1, x2, z2, x0, z0, n, ecm_inf);
            /* P <- 2 * P */
            fmpz_factor_ecm_double(x2, z2, x2, z2, n, ecm_inf);
        }
        else
        {
            /* P <- P + Q */
            fmpz_factor_ecm_add(x2, z2, x1, z1, x2, z2, x0, z0, n, ecm_inf);
            /* Q <- 2 * Q */
            fmpz_factor_ecm_double(x1, z1, x1, z1, n, ecm_inf);
        }

        if (len == 0)
            break;
        else
            len -= 1;
    }

    flint_mpn_copyi(x, x1, ecm_inf->n_size);
    flint_mpn_copyi(z, z1, ecm_inf->n_size);

    TMP_END;
}
