/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_is_probabprime, state)
{
    int i, result;
    ulong d;
    mpz_t d_m;
    slong pow;
    ulong bits;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that primes pass the test */
    {
        mpz_init(d_m);

        do
        {
            d = n_randtest(state) | 1;
            flint_mpz_set_ui(d_m, d);
            mpz_nextprime(d_m, d_m);
            d = flint_mpz_get_ui(d_m);
        } while (mpz_size(d_m) > 1);

        result = n_is_probabprime(d);
        if (!result)
            TEST_FUNCTION_FAIL("d = %wu is declared composite\n", d);

        mpz_clear(d_m);
    }

    for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that composites do not pass */
    {
        mpz_init(d_m);

        do
        {
            d = n_randtest(state) | 1;
            if (d == UWORD(1)) d++;
            flint_mpz_set_ui(d_m, d);
        } while (mpz_probab_prime_p(d_m, 12));

        result = !n_is_probabprime(d);
        if (!result)
            TEST_FUNCTION_FAIL("d = %wu is declared prime\n", d);

        mpz_clear(d_m);
    }

    /* Test that powers do not pass */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        pow = n_randint(state, 6) + 2;
        bits = n_randint(state, FLINT_BITS) + 1;
        bits /= pow;

        d = n_randbits(state, bits);
        d = n_pow(d, pow);

        result = !n_is_probabprime(d);
        if (!result)
            TEST_FUNCTION_FAIL("Perfect power d = %wu is declared prime\n", d);
    }

    /* Regression test, check certain composites do not pass */
#if FLINT64
    {
        d = UWORD(2007193456621);

        result = !n_is_probabprime(d);
        if (!result)
            TEST_FUNCTION_FAIL("Known composite d = %wu is declared prime\n", d);
    }
#endif

    TEST_FUNCTION_END(state);
}
