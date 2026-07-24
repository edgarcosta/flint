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
#include "acb_poly.h"
#include "acb_ppav.h"

/* f = prod_i (x - roots[i]); integer roots keep every coefficient an
   exact integer at any prec, so the returned balls can be tested for
   containment against exact-integer targets. */
static void
_weierstrass_poly_from_int_roots(acb_poly_t f, const slong * roots, slong d,
    slong prec)
{
    acb_poly_t lin;
    slong i;

    acb_poly_init(lin);
    acb_poly_set_coeff_si(lin, 1, 1);
    acb_poly_one(f);
    for (i = 0; i < d; i++)
    {
        acb_poly_set_coeff_si(lin, 0, -roots[i]);
        acb_poly_mul(f, f, lin, prec);
    }
    acb_poly_clear(lin);
}

/* Bijection between the d returned balls and the d integer roots: each
   root sits in exactly one ball and each ball holds exactly one root.
   Compared as sets, since find_roots output order is not canonical.
   Returns 1 on a clean set-match, 0 otherwise. */
static int
_weierstrass_set_match(acb_srcptr w, const slong * roots, slong d)
{
    acb_t t;
    slong i, j;
    int ok = 1;

    acb_init(t);
    for (i = 0; i < d; i++)
    {
        slong cnt = 0;
        acb_set_si(t, roots[i]);
        for (j = 0; j < d; j++)
            if (acb_contains(&w[j], t))
                cnt++;
        if (cnt != 1)
            ok = 0;
    }
    for (j = 0; j < d; j++)
    {
        slong cnt = 0;
        for (i = 0; i < d; i++)
        {
            acb_set_si(t, roots[i]);
            if (acb_contains(&w[j], t))
                cnt++;
        }
        if (cnt != 1)
            ok = 0;
    }
    acb_clear(t);
    return ok;
}

/* Minimum relative accuracy (bits) over the d returned balls. */
static slong
_weierstrass_min_rel_acc(acb_srcptr w, slong d)
{
    slong i, m = acb_rel_accuracy_bits(&w[0]);

    for (i = 1; i < d; i++)
    {
        slong a = acb_rel_accuracy_bits(&w[i]);
        if (a < m)
            m = a;
    }
    return m;
}

TEST_FUNCTION_START(acb_ppav_weierstrass, state)
{
    slong precs[3] = {64, 128, 256};
    const slong rootsA[6] = {-1, 1, 2, 3, 4, 6};    /* (A) g=2 squarefree sextic */
    const slong rootsB[3] = {-1, 2, 5};             /* (B) g=1 cubic */
    const slong rootsC[6] = {1, 1, 2, 3, 4, 5};     /* (C) double root, inseparable */
    const slong rootsD5[5] = {1, 2, 3, 4, 5};       /* (D i) degree 5, wrong for g=2 */
    const slong rootsD7[7] = {1, 2, 3, 4, 5, 6, 7}; /* (D ii) degree 7, overrun canary */
    slong iter;

    for (iter = 0; iter < 3; iter++)
    {
        slong prec = precs[iter];

        /* (A) g=2 squarefree sextic: all six roots isolated and
           set-matched to {-1,1,2,3,4,6}, no zero root so relative
           accuracy is meaningful and must scale with prec. */
        {
            acb_poly_t f;
            acb_ptr w = _acb_vec_init(6);

            acb_poly_init(f);
            _weierstrass_poly_from_int_roots(f, rootsA, 6, prec);
            if (acb_ppav_weierstrass(w, f, 2, prec) != 1)
                TEST_FUNCTION_FAIL("(A) prec=%wd: expected return 1\n", prec);
            if (!_weierstrass_set_match(w, rootsA, 6))
                TEST_FUNCTION_FAIL("(A) prec=%wd: roots do not set-match "
                    "{-1,1,2,3,4,6}\n", prec);
            if (_weierstrass_min_rel_acc(w, 6) < prec - 30)
                TEST_FUNCTION_FAIL("(A) prec=%wd: rel accuracy %wd < %wd\n",
                    prec, _weierstrass_min_rel_acc(w, 6), prec - 30);
            acb_poly_clear(f);
            _acb_vec_clear(w, 6);
        }

        /* (B) g=1 cubic: three roots set-matched to {-1,2,5}. */
        {
            acb_poly_t f;
            acb_ptr w = _acb_vec_init(3);

            acb_poly_init(f);
            _weierstrass_poly_from_int_roots(f, rootsB, 3, prec);
            if (acb_ppav_weierstrass(w, f, 1, prec) != 1)
                TEST_FUNCTION_FAIL("(B) prec=%wd: expected return 1\n", prec);
            if (!_weierstrass_set_match(w, rootsB, 3))
                TEST_FUNCTION_FAIL("(B) prec=%wd: roots do not set-match "
                    "{-1,2,5}\n", prec);
            if (_weierstrass_min_rel_acc(w, 3) < prec - 30)
                TEST_FUNCTION_FAIL("(B) prec=%wd: rel accuracy %wd < %wd\n",
                    prec, _weierstrass_min_rel_acc(w, 3), prec - 30);
            acb_poly_clear(f);
            _acb_vec_clear(w, 3);
        }

        /* (C) non-squarefree sextic: the double root at 1 is an
           inseparable cluster that find_roots cannot split at any finite
           prec, so nothing is certified -> return 0. */
        {
            acb_poly_t f;
            acb_ptr w = _acb_vec_init(6);

            acb_poly_init(f);
            _weierstrass_poly_from_int_roots(f, rootsC, 6, prec);
            if (acb_ppav_weierstrass(w, f, 2, prec) != 0)
                TEST_FUNCTION_FAIL("(C) prec=%wd: double root must return 0\n",
                    prec);
            acb_poly_clear(f);
            _acb_vec_clear(w, 6);
        }

        /* (D i) degree 5 with g=2 (expected degree 6): degree guard -> 0. */
        {
            acb_poly_t f;
            acb_ptr w = _acb_vec_init(6);

            acb_poly_init(f);
            _weierstrass_poly_from_int_roots(f, rootsD5, 5, prec);
            if (acb_ppav_weierstrass(w, f, 2, prec) != 0)
                TEST_FUNCTION_FAIL("(D i) prec=%wd: degree 5 with g=2 must "
                    "return 0\n", prec);
            acb_poly_clear(f);
            _acb_vec_clear(w, 6);
        }

        /* (D ii) degree 7 with g=2: the guard must reject before
           find_roots writes deg(f)=7 entries into a length-6 buffer.
           A sentinel one slot past w[5] must survive untouched. */
        {
            acb_poly_t f;
            acb_ptr w = _acb_vec_init(7);
            acb_t canary;

            acb_init(canary);
            acb_set_si(canary, 987654321);
            acb_mul_onei(canary, canary);       /* pure-imaginary, never a root */
            acb_set(&w[6], canary);

            acb_poly_init(f);
            _weierstrass_poly_from_int_roots(f, rootsD7, 7, prec);
            if (acb_ppav_weierstrass(w, f, 2, prec) != 0)
                TEST_FUNCTION_FAIL("(D ii) prec=%wd: degree 7 with g=2 must "
                    "return 0\n", prec);
            if (!acb_equal(&w[6], canary))
                TEST_FUNCTION_FAIL("(D ii) prec=%wd: buffer overrun past "
                    "w[5]: %{acb}\n", prec, &w[6]);
            acb_clear(canary);
            acb_poly_clear(f);
            _acb_vec_clear(w, 7);
        }

        /* (D iii) invalid genus g=3 on a valid degree-6 poly: genus
           guard -> 0. */
        {
            acb_poly_t f;
            acb_ptr w = _acb_vec_init(6);

            acb_poly_init(f);
            _weierstrass_poly_from_int_roots(f, rootsA, 6, prec);
            if (acb_ppav_weierstrass(w, f, 3, prec) != 0)
                TEST_FUNCTION_FAIL("(D iii) prec=%wd: genus 3 must return 0\n",
                    prec);
            acb_poly_clear(f);
            _acb_vec_clear(w, 6);
        }
    }

    TEST_FUNCTION_END(state);
}
