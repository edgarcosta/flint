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

TEST_FUNCTION_START(acb_ppav_g2_theta2, state)
{
    /* The ten even acb_theta indices; the other six are odd (identically 0). */
    const slong evk[10] = {0, 1, 2, 3, 4, 6, 8, 9, 12, 15};
    const slong oddk[6] = {5, 7, 10, 11, 13, 14};
    /* Exact rational squared-theta table for (l, m, n) = (2, 3, 5):
       th4[k] = th2[k]^2.  The th2 quotients are the square roots of the th4
       quotients, so squaring is exact up to enclosure.  Panel-verified;
       re-derivable from Thomae's formulae. */
    const slong t4num[10] = {1, 1, 3, 8, 2, 1, 3, 1, 1, 1};
    const slong t4den[10] = {1, 2, 5, 15, 5, 6, 10, 15, 3, 10};
    const slong precs[2] = {128, 256};
    slong ip, i, k;

    /* (A) Squaring self-consistency on the exact rational table, two precisions.
       Decoupled from the g2 fourth-power routine: th4 is hand-built here, never
       produced by another entry point.  With neg = signs = 0 every fundamental
       takes the principal root, so th2[k] = sqrt(th4[k]) and sqr(th2[k]) = th4[k]
       exactly up to enclosure.  Catches a wrong Thomae quotient or an A[]
       misalignment with no tau and no oracle call. */
    for (ip = 0; ip < 2; ip++)
    {
        slong prec = precs[ip];
        acb_ptr th4 = _acb_vec_init(16);
        acb_ptr th2 = _acb_vec_init(16);
        acb_ptr ros = _acb_vec_init(3);
        acb_t sq, one;

        acb_init(sq);
        acb_init(one);
        acb_one(one);

        for (i = 0; i < 10; i++)
        {
            acb_set_si(&th4[evk[i]], t4num[i]);
            acb_div_si(&th4[evk[i]], &th4[evk[i]], t4den[i], prec);
        }
        acb_set_si(&ros[0], 2);
        acb_set_si(&ros[1], 3);
        acb_set_si(&ros[2], 5);

        acb_ppav_g2_theta2(th2, th4, ros, 0, 0, prec);

        if (acb_contains_zero(&th2[0]))
            TEST_FUNCTION_FAIL("(A) prec %wd: th2[0] must be certainly nonzero: %{acb}\n",
                prec, &th2[0]);

        /* pin the value +1, not merely nonzero: the routine sets th2[0] = 1 */
        if (!acb_contains(&th2[0], one))
            TEST_FUNCTION_FAIL("(A) prec %wd: th2[0] must contain +1: %{acb}\n",
                prec, &th2[0]);

        for (i = 0; i < 10; i++)
        {
            k = evk[i];
            if (!acb_is_finite(&th2[k]))
                TEST_FUNCTION_FAIL("(A) prec %wd: th2[%wd] not finite: %{acb}\n",
                    prec, k, &th2[k]);

            acb_sqr(sq, &th2[k], prec);
            if (!acb_overlaps(sq, &th4[k]))
                TEST_FUNCTION_FAIL("(A) prec %wd: sqr(th2[%wd]) = %{acb} does not overlap th4 = %{acb}\n",
                    prec, k, sq, &th4[k]);

            /* nonvacuity floor: a fat ball must not pass */
            if (acb_rel_accuracy_bits(&th2[k]) < prec - 30)
                TEST_FUNCTION_FAIL("(A) prec %wd: th2[%wd] fat, acc %wd < %wd: %{acb}\n",
                    prec, k, acb_rel_accuracy_bits(&th2[k]), prec - 30, &th2[k]);
        }

        for (i = 0; i < 6; i++)
        {
            if (!acb_is_zero(&th2[oddk[i]]))
                TEST_FUNCTION_FAIL("(A) prec %wd: odd th2[%wd] must be exactly zero: %{acb}\n",
                    prec, oddk[i], &th2[oddk[i]]);
        }

        acb_clear(sq);
        acb_clear(one);
        _acb_vec_clear(th4, 16);
        _acb_vec_clear(th2, 16);
        _acb_vec_clear(ros, 3);
    }

    /* (B) Sign/neg selection on synthetic th4 (teeth for the sqrt machinery).
       ros = (2, 3, 5) is kept consistent-enough only so the derived entries do
       not divide by zero; the derived entries are NOT asserted on here, only the
       fundamentals we directly control. */
    {
        slong prec = 128;
        acb_ptr th4 = _acb_vec_init(16);
        acb_ptr th2 = _acb_vec_init(16);
        acb_ptr ros = _acb_vec_init(3);
        acb_t saved, psqrt, sq;
        arb_t err;

        acb_init(saved);
        acb_init(psqrt);
        acb_init(sq);
        arb_init(err);

        for (i = 0; i < 10; i++)
        {
            acb_set_si(&th4[evk[i]], t4num[i]);
            acb_div_si(&th4[evk[i]], &th4[evk[i]], t4den[i], prec);
        }
        acb_set_si(&ros[0], 2);
        acb_set_si(&ros[1], 3);
        acb_set_si(&ros[2], 5);

        /* (B)(i) Fundamental kf = 1 is th4[1] = 1/2 > 0, neg = 0.
           signs = 0: principal root, certainly-positive real part.
           signs = 1 (bit j = 0 controls kf[0] = 1): the negation. */
        acb_ppav_g2_theta2(th2, th4, ros, 0, 0, prec);
        /* Guard the fundamental before the real-part sign tests below.
           arb_is_positive/negative constrain only the real part, so a fat
           imaginary component could make the negation overlap vacuous; the
           acb-level rel-accuracy floor accounts for both components and closes
           that gap.  acb_sqrt of this exactly-real positive fixture keeps the
           imaginary part exactly zero, so th2[1] is also certainly real. */
        if (!acb_is_finite(&th2[1]))
            TEST_FUNCTION_FAIL("(B)(i) neg=0 signs=0: th2[1] not finite: %{acb}\n",
                &th2[1]);
        if (acb_rel_accuracy_bits(&th2[1]) < prec - 30)
            TEST_FUNCTION_FAIL("(B)(i) neg=0 signs=0: th2[1] fat, acc %wd < %wd: %{acb}\n",
                acb_rel_accuracy_bits(&th2[1]), prec - 30, &th2[1]);
        if (!acb_is_real(&th2[1]))
            TEST_FUNCTION_FAIL("(B)(i) neg=0 signs=0: th2[1] not real: %{acb}\n",
                &th2[1]);
        if (!arb_is_positive(acb_realref(&th2[1])))
            TEST_FUNCTION_FAIL("(B)(i) neg=0 signs=0: root at 1 not certainly positive: %{acb}\n",
                &th2[1]);
        acb_set(saved, &th2[1]);

        acb_ppav_g2_theta2(th2, th4, ros, 0, 1, prec);
        if (!acb_is_finite(&th2[1]))
            TEST_FUNCTION_FAIL("(B)(i) neg=0 signs=1: th2[1] not finite: %{acb}\n",
                &th2[1]);
        if (acb_rel_accuracy_bits(&th2[1]) < prec - 30)
            TEST_FUNCTION_FAIL("(B)(i) neg=0 signs=1: th2[1] fat, acc %wd < %wd: %{acb}\n",
                acb_rel_accuracy_bits(&th2[1]), prec - 30, &th2[1]);
        if (!acb_is_real(&th2[1]))
            TEST_FUNCTION_FAIL("(B)(i) neg=0 signs=1: th2[1] not real: %{acb}\n",
                &th2[1]);
        if (!arb_is_negative(acb_realref(&th2[1])))
            TEST_FUNCTION_FAIL("(B)(i) neg=0 signs=1: root at 1 not certainly negative: %{acb}\n",
                &th2[1]);
        acb_neg(saved, saved);
        if (!acb_overlaps(&th2[1], saved))
            TEST_FUNCTION_FAIL("(B)(i) neg=0 signs=1: root at 1 = %{acb} does not overlap negation %{acb}\n",
                &th2[1], saved);

        /* (B)(ii) Straddling branch-cut entry at kf = 2 (so j = 1):
           th4[2] = -1 + i*[-eps, eps] with eps = 2^-prec.  eps is scaled to prec
           (not a fixed 1e-10) so the rotated root can keep ~prec bits: a fixed
           1e-10 would cap the input relative accuracy near 33 bits, far below
           prec - 30, and the rotated-root accuracy assertion could not hold.
           The center is a nonzero negative real; the imaginary interval contains
           0, so the principal acb_sqrt straddles the branch cut and goes fat,
           while i*sqrt(-z) stays away from it and stays sharp. */
        acb_set_si(&th4[2], -1);
        arb_one(err);
        arb_mul_2exp_si(err, err, -prec);
        arb_add_error(acb_imagref(&th4[2]), err);

        acb_sqrt(psqrt, &th4[2], prec);
        if (acb_rel_accuracy_bits(psqrt) >= prec - 30)
            TEST_FUNCTION_FAIL("(B)(ii) principal sqrt unexpectedly sharp, acc %wd >= %wd: %{acb}\n",
                acb_rel_accuracy_bits(psqrt), prec - 30, psqrt);

        /* bit POSITION kf = 2 of neg means value 1 << 2 = 4 */
        acb_ppav_g2_theta2(th2, th4, ros, 4, 0, prec);
        if (!arb_is_nonnegative(acb_imagref(&th2[2])))
            TEST_FUNCTION_FAIL("(B)(ii) neg=4: rotated root at 2 imag not certainly >= 0: %{acb}\n",
                &th2[2]);
        if (acb_rel_accuracy_bits(&th2[2]) < prec - 30)
            TEST_FUNCTION_FAIL("(B)(ii) neg=4: rotated root at 2 fat, acc %wd < %wd: %{acb}\n",
                acb_rel_accuracy_bits(&th2[2]), prec - 30, &th2[2]);
        acb_sqr(sq, &th2[2], prec);
        if (!acb_overlaps(sq, &th4[2]))
            TEST_FUNCTION_FAIL("(B)(ii) neg=4: sqr(rotated root) = %{acb} does not overlap th4[2] = %{acb}\n",
                sq, &th4[2]);

        arb_clear(err);
        acb_clear(saved);
        acb_clear(psqrt);
        acb_clear(sq);
        _acb_vec_clear(th4, 16);
        _acb_vec_clear(th2, 16);
        _acb_vec_clear(ros, 3);
    }

    TEST_FUNCTION_END(state);
}
