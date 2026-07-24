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
#include "acb_mat.h"
#include "acb_theta.h"
#include "acb_ppav.h"
#include "acb_ppav_thomae_ref.h"

/* Exact rational ball p/q: p and q are exact, one rounding in the division. */
static void
_acb_set_si_div_si(acb_t res, slong p, slong q, slong prec)
{
    acb_set_si(res, p);
    acb_div_si(res, res, q, prec);
}

/* Case (C)'s input-only conditioning gate, leg (b').  The nine Thomae factors
   {l, m, n, l-1, m-1, n-1, l-m, l-n, m-n} recovered from ros must each be
   finite, certainly nonzero, and tight (rel_acc >= prec-25).  Relative accuracy,
   not absolute magnitude, is the right control: it is exactly what the th4 floor
   measures, so a large-magnitude difference that is only relatively near zero
   (cancellation) is correctly rejected while a tiny-but-sharp one is kept.  Reads
   only ros (test-side); never th4 or theta4's return value. */
static int
_factors_well_conditioned(acb_srcptr ros, slong prec)
{
    acb_t f;
    slong i;
    int ok = 1;

    acb_init(f);
    for (i = 0; i < 9 && ok; i++)
    {
        switch (i)
        {
            case 0: acb_set(f, &ros[0]); break;                /* l   */
            case 1: acb_set(f, &ros[1]); break;                /* m   */
            case 2: acb_set(f, &ros[2]); break;                /* n   */
            case 3: acb_sub_si(f, &ros[0], 1, prec); break;    /* l-1 */
            case 4: acb_sub_si(f, &ros[1], 1, prec); break;    /* m-1 */
            case 5: acb_sub_si(f, &ros[2], 1, prec); break;    /* n-1 */
            case 6: acb_sub(f, &ros[0], &ros[1], prec); break; /* l-m */
            case 7: acb_sub(f, &ros[0], &ros[2], prec); break; /* l-n */
            case 8: acb_sub(f, &ros[1], &ros[2], prec); break; /* m-n */
        }
        ok = acb_is_finite(f) && !acb_contains_zero(f)
             && acb_rel_accuracy_bits(f) >= prec - 25;
    }
    acb_clear(f);
    return ok;
}

TEST_FUNCTION_START(acb_ppav_g2_theta4, state)
{
    slong prec = 200;

    /* (A) placement + normalization and (B) tau-free Thomae oracle, both at the
       exact rational Rosenhain triple (l, m, n) = (2, 3, 5).  The ten even
       acb_theta indices carry exactly (hand-checked against hdme thomae_theta4.c
       with the A[] label map):
         k      0    1    2     3    4    6     8    9    12   15
         th4   1   1/2  3/5  8/15  2/5 1/6  3/10 1/15 1/3  1/10  */
    {
        acb_ptr ros = _acb_vec_init(3);
        acb_ptr th4 = _acb_vec_init(16);
        acb_t one, expected;
        slong k;
        const slong idx[10] = { 0, 1, 2,  3, 4, 6,  8,  9, 12, 15 };
        const slong num[10] = { 1, 1, 3,  8, 2, 1,  3,  1,  1,  1 };
        const slong den[10] = { 1, 2, 5, 15, 5, 6, 10, 15,  3, 10 };

        acb_init(one);
        acb_init(expected);
        acb_one(one);

        acb_set_si(&ros[0], 2);
        acb_set_si(&ros[1], 3);
        acb_set_si(&ros[2], 5);

        acb_ppav_g2_theta4(th4, ros, prec);

        /* (A) th4[0] contains 1 (acb_contains, NOT acb_contains_int) */
        if (!acb_contains(&th4[0], one))
        {
            flint_printf("FAIL (A): th4[0] does not contain 1\n");
            flint_abort();
        }

        /* (A) even indices finite and certainly nonzero, odd indices exactly 0 */
        for (k = 0; k < 16; k++)
        {
            if (acb_theta_char_is_even(k, 2))
            {
                if (!acb_is_finite(&th4[k]) || acb_contains_zero(&th4[k]))
                {
                    flint_printf("FAIL (A): even th4[%wd] not finite-nonzero\n", k);
                    flint_abort();
                }
            }
            else if (!acb_is_zero(&th4[k]))
            {
                flint_printf("FAIL (A): odd th4[%wd] not exactly zero\n", k);
                flint_abort();
            }
        }

        /* (B) exact ten-value oracle + accuracy (straight-line field arithmetic
           on exact input must not decay) */
        for (k = 0; k < 10; k++)
        {
            _acb_set_si_div_si(expected, num[k], den[k], prec);
            if (!acb_contains(&th4[idx[k]], expected))
            {
                flint_printf("FAIL (B): th4[%wd] != %wd/%wd\n",
                    idx[k], num[k], den[k]);
                flint_abort();
            }
            if (acb_rel_accuracy_bits(&th4[idx[k]]) < prec - 30)
            {
                flint_printf("FAIL (B): th4[%wd] lost accuracy\n", idx[k]);
                flint_abort();
            }
        }

        acb_clear(one);
        acb_clear(expected);
        _acb_vec_clear(ros, 3);
        _acb_vec_clear(th4, 16);
    }

    /* (C) acb_theta cross-check, projective, deterministic random reduced tau.
       theta^2 via acb_theta_all(sqr=1); recover (l,m,n) with the shared
       inverse-Thomae helper; call theta4; then check the ten even entries
       satisfy th4[k]/th4[0] == (th2[k]/th2[0])^2 (cross-multiplied, so no
       division by an uncertified ball).  Squaring th2 once yields theta^4 since
       acb_theta_all(sqr=1) already returns theta-SQUARED.  Both cross-multiplied
       sides rigorously enclose theta[k]^4, so a correct theta4 always overlaps
       regardless of conditioning; the overlap bites only a wrong formula (teeth).

       CONDITIONING GATE (input-only cutoff, not a correctness claim).  A reduced
       draw is gate-accepted iff (a) all ten even reference th2[k] are finite and
       certainly nonzero, and (b') each of the nine Thomae factors
       {l,m,n,l-1,m-1,n-1,l-m,l-n,m-n} recovered from ros is finite, certainly
       nonzero and tight to prec-25 (leg (b'), the helper above).  Calibration of
       the prec-25 bound: every even th4[k] is a product/quotient of these factors,
       and ball mul/div roughly adds relative errors, but the chained derived
       labels (3, 6, 9, 12, 15) accumulate the summed errors of several factors
       (worst case, label 15), so the analytic derivation margin at the prec-30
       th4 floor is near zero -- NOT a uniform ~2 bits.  prec-25 is chosen so the
       worst-case chained label just clears prec-30; the actual slack on the
       pinned seed is the measured +2 bits (th4 side).  The gate reads only ros
       and th2 -- never th4 or theta4's return value, so a fat or wrong theta4
       cannot buy itself a skip. */
    {
        slong ntry, n_accepted = 0;
        acb_mat_t tau;
        acb_ptr z, th2, ros, th4;
        acb_t lhs, rhs, t0sq, tksq;
        const slong even_idx[10] = { 0, 1, 2, 3, 4, 6, 8, 9, 12, 15 };

        acb_mat_init(tau, 2, 2);
        z = _acb_vec_init(2);
        th2 = _acb_vec_init(16);
        ros = _acb_vec_init(3);
        th4 = _acb_vec_init(16);
        acb_init(lhs);
        acb_init(rhs);
        acb_init(t0sq);
        acb_init(tksq);

        for (ntry = 0; ntry < 50; ntry++)
        {
            slong j;
            int gate = 1;

            /* LANDMINE: acb_siegel_randtest_reduced violates its own -20
               precondition ~30% of the time; re-check and skip failures. */
            acb_siegel_randtest_reduced(tau, state, prec, 2);
            if (!acb_siegel_is_reduced(tau, -20, prec))
                continue;

            acb_theta_all(th2, z, tau, 1, prec);  /* theta-SQUARED */

            /* Gate leg (a): ten even reference th2[k] finite and certainly
               nonzero (pins the reference side of every overlap). */
            for (j = 0; j < 10 && gate; j++)
                gate = acb_is_finite(&th2[even_idx[j]])
                       && !acb_contains_zero(&th2[even_idx[j]]);

            _ros_from_theta2(ros, th2, prec);

            /* Gate leg (b'): nine recovered Thomae factors tight to prec-25. */
            if (gate)
                gate = _factors_well_conditioned(ros, prec);

            acb_ppav_g2_theta4(th4, ros, prec);

            /* Skip an undefined projective comparison only on NON-gated
               (bonus-coverage) draws whose k=0 normalizer th2[0] or th4[0] is
               non-finite or contains zero.  The `!gate` guard is essential: on a
               gate-accepted draw nothing is skipped -- leg (a) already guarantees
               th2[0] finite and certainly nonzero, and the hard-floor loop below
               (which covers k = 0) aborts on any bad th4[0].  So the normalizer
               preconditions are subsumed by the floors on the gated path, and
               reading th4[0] here can never let a broken theta4 buy itself a skip
               (the output-buys-a-skip circularity the input-only gate removes). */
            if (!gate && (!acb_is_finite(&th2[0]) || acb_contains_zero(&th2[0])
                || !acb_is_finite(&th4[0]) || acb_contains_zero(&th4[0])))
            {
                continue;
            }

            /* Gate-accepted: hard floor on every even th4[k] (finite, certainly
               nonzero, tight).  A well-conditioned input must give a tight
               theta4; anything less is a theta4 accuracy bug, never a skip. */
            if (gate)
            {
                for (j = 0; j < 10; j++)
                {
                    slong k = even_idx[j];
                    slong acc = acb_rel_accuracy_bits(&th4[k]);

                    if (!acb_is_finite(&th4[k]) || acb_contains_zero(&th4[k])
                        || acc < prec - 30)
                    {
                        flint_printf("FAIL (C): th4[%wd] not finite-nonzero-tight "
                            "(rel_acc %wd < %wd)\n", k, acc, prec - 30);
                        acb_mat_printd(tau, 5);
                        flint_abort();
                    }
                }
            }

            /* Projective overlap on the ten even k.  Asserted on EVERY reduced
               draw with good normalizers: both sides enclose theta[k]^4, so fat
               balls only make the overlap easier, widening correctness coverage
               safely.  On gate-accepted draws BOTH cross-multiplied sides must
               also clear a vacuity floor (below) before the overlap, so no
               accepted overlap passes vacuously. */
            acb_sqr(t0sq, &th2[0], prec);
            for (j = 0; j < 10; j++)
            {
                slong k = even_idx[j];

                /* th4[k] * th2[0]^2  overlaps  th4[0] * th2[k]^2 */
                acb_sqr(tksq, &th2[k], prec);
                acb_mul(lhs, &th4[k], t0sq, prec);
                acb_mul(rhs, &th4[0], tksq, prec);

                if (gate)
                {
                    slong al = acb_rel_accuracy_bits(lhs);
                    slong ar = acb_rel_accuracy_bits(rhs);

                    /* Vacuity guard, NOT the accuracy certification (that is the
                       th4-side floor above at prec-30, left untouched).  Set
                       deliberately below it at prec-40 to be robust to platform
                       and library drift: a prec-40-accurate ball is still nowhere
                       near vacuous, so the overlap stays a real agreement. */
                    if (al < prec - 40 || ar < prec - 40)
                    {
                        flint_printf("FAIL (C): th4[%wd] cross-mult side lost "
                            "accuracy (lhs %wd, rhs %wd, floor %wd)\n",
                            k, al, ar, prec - 40);
                        acb_mat_printd(tau, 5);
                        flint_abort();
                    }
                }

                if (!acb_overlaps(lhs, rhs))
                {
                    flint_printf("FAIL (C): th4[%wd] disagrees with acb_theta\n", k);
                    acb_mat_printd(tau, 5);
                    flint_abort();
                }
            }
            if (gate)
                n_accepted++;
        }

        /* Guard against a vacuous pass that gate-accepted too few draws. */
        if (n_accepted < 5)
        {
            flint_printf("FAIL (C): too few accepted draws (%wd)\n", n_accepted);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, 2);
        _acb_vec_clear(th2, 16);
        _acb_vec_clear(ros, 3);
        _acb_vec_clear(th4, 16);
        acb_clear(lhs);
        acb_clear(rhs);
        acb_clear(t0sq);
        acb_clear(tksq);
    }

    /* (D) branch-cut flag.  Fixture (l, m, n) = (2, 5, 3) makes label 2 the
       exact negative rational
         th4[1] = m(l-1)(n-m) / (l(m-1)(n-l))
                = 5*1*(-2) / (2*4*1) = -10/8 = -5/4  < 0,
       a certain negative real (imag exactly 0), so bit 1 must be SET; and
         th4[8] = m/(l n) = 5/6 > 0,
       a certain positive real, so bit 8 must be CLEAR. */
    {
        acb_ptr ros = _acb_vec_init(3);
        acb_ptr th4 = _acb_vec_init(16);
        acb_t neg;
        slong flag;

        acb_init(neg);

        acb_set_si(&ros[0], 2);
        acb_set_si(&ros[1], 5);
        acb_set_si(&ros[2], 3);

        flag = acb_ppav_g2_theta4(th4, ros, prec);

        _acb_set_si_div_si(neg, -5, 4, prec);
        if (!acb_contains(&th4[1], neg))
        {
            flint_printf("FAIL (D): th4[1] does not contain -5/4\n");
            flint_abort();
        }
        if (!(flag & (WORD(1) << 1)))
        {
            flint_printf("FAIL (D): branch-cut bit 1 (negative real) not set\n");
            flint_abort();
        }
        if (flag & (WORD(1) << 8))
        {
            flint_printf("FAIL (D): branch-cut bit 8 (positive real) is set\n");
            flint_abort();
        }

        acb_clear(neg);
        _acb_vec_clear(ros, 3);
        _acb_vec_clear(th4, 16);
    }

    TEST_FUNCTION_END(state);
}
