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
#include "acb_mat.h"
#include "acb_theta.h"
#include "acb_ppav.h"
#include "acb_ppav_thomae_ref.h"

/* Secondary check, one tau branch: recompute the squared theta constants with
   acb_theta_all (sqr = 1, so the output already IS theta-SQUARED, matching th2;
   no extra squaring on either side) at z = 0 and compare projectively with the
   candidate's th2 on the ten even characteristics.  Returns 1 iff, on those ten
   even characteristics, both th_fwd and th2 are finite and accurate to
   prec - 30 bits, both k = 0 normalizers are certainly nonzero, and all ten
   cross-multiplied overlaps hold.  No division anywhere: the projective compare
   is th_fwd[k] th2[0] =?= th2[k] th_fwd[0]. */
static int
_theta_all_matches(const acb_mat_t tau, acb_srcptr th2, slong prec)
{
    static const slong evenk[10] = {0, 1, 2, 3, 4, 6, 8, 9, 12, 15};
    acb_ptr z, th_fwd;
    acb_t lhs, rhs;
    slong j, k;
    int ok = 1;

    z = _acb_vec_init(2);              /* z = 0 */
    th_fwd = _acb_vec_init(16);
    acb_init(lhs);
    acb_init(rhs);

    acb_theta_all(th_fwd, z, tau, 1, prec);

    /* Both sides' ten even entries must be finite and accurate to prec - 30
       bits, else this branch does not count as a match: an overlap against an
       indeterminate ball would be vacuous. */
    for (j = 0; ok && j < 10; j++)
    {
        k = evenk[j];
        if (!acb_is_finite(&th_fwd[k])
            || acb_rel_accuracy_bits(&th_fwd[k]) < prec - 30
            || !acb_is_finite(&th2[k])
            || acb_rel_accuracy_bits(&th2[k]) < prec - 30)
            ok = 0;
    }

    /* both k = 0 normalizers certainly nonzero BEFORE any projective use */
    if (ok && (acb_contains_zero(&th_fwd[0]) || acb_contains_zero(&th2[0])))
        ok = 0;

    for (j = 0; ok && j < 10; j++)
    {
        k = evenk[j];
        acb_mul(lhs, &th_fwd[k], &th2[0], prec);
        acb_mul(rhs, &th2[k], &th_fwd[0], prec);
        if (!acb_overlaps(lhs, rhs))
            ok = 0;
    }

    _acb_vec_clear(z, 2);
    _acb_vec_clear(th_fwd, 16);
    acb_clear(lhs);
    acb_clear(rhs);
    return ok;
}

/* Secondary check for one collected candidate.  Reconstruct tau from the
   periods_from_theta2 output res in the encoding confirmed in
   periods_from_theta2.c: res00 = tau11, res11 = tau22, res01 = tau12^2 (the
   square, sign-free).  tau12 is then +/- sqrt(res01); build both symmetric
   candidates and prefer the acb_siegel_is_reduced(-20) branch, trying both only
   when neither is reduced.  theta^2 is invariant under tau12 -> -tau12, so both
   branches reproduce th2; reducedness is what selects the intended
   representative.  Passes iff SOME tried branch matches all ten even
   characteristics. */
static int
_secondary_ok(const acb_mat_t res, acb_srcptr th2, slong prec)
{
    acb_mat_t tauP, tauM;
    acb_t w;
    int redP, redM, ok;

    acb_mat_init(tauP, 2, 2);
    acb_mat_init(tauM, 2, 2);
    acb_init(w);

    acb_sqrt(w, acb_mat_entry(res, 0, 1), prec);

    acb_set(acb_mat_entry(tauP, 0, 0), acb_mat_entry(res, 0, 0));
    acb_set(acb_mat_entry(tauP, 1, 1), acb_mat_entry(res, 1, 1));
    acb_set(acb_mat_entry(tauP, 0, 1), w);
    acb_set(acb_mat_entry(tauP, 1, 0), w);

    acb_neg(w, w);
    acb_set(acb_mat_entry(tauM, 0, 0), acb_mat_entry(res, 0, 0));
    acb_set(acb_mat_entry(tauM, 1, 1), acb_mat_entry(res, 1, 1));
    acb_set(acb_mat_entry(tauM, 0, 1), w);
    acb_set(acb_mat_entry(tauM, 1, 0), w);

    redP = acb_siegel_is_reduced(tauP, -20, prec);
    redM = acb_siegel_is_reduced(tauM, -20, prec);

    ok = 0;
    if (redP || redM)
    {
        if (redP)
            ok = _theta_all_matches(tauP, th2, prec);
        if (!ok && redM)
            ok = _theta_all_matches(tauM, th2, prec);
    }
    else
    {
        /* neither branch reduced: try both */
        ok = _theta_all_matches(tauP, th2, prec);
        if (!ok)
            ok = _theta_all_matches(tauM, th2, prec);
    }

    acb_mat_clear(tauP);
    acb_mat_clear(tauM);
    acb_clear(w);
    return ok;
}

TEST_FUNCTION_START(acb_ppav_thomae_chain, state)
{
    /* End-to-end reconciliation of the Thomae chain against acb_theta at a
       pinned, deterministic precision (no random draws anywhere).  This loop
       previews the search that PR-3's acb_ppav_g2_periods_lowprec will run, but
       it is test-only harness code, not that production routine. */
    slong prec = 128;
    slong perm, signs, i, k;
    slong n_valid = 0;
    acb_poly_t f;
    acb_ptr w, ros, th4, th2, ros_rec;
    acb_mat_t res;

    acb_poly_init(f);
    w = _acb_vec_init(6);
    ros = _acb_vec_init(3);
    th4 = _acb_vec_init(16);
    th2 = _acb_vec_init(16);
    ros_rec = _acb_vec_init(3);
    acb_mat_init(res, 2, 2);

    /* Fixture (exact).  f = 4 g + h^2 from the model y^2 + h y = g with
       h = x^3 + x + 1 and g = -x^2 - x, i.e.
           f = x^6 + 2 x^4 + 2 x^3 - 3 x^2 - 2 x + 1.
       f mod 3 is separable, hence f is squarefree (the only property relied on).
       Coefficients constant -> leading: 1, -2, -3, 2, 2, 0, 1. */
    acb_poly_set_coeff_si(f, 0, 1);
    acb_poly_set_coeff_si(f, 1, -2);
    acb_poly_set_coeff_si(f, 2, -3);
    acb_poly_set_coeff_si(f, 3, 2);
    acb_poly_set_coeff_si(f, 4, 2);
    acb_poly_set_coeff_si(f, 5, 0);
    acb_poly_set_coeff_si(f, 6, 1);

    if (acb_ppav_weierstrass(w, f, 2, prec) != 1)
        TEST_FUNCTION_FAIL("weierstrass failed to isolate the six roots of the fixture\n");

    /* Full perm in [0, 720) x signs in [0, 16) grid.  Rosenhain and theta4 depend
       only on perm; theta2 also on signs. */
    for (perm = 0; perm < 720; perm++)
    {
        slong neg;

        acb_ppav_g2_rosenhain(ros, w, perm, prec);
        neg = acb_ppav_g2_theta4(th4, ros, prec);

        for (signs = 0; signs < 16; signs++)
        {
            int discard, code;

            acb_ppav_g2_theta2(th2, th4, ros, neg, signs, prec);

            /* Positivity screen, mirroring hdme igusa/thomae_discard.c:10-13
               EXACTLY: discard iff some k in {1, 2, 3} has certainly-nonpositive
               real part (acb_theta indices; index 3 is a derived quotient, not a
               fundamental -- HDME's comment says "fundamental" but its code
               checks 1, 2, 3).  Certain nonpositivity only: a ball straddling 0
               is KEPT, since requiring certain positivity would be strictly
               stronger and could drop the valid candidate. */
            discard = 0;
            for (k = 1; k < 4; k++)
                if (arb_is_nonpositive(acb_realref(&th2[k])))
                    discard = 1;
            if (discard)
                continue;

            code = acb_ppav_periods_from_theta2(res, th2, 2, prec);
            if (code != 1)
                continue;

            n_valid++;

            /* PRIMARY ordering witness.  Recover the Rosenhain triple from th2
               with the HARD-CODED-index inverse-Thomae helper and require it to
               overlap the ros that fed theta2 for THIS candidate.  code == 1
               alone is not a witness: the S_6 action on Weierstrass points
               mirrors the Sp(4,F_2) action on characteristics, so a wrong label
               map could be absorbed by a compensating perm; reading th2 at fixed
               indices (never through a label map) closes that hole. */

            /* arb's overlaps returns 1 on indeterminate input, so an unguarded
               recovery would overlap anything and the witness would be vacuous.
               An accepted code == 1 candidate with unhealthy theta nulls is a
               test failure, not a skippable draw: hard-assert, never skip. */
            if (!acb_is_finite(&th2[0]) || !acb_is_finite(&th2[1])
                || !acb_is_finite(&th2[4]) || !acb_is_finite(&th2[8])
                || !acb_is_finite(&th2[9]) || !acb_is_finite(&th2[12]))
                TEST_FUNCTION_FAIL("perm %wd signs %wd: a theta null read by the "
                    "inverse Thomae is not finite\n", perm, signs);
            if (acb_contains_zero(&th2[8]) || acb_contains_zero(&th2[9])
                || acb_contains_zero(&th2[12]))
                TEST_FUNCTION_FAIL("perm %wd signs %wd: an inverse-Thomae "
                    "denominator (th2[8], th2[9], th2[12]) is not certainly "
                    "nonzero\n", perm, signs);

            _ros_from_theta2(ros_rec, th2, prec);

            for (i = 0; i < 3; i++)
                if (!acb_is_finite(&ros_rec[i])
                    || acb_rel_accuracy_bits(&ros_rec[i]) < prec - 30)
                    TEST_FUNCTION_FAIL("perm %wd signs %wd: recovered Rosenhain "
                        "entry %wd is not finite to prec - 30 bits\n",
                        perm, signs, i);

            for (i = 0; i < 3; i++)
                if (!acb_overlaps(&ros_rec[i], &ros[i]))
                    TEST_FUNCTION_FAIL("perm %wd signs %wd: recovered Rosenhain "
                        "entry %wd does not overlap the input\n", perm, signs, i);

            /* SECONDARY check: acb_theta_all at the reconstructed tau reproduces
               the candidate's squared theta constants projectively on the ten
               even characteristics. */
            if (!_secondary_ok(res, th2, prec))
                TEST_FUNCTION_FAIL("perm %wd signs %wd: acb_theta_all at the "
                    "reconstructed tau does not match th2\n", perm, signs);
        }
    }

    /* The chain yields a valid period matrix at all (guards a vacuous pass). */
    if (n_valid < 1)
        TEST_FUNCTION_FAIL("no (perm, signs) on the grid yielded a certified "
            "period matrix\n");

    acb_poly_clear(f);
    _acb_vec_clear(w, 6);
    _acb_vec_clear(ros, 3);
    _acb_vec_clear(th4, 16);
    _acb_vec_clear(th2, 16);
    _acb_vec_clear(ros_rec, 3);
    acb_mat_clear(res);

    TEST_FUNCTION_END(state);
}
