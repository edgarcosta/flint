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

/* Shared deterministic fixture (used by the regression test below and by
   Task 12): a fixed, concrete, asymmetric reduced tau0. Asymmetric means
   tau11 != tau22 and tau12 != 0, so a tau11 <-> tau22 coordinate swap cannot
   hide: such a swap is P tau P^T, which leaves tau12^2 (the (0,1) output)
   unchanged, and equal diagonals would make the swap invisible in the direct
   entrywise checks. The midpoints are the exact doubles below (arb_set_d is
   exact on a dyadic double, radius 0), so acb_contains, not acb_overlaps, is
   the right predicate for the reconstruction. Callers assert
   acb_siegel_is_reduced on this matrix. */
static void
_ppav_test_tau0(acb_mat_t tau)
{
    arb_set_d(acb_realref(acb_mat_entry(tau, 0, 0)), 0.1);
    arb_set_d(acb_imagref(acb_mat_entry(tau, 0, 0)), 1.1);
    arb_set_d(acb_realref(acb_mat_entry(tau, 1, 1)), -0.15);
    arb_set_d(acb_imagref(acb_mat_entry(tau, 1, 1)), 1.4);
    arb_set_d(acb_realref(acb_mat_entry(tau, 0, 1)), 0.2);
    arb_set_d(acb_imagref(acb_mat_entry(tau, 0, 1)), 0.1);
    acb_set(acb_mat_entry(tau, 1, 0), acb_mat_entry(tau, 0, 1));
}

/* Round-trip test. From a random reduced tau, forward-evaluate the 16 squared
   theta constants with acb_theta_all (sqr = 1), invert with
   acb_ppav_periods_from_theta2, and check we recover
   r = [[tau11, tau12^2], [tau12^2, tau22]].

   Three properties are checked on every verified round-trip:

   - Convention: the direct entrywise overlaps res00 ~ tau00, res11 ~ tau11,
     res01 ~ tau01^2. These pin the coordinate convention and catch a
     tau11 <-> tau22 swap or an off-diagonal mismatch that the up-to-Sp(4,Z)
     comparison is blind to (a coordinate swap P leaves P tau P^T with the same
     tau12^2). randtest tau entries carry a small nonzero radius
     (arb_randtest_precise), so these are acb_overlaps, not acb_contains.

   - Equivalence up to Sp(4,Z), via _acb_ppav_periods_overlap_sp4 on a rebuilt
     tb. res01 is only tau12^2, so we pick the square root nearer the input
     tau12 (branch choice done off the input we already know): tb then sits on
     the same reduced representative as tau, which keeps the check stable on
     boundary draws where tau12 is real. The branch-independent convention pin
     is the direct res01 overlap above.

   - Accuracy: the diagonal keeps close to prec relative bits. This only holds
     when the input is well conditioned: at large imaginary part some even
     theta^2 constants are tiny, so at absolute precision 2^{-prec} they lose
     relative accuracy, capping the output. We therefore assert
     acb_rel_accuracy_bits(res_ii) >= prec - 30 only when every even theta^2
     input still carries >= prec - 15 relative bits (the odd characteristics
     vanish at z = 0 and are excluded). Empirically that gate leaves a comfortable
     margin, and it is exercised on a healthy fraction of draws.

   Two facts about the generator shape the loop. acb_siegel_randtest_reduced
   scales the imaginary part AFTER its own looser is_reduced(-1) check with no
   re-check, so a fraction of draws are not -20-reduced and violate the
   contract's precondition (skipped). And at prec = 300 a reduced tau with very
   large imaginary part pushes some theta^2 below the precision floor or onto the
   reduction boundary; the function then correctly returns 2 (insufficient
   precision), which we skip. A reduced input returning 0 is a hard failure, and
   so is failing to accumulate the required number of verified round-trips (this
   guards against an all-2 or all-0 regression). */

TEST_FUNCTION_START(acb_ppav_periods_from_theta2, state)
{
    slong prec = 300;
    slong target = 20;         /* verified round-trips required */
    slong accmin = 3;          /* of which, accuracy-exercised */
    slong maxdraws = 200;
    slong nsucc = 0, nacc = 0, ndraw = 0;
    /* even characteristics (nonzero at z = 0); the AGM uses only these */
    const slong even[10] = {0, 1, 2, 3, 4, 6, 8, 9, 12, 15};

    while ((nsucc < target || nacc < accmin) && ndraw < maxdraws)
    {
        acb_mat_t tau, res, tb;
        acb_ptr th2, z;
        acb_t t12sq;
        slong k, inacc;
        int code;

        ndraw++;
        acb_mat_init(tau, 2, 2);
        acb_mat_init(res, 2, 2);
        acb_mat_init(tb, 2, 2);
        th2 = _acb_vec_init(16);
        z = _acb_vec_init(2);
        acb_init(t12sq);

        acb_siegel_randtest_reduced(tau, state, prec, 4);

        /* precondition: only -20-reduced tau are in scope */
        if (!acb_siegel_is_reduced(tau, -20, prec))
            goto next;

        acb_theta_all(th2, z, tau, 1, prec);
        code = acb_ppav_periods_from_theta2(res, th2, 2, prec);

        if (code == 2)          /* insufficient precision for this tau: skip */
            goto next;
        if (code != 1)          /* 0 on a valid reduced input is a bug */
            TEST_FUNCTION_FAIL("draw %wd: reduced input got code %d\n", ndraw, code);

        /* convention: direct entrywise overlaps */
        if (!acb_overlaps(acb_mat_entry(res, 0, 0), acb_mat_entry(tau, 0, 0)))
            TEST_FUNCTION_FAIL("draw %wd: res00 does not overlap tau00\n", ndraw);
        if (!acb_overlaps(acb_mat_entry(res, 1, 1), acb_mat_entry(tau, 1, 1)))
            TEST_FUNCTION_FAIL("draw %wd: res11 does not overlap tau11\n", ndraw);
        acb_sqr(t12sq, acb_mat_entry(tau, 0, 1), prec);
        if (!acb_overlaps(acb_mat_entry(res, 0, 1), t12sq))
            TEST_FUNCTION_FAIL("draw %wd: res01 does not overlap tau01^2\n", ndraw);

        /* equivalence up to Sp(4,Z); pick the root of res01 nearer input tau12 */
        acb_set(acb_mat_entry(tb, 0, 0), acb_mat_entry(res, 0, 0));
        acb_set(acb_mat_entry(tb, 1, 1), acb_mat_entry(res, 1, 1));
        acb_sqrt(acb_mat_entry(tb, 0, 1), acb_mat_entry(res, 0, 1), prec);
        {
            acb_t d;
            arb_t nd, np;
            acb_init(d);
            arb_init(nd);
            arb_init(np);
            acb_sub(d, acb_mat_entry(tb, 0, 1), acb_mat_entry(tau, 0, 1), prec);
            acb_abs(nd, d, prec);
            acb_add(d, acb_mat_entry(tb, 0, 1), acb_mat_entry(tau, 0, 1), prec);
            acb_abs(np, d, prec);
            if (arb_gt(nd, np))         /* -w is nearer tau12 than w */
                acb_neg(acb_mat_entry(tb, 0, 1), acb_mat_entry(tb, 0, 1));
            acb_clear(d);
            arb_clear(nd);
            arb_clear(np);
        }
        acb_set(acb_mat_entry(tb, 1, 0), acb_mat_entry(tb, 0, 1));
        if (!_acb_ppav_periods_overlap_sp4(tb, tau, prec))
            TEST_FUNCTION_FAIL("draw %wd: rebuilt tau not Sp(4,Z)-equivalent\n",
                ndraw);

        /* accuracy, gated on the input being well conditioned */
        inacc = prec + 100;
        for (k = 0; k < 10; k++)
        {
            slong a = acb_rel_accuracy_bits(&th2[even[k]]);
            if (a < inacc)
                inacc = a;
        }
        if (inacc >= prec - 15)
        {
            if (acb_rel_accuracy_bits(acb_mat_entry(res, 0, 0)) < prec - 30)
                TEST_FUNCTION_FAIL("draw %wd: res00 lost bits (%wd)\n", ndraw,
                    acb_rel_accuracy_bits(acb_mat_entry(res, 0, 0)));
            if (acb_rel_accuracy_bits(acb_mat_entry(res, 1, 1)) < prec - 30)
                TEST_FUNCTION_FAIL("draw %wd: res11 lost bits (%wd)\n", ndraw,
                    acb_rel_accuracy_bits(acb_mat_entry(res, 1, 1)));
            nacc++;
        }

        nsucc++;

next:
        acb_clear(t12sq);
        _acb_vec_clear(th2, 16);
        _acb_vec_clear(z, 2);
        acb_mat_clear(tau);
        acb_mat_clear(res);
        acb_mat_clear(tb);
    }

    if (nsucc < target)
        TEST_FUNCTION_FAIL("only %wd/%wd verified round-trips in %wd draws\n",
            nsucc, target, ndraw);
    if (nacc < accmin)
        TEST_FUNCTION_FAIL("only %wd/%wd accuracy-exercised round-trips\n",
            nacc, accmin);

    TEST_FUNCTION_END(state);
}

/* Deterministic regression on the single fixed asymmetric reduced fixture
   _ppav_test_tau0 (no randomness). Because the fixture midpoints are exact and
   the fixture is reduced and well conditioned by construction, every check is
   unconditional: the inverse must return 1, each output entry must CONTAIN the
   exact input entry (res00 ~ tau11, res11 ~ tau22, res01 ~ tau12^2), the two
   diagonal entries must keep at least prec - 30 relative bits (no gating), and
   the rebuilt matrix must be Sp(4,Z)-equivalent to tau0. The asymmetry is what
   gives the entrywise checks teeth: swapping the two diagonal assertions
   (res00 vs tau22) fails, since tau11 != tau22 and the up-to-Sp(4,Z) check is
   blind to the swap. */
TEST_FUNCTION_START(acb_ppav_periods_from_theta2_regression, state)
{
    slong prec = 300;
    acb_mat_t tau0, res, tb;
    acb_ptr th2, z;
    acb_t t12sq, w;
    int code;

    acb_mat_init(tau0, 2, 2);
    acb_mat_init(res, 2, 2);
    acb_mat_init(tb, 2, 2);
    th2 = _acb_vec_init(16);
    z = _acb_vec_init(2);
    acb_init(t12sq);
    acb_init(w);

    _ppav_test_tau0(tau0);

    /* self-validating fixture: the precondition is that tau0 is -20-reduced */
    if (!acb_siegel_is_reduced(tau0, -20, prec))
        TEST_FUNCTION_FAIL("fixture tau0 is not -20-reduced\n");

    /* forward: the 16 squared theta constants at z = 0 */
    acb_theta_all(th2, z, tau0, 1, prec);

    /* inverse: a reduced, well-conditioned input must return 1 */
    code = acb_ppav_periods_from_theta2(res, th2, 2, prec);
    if (code != 1)
        TEST_FUNCTION_FAIL("periods_from_theta2 returned %d, expected 1\n", code);

    /* direct entrywise containment (tau0 entries are exact points, radius 0) */
    if (!acb_contains(acb_mat_entry(res, 0, 0), acb_mat_entry(tau0, 0, 0)))
        TEST_FUNCTION_FAIL("res00 does not contain tau0_00\n");
    if (!acb_contains(acb_mat_entry(res, 1, 1), acb_mat_entry(tau0, 1, 1)))
        TEST_FUNCTION_FAIL("res11 does not contain tau0_11\n");
    /* tau0_01 is a dyadic point, so its square is exact at this precision */
    acb_sqr(t12sq, acb_mat_entry(tau0, 0, 1), prec);
    if (!acb_contains(acb_mat_entry(res, 0, 1), t12sq))
        TEST_FUNCTION_FAIL("res01 does not contain tau0_01^2\n");

    /* unconditional relative accuracy: the fixture is well conditioned, and
       tau12 != 0 so res01 = tau12^2 has well-defined relative accuracy */
    if (acb_rel_accuracy_bits(acb_mat_entry(res, 0, 0)) < prec - 30)
        TEST_FUNCTION_FAIL("res00 lost bits (%wd)\n",
            acb_rel_accuracy_bits(acb_mat_entry(res, 0, 0)));
    if (acb_rel_accuracy_bits(acb_mat_entry(res, 1, 1)) < prec - 30)
        TEST_FUNCTION_FAIL("res11 lost bits (%wd)\n",
            acb_rel_accuracy_bits(acb_mat_entry(res, 1, 1)));
    if (acb_rel_accuracy_bits(acb_mat_entry(res, 0, 1)) < prec - 30)
        TEST_FUNCTION_FAIL("res01 lost bits (%wd)\n",
            acb_rel_accuracy_bits(acb_mat_entry(res, 0, 1)));

    /* equivalence up to Sp(4,Z) on the rebuilt matrix; tau0_12 has positive
       imaginary part, so the Im >= 0 root of res01 recovers it */
    acb_set(acb_mat_entry(tb, 0, 0), acb_mat_entry(res, 0, 0));
    acb_set(acb_mat_entry(tb, 1, 1), acb_mat_entry(res, 1, 1));
    acb_sqrt(w, acb_mat_entry(res, 0, 1), prec);
    if (arf_sgn(arb_midref(acb_imagref(w))) < 0)
        acb_neg(w, w);
    acb_set(acb_mat_entry(tb, 0, 1), w);
    acb_set(acb_mat_entry(tb, 1, 0), w);
    if (!_acb_ppav_periods_overlap_sp4(tb, tau0, prec))
        TEST_FUNCTION_FAIL("rebuilt tau not Sp(4,Z)-equivalent to tau0\n");

    acb_clear(t12sq);
    acb_clear(w);
    _acb_vec_clear(th2, 16);
    _acb_vec_clear(z, 2);
    acb_mat_clear(tau0);
    acb_mat_clear(res);
    acb_mat_clear(tb);

    TEST_FUNCTION_END(state);
}

/* Scalar invariance: the core contract. th is only defined up to a common
   nonzero scalar (acb_ppav.rst), so periods_from_theta2 must return the same
   reduced period data for th and for s th, any nonzero s. Both forward oracles
   feed theta00^2 near 1, so an impl that silently assumes a normalized first
   entry, or rechecks the recomputed theta^2 against th absolutely instead of
   projectively, passes every other test and fails only this one.

   From each random reduced, code-1 draw we rescale the whole 16-vector by four
   deliberately awkward scalars and require code 1 with output overlapping the
   unscaled res0 entrywise:
     s1 = 3.7 e^{1.3 i}    generic magnitude and argument (rules out a real- or
                           unit-scalar-only normalizer);
     s2 = 10^6             large real (rules out one keyed to magnitude near 1);
     s3 = 2^-40            tiny real (guards one that divides by an entry only
                           when it looks "big");
     s4 = 2^-(prec+100)    far below the target precision. The AGM normalizes the
                           common magnitude by an exact power of two before it
                           derives its step counts, so the after-quad-conv count
                           log2(log2(M0/7)+prec+1) stays well-posed and the
                           residual widening is scale-relative even here. Without
                           that normalization M0 ~ 2^-(prec+100) sends the inner
                           log2 nonpositive and the AGM aborts.
   Overlap is guaranteed, not just likely: periods_from_theta2 returns a rigorous
   enclosure of the same mathematical r in every case, so res0 and each res both
   contain that r. The scalars are exact (s2, s3, s4) or a precise nonzero ball
   (s1), and none contains zero, so rescaling preserves the certified nonzero
   magnitudes the AGM and its projective recheck rely on; scaling by an exact
   power of two keeps the zero/nonzero status of every entry byte for byte. The
   draw-and-skip discipline matches the round-trip above (skip non -20-reduced
   draws and honest 2-returns; a reduced input returning 0 is a hard failure), and
   we require a few verified draws so the check is exercised on more than one tau
   shape.

   On every scaled call we also require each assembled entry to keep at least
   prec - 30 relative bits. The entries are scale-free ratios of degree-1
   Borchardt means, so an exact rescaling of the input cannot degrade them; a
   scale-dependent accuracy loss (which overlap alone does not see) fails here.
   res01 = tau12^2 is nonzero on these reduced draws, so its relative accuracy is
   well-defined. */
TEST_FUNCTION_START(acb_ppav_periods_from_theta2_scalar, state)
{
    slong prec = 300;
    slong target = 3;          /* verified scalar-invariant draws required */
    slong maxdraws = 200;
    slong nsucc = 0, ndraw = 0;
    acb_ptr s = _acb_vec_init(3);

    /* s1 = 3.7 e^{1.3 i} = 3.7 (cos 1.3 + i sin 1.3), built via arb trig */
    {
        arb_t ang, c, sn, mag;
        arb_init(ang);
        arb_init(c);
        arb_init(sn);
        arb_init(mag);
        arb_set_d(ang, 1.3);
        arb_sin_cos(sn, c, ang, prec);
        arb_set_d(mag, 3.7);
        arb_mul(acb_realref(&s[0]), mag, c, prec);
        arb_mul(acb_imagref(&s[0]), mag, sn, prec);
        arb_clear(ang);
        arb_clear(c);
        arb_clear(sn);
        arb_clear(mag);
    }
    acb_set_ui(&s[1], 1000000);        /* s2 = 10^6, exact */
    acb_one(&s[2]);                    /* s3 = 2^-40, exact */
    acb_mul_2exp_si(&s[2], &s[2], -40);

    while (nsucc < target && ndraw < maxdraws)
    {
        acb_mat_t tau, res0, res;
        acb_ptr th2, th2s, z;
        slong si, k;
        int code;

        ndraw++;
        acb_mat_init(tau, 2, 2);
        acb_mat_init(res0, 2, 2);
        acb_mat_init(res, 2, 2);
        th2 = _acb_vec_init(16);
        th2s = _acb_vec_init(16);
        z = _acb_vec_init(2);

        acb_siegel_randtest_reduced(tau, state, prec, 4);

        /* precondition: only -20-reduced tau are in scope */
        if (!acb_siegel_is_reduced(tau, -20, prec))
            goto next;

        acb_theta_all(th2, z, tau, 1, prec);
        code = acb_ppav_periods_from_theta2(res0, th2, 2, prec);

        if (code == 2)          /* insufficient precision for this tau: skip */
            goto next;
        if (code != 1)          /* 0 on a valid reduced input is a bug */
            TEST_FUNCTION_FAIL("draw %wd: unscaled reduced input got code %d\n",
                ndraw, code);

        for (si = 0; si < 4; si++)
        {
            if (si < 3)
                _acb_vec_scalar_mul(th2s, th2, 16, &s[si], prec);
            else
                /* s4 = 2^-(prec+100), applied exactly per entry */
                _acb_vec_scalar_mul_2exp_si(th2s, th2, 16, -(prec + 100));
            code = acb_ppav_periods_from_theta2(res, th2s, 2, prec);
            if (code != 1)
                TEST_FUNCTION_FAIL("draw %wd scalar %wd: got code %d, want 1\n",
                    ndraw, si, code);
            for (k = 0; k < 4; k++)
                if (!acb_overlaps(acb_mat_entry(res, k / 2, k % 2),
                        acb_mat_entry(res0, k / 2, k % 2)))
                    TEST_FUNCTION_FAIL("draw %wd scalar %wd: entry (%wd,%wd) "
                        "does not overlap res0\n", ndraw, si, k / 2, k % 2);
            /* scale-free ratios: an exact rescaling must not cost accuracy */
            if (acb_rel_accuracy_bits(acb_mat_entry(res, 0, 0)) < prec - 30)
                TEST_FUNCTION_FAIL("draw %wd scalar %wd: res00 lost bits (%wd)\n",
                    ndraw, si, acb_rel_accuracy_bits(acb_mat_entry(res, 0, 0)));
            if (acb_rel_accuracy_bits(acb_mat_entry(res, 1, 1)) < prec - 30)
                TEST_FUNCTION_FAIL("draw %wd scalar %wd: res11 lost bits (%wd)\n",
                    ndraw, si, acb_rel_accuracy_bits(acb_mat_entry(res, 1, 1)));
            if (acb_rel_accuracy_bits(acb_mat_entry(res, 0, 1)) < prec - 30)
                TEST_FUNCTION_FAIL("draw %wd scalar %wd: res01 lost bits (%wd)\n",
                    ndraw, si, acb_rel_accuracy_bits(acb_mat_entry(res, 0, 1)));
        }

        nsucc++;

next:
        _acb_vec_clear(th2, 16);
        _acb_vec_clear(th2s, 16);
        _acb_vec_clear(z, 2);
        acb_mat_clear(tau);
        acb_mat_clear(res0);
        acb_mat_clear(res);
    }

    if (nsucc < target)
        TEST_FUNCTION_FAIL("only %wd/%wd verified scalar-invariant draws in "
            "%wd draws\n", nsucc, target, ndraw);

    _acb_vec_clear(s, 3);

    TEST_FUNCTION_END(state);
}

/* Certain-invalid negative test: periods_from_theta2 must return 0 (certainly
   not a valid theta^2 vector) on decidably invalid input, never 1 and never a
   hedging 2. Three inputs, all certainly invalid at 300 bits:

   (1) A valid theta^2 vector for the fixture tau0 with one ODD entry forced
       certainly nonzero. The six odd characteristics k in {5,7,10,11,13,14}
       have theta_{a,b}(0,tau) identically zero, so th2[5] = 0.3 cannot arise
       from any tau. The four Borchardt sequences (SEQ) consume only even
       characteristics, so the assembled tau is still the correct tau0 and only
       the recheck can catch the planted entry: a code-1 return would mean the
       recheck skips odd positions (a production gap), a code-2 return would mean
       it hedges on a decidable mismatch. Both are failures by design.

   (2) The coarse all-ones vector: the odd positions are then certainly nonzero
       too, so again certainly invalid.

   (3) A certain invalidity that lives in a LATER Borchardt sequence than the
       first undecidable one, exercising the "run all four sequences" path (HDME
       theta2_invalid). SEQ0 = {th0,th2,th1,th3} = {1, [0+/-1/2], 1, 1} is
       undecidable (the fat zero kills the certain-yes half-plane and the three
       1's do not cover), so its AGM returns 2. SEQ1 = {th8,th0,th9,th1} =
       {-1, 1, 1, 1} certainly leaves every half-plane, so its AGM returns 0.
       Stopping at the first non-1 sequence would return 2; because a later
       sequence certifies invalidity, the answer must be 0. */
TEST_FUNCTION_START(acb_ppav_periods_from_theta2_invalid, state)
{
    slong prec = 300;
    acb_mat_t tau0, res;
    acb_ptr th2, z;
    slong k;
    int code;

    acb_mat_init(tau0, 2, 2);
    acb_mat_init(res, 2, 2);
    th2 = _acb_vec_init(16);
    z = _acb_vec_init(2);

    _ppav_test_tau0(tau0);

    /* (1) valid even entries, one certainly-nonzero odd entry */
    acb_theta_all(th2, z, tau0, 1, prec);
    acb_set_d(&th2[5], 0.3);
    code = acb_ppav_periods_from_theta2(res, th2, 2, prec);
    if (code != 0)
        TEST_FUNCTION_FAIL("odd entry th2[5]=0.3: got code %d, expected 0\n", code);

    /* (2) coarsely invalid all-ones */
    for (k = 0; k < 16; k++)
        acb_one(&th2[k]);
    code = acb_ppav_periods_from_theta2(res, th2, 2, prec);
    if (code != 0)
        TEST_FUNCTION_FAIL("all-ones th2: got code %d, expected 0\n", code);

    /* (3) invalidity in SEQ1 while SEQ0 is only undecidable (2). SEQ0 =
       {th0,th2,th1,th3} = {1, fat-zero, 1, 1} -> AGM 2; SEQ1 = {th8,th0,th9,th1}
       = {-1, 1, 1, 1} -> AGM 0. Only reachable by continuing past SEQ0's 2. */
    for (k = 0; k < 16; k++)
        acb_zero(&th2[k]);
    acb_one(&th2[0]);
    acb_one(&th2[1]);
    acb_one(&th2[3]);
    arb_add_error_2exp_si(acb_realref(&th2[2]), -1);   /* th2[2] = [0 +/- 1/2] */
    acb_set_si(&th2[8], -1);
    acb_one(&th2[9]);
    code = acb_ppav_periods_from_theta2(res, th2, 2, prec);
    if (code != 0)
        TEST_FUNCTION_FAIL("late-sequence invalid: got code %d, expected 0\n", code);

    _acb_vec_clear(th2, 16);
    _acb_vec_clear(z, 2);
    acb_mat_clear(tau0);
    acb_mat_clear(res);

    TEST_FUNCTION_END(state);
}

/* Near-boundary escalation test. acb_siegel_randtest_compact draws reduced tau
   closer to the reduction boundary than randtest_reduced (bounded imaginary
   part, exact entries), exercising the half-plane and square-root sign choices
   that deep-interior draws never reach. For each kept tau we escalate the working
   precision through {300, 600, 1200}: code 0 is impossible on a valid input, code
   2 is an honest "insufficient precision" and escalates, and code 1 must be sound
   (direct entrywise overlaps against tau plus Sp(4,Z)-equivalence, as in the
   round-trip). An honest 2 at low precision is allowed, but a perpetual 2 (still
   2 after 1200 bits) is a failure: the contract must reach a certified 1 on these
   sign-machinery-exercising inputs. */
TEST_FUNCTION_START(acb_ppav_periods_from_theta2_boundary, state)
{
    slong genprec = 300;
    slong target = 4;          /* completed (eventually code-1) tau required */
    slong maxdraws = 100;
    slong nsucc = 0, ndraw = 0;
    const slong precs[3] = {300, 600, 1200};

    /* Deterministic 2 -> 1 witness on a single fixed input. tau = diag(96i, 97i)
       with tau12 = 2^-4 is -20-reduced but has even theta^2 constants small
       enough (~exp(-96 pi)) that at 300 bits the AGM cannot certify convergence
       and periods_from_theta2 returns 2; by 600 bits it succeeds with 1. All
       entries are exact dyadics, so both codes are deterministic. This positively
       exercises the caller-side 2 -> 1 escalation on one input, which the random
       ladder below only drives to eventual success in aggregate. Empirically the
       transition sits near 440 bits, so 300 (solid 2) and 600 (solid 1) both have
       margin. */
    {
        acb_mat_t wtau, wres;
        acb_ptr wth2, wz;
        acb_t wt12sq;
        int wcode;

        acb_mat_init(wtau, 2, 2);
        acb_mat_init(wres, 2, 2);
        wth2 = _acb_vec_init(16);
        wz = _acb_vec_init(2);
        acb_init(wt12sq);

        arb_set_si(acb_imagref(acb_mat_entry(wtau, 0, 0)), 96);
        arb_set_si(acb_imagref(acb_mat_entry(wtau, 1, 1)), 97);
        arb_set_d(acb_realref(acb_mat_entry(wtau, 0, 1)), 0.0625);   /* 2^-4 */
        acb_set(acb_mat_entry(wtau, 1, 0), acb_mat_entry(wtau, 0, 1));

        if (!acb_siegel_is_reduced(wtau, -20, genprec))
            TEST_FUNCTION_FAIL("2->1 witness tau is not -20-reduced\n");

        /* 300 bits: honest insufficient precision */
        acb_theta_all(wth2, wz, wtau, 1, 300);
        wcode = acb_ppav_periods_from_theta2(wres, wth2, 2, 300);
        if (wcode != 2)
            TEST_FUNCTION_FAIL("2->1 witness: expected code 2 at 300, got %d\n", wcode);

        /* 600 bits: certified 1, and sound in the res encoding */
        acb_theta_all(wth2, wz, wtau, 1, 600);
        wcode = acb_ppav_periods_from_theta2(wres, wth2, 2, 600);
        if (wcode != 1)
            TEST_FUNCTION_FAIL("2->1 witness: expected code 1 at 600, got %d\n", wcode);
        if (!acb_overlaps(acb_mat_entry(wres, 0, 0), acb_mat_entry(wtau, 0, 0)))
            TEST_FUNCTION_FAIL("2->1 witness: res00 does not overlap tau00\n");
        if (!acb_overlaps(acb_mat_entry(wres, 1, 1), acb_mat_entry(wtau, 1, 1)))
            TEST_FUNCTION_FAIL("2->1 witness: res11 does not overlap tau11\n");
        acb_sqr(wt12sq, acb_mat_entry(wtau, 0, 1), 600);
        if (!acb_overlaps(acb_mat_entry(wres, 0, 1), wt12sq))
            TEST_FUNCTION_FAIL("2->1 witness: res01 does not overlap tau12^2\n");

        acb_clear(wt12sq);
        _acb_vec_clear(wth2, 16);
        _acb_vec_clear(wz, 2);
        acb_mat_clear(wtau);
        acb_mat_clear(wres);
    }

    while (nsucc < target && ndraw < maxdraws)
    {
        acb_mat_t tau;
        slong pi;
        int resolved = 0;

        ndraw++;
        acb_mat_init(tau, 2, 2);
        acb_siegel_randtest_compact(tau, state, 1, genprec);

        /* precondition: only -20-reduced tau are in scope. tau is exact, so
           certifying reducedness once at genprec suffices for higher precs. */
        if (!acb_siegel_is_reduced(tau, -20, genprec))
        {
            acb_mat_clear(tau);
            continue;
        }

        for (pi = 0; pi < 3 && !resolved; pi++)
        {
            slong prec = precs[pi];
            acb_mat_t res, tb;
            acb_ptr th2, z;
            acb_t t12sq;
            int code;

            acb_mat_init(res, 2, 2);
            acb_mat_init(tb, 2, 2);
            th2 = _acb_vec_init(16);
            z = _acb_vec_init(2);
            acb_init(t12sq);

            acb_theta_all(th2, z, tau, 1, prec);
            code = acb_ppav_periods_from_theta2(res, th2, 2, prec);

            if (code == 0)
                TEST_FUNCTION_FAIL("draw %wd prec %wd: valid input got code 0\n",
                    ndraw, prec);
            if (code == 1)
            {
                /* convention: direct entrywise overlaps against the input tau */
                if (!acb_overlaps(acb_mat_entry(res, 0, 0), acb_mat_entry(tau, 0, 0)))
                    TEST_FUNCTION_FAIL("draw %wd: res00 does not overlap tau00\n", ndraw);
                if (!acb_overlaps(acb_mat_entry(res, 1, 1), acb_mat_entry(tau, 1, 1)))
                    TEST_FUNCTION_FAIL("draw %wd: res11 does not overlap tau11\n", ndraw);
                acb_sqr(t12sq, acb_mat_entry(tau, 0, 1), prec);
                if (!acb_overlaps(acb_mat_entry(res, 0, 1), t12sq))
                    TEST_FUNCTION_FAIL("draw %wd: res01 does not overlap tau01^2\n", ndraw);

                /* equivalence up to Sp(4,Z); pick root of res01 nearer input tau12 */
                acb_set(acb_mat_entry(tb, 0, 0), acb_mat_entry(res, 0, 0));
                acb_set(acb_mat_entry(tb, 1, 1), acb_mat_entry(res, 1, 1));
                acb_sqrt(acb_mat_entry(tb, 0, 1), acb_mat_entry(res, 0, 1), prec);
                {
                    acb_t d;
                    arb_t nd, np;
                    acb_init(d);
                    arb_init(nd);
                    arb_init(np);
                    acb_sub(d, acb_mat_entry(tb, 0, 1), acb_mat_entry(tau, 0, 1), prec);
                    acb_abs(nd, d, prec);
                    acb_add(d, acb_mat_entry(tb, 0, 1), acb_mat_entry(tau, 0, 1), prec);
                    acb_abs(np, d, prec);
                    if (arb_gt(nd, np))
                        acb_neg(acb_mat_entry(tb, 0, 1), acb_mat_entry(tb, 0, 1));
                    acb_clear(d);
                    arb_clear(nd);
                    arb_clear(np);
                }
                acb_set(acb_mat_entry(tb, 1, 0), acb_mat_entry(tb, 0, 1));
                if (!_acb_ppav_periods_overlap_sp4(tb, tau, prec))
                    TEST_FUNCTION_FAIL("draw %wd: rebuilt tau not Sp(4,Z)-equivalent\n",
                        ndraw);

                resolved = 1;
            }
            /* code == 2: escalate to the next precision */

            acb_clear(t12sq);
            _acb_vec_clear(th2, 16);
            _acb_vec_clear(z, 2);
            acb_mat_clear(res);
            acb_mat_clear(tb);
        }

        if (!resolved)
            TEST_FUNCTION_FAIL("draw %wd: still code 2 after 1200 bits\n", ndraw);

        nsucc++;
        acb_mat_clear(tau);
    }

    if (nsucc < target)
        TEST_FUNCTION_FAIL("only %wd/%wd boundary cases resolved in %wd draws\n",
            nsucc, target, ndraw);

    TEST_FUNCTION_END(state);
}

/* Dump the recorded per-precision escalation outcomes. Called only on a failure
   path (a passing run stays silent per FLINT convention); rel accuracy is shown
   only for a code-1 rung, where res is fully assembled. codes[i] < 0 marks a
   rung the ladder did not reach before failing. */
static void
_escalation_report(const slong * precs, const int * codes,
    const slong * acc00, const slong * acc11, slong n)
{
    slong i;
    flint_printf("escalation outcomes:");
    for (i = 0; i < n; i++)
    {
        if (codes[i] < 0)
            flint_printf(" %wd:-", precs[i]);
        else if (codes[i] == 1)
            flint_printf(" %wd:1(acc %wd,%wd)", precs[i], acc00[i], acc11[i]);
        else
            flint_printf(" %wd:%d", precs[i], codes[i]);
    }
    flint_printf("\n");
}

/* Precision-escalation soundness and convergence on the shared deterministic
   fixture _ppav_test_tau0. The fixture is exact and reduced, so it is a valid
   theta^2 input at every precision; the return contract is exercised as the
   working precision climbs the ladder {24, 48, 96, 192, 384}.

   The whole pipeline runs at each precision: th2 is recomputed with
   acb_theta_all at that precision, since feeding a fixed high-precision th2 to a
   low-precision inverse would test a different question. At each rung:

     - code 0 is always a hard failure: the input is valid at every precision, so
       0 (certainly invalid) is never correct here.
     - code 1 must be SOUND in the res encoding (see _ppav_test_tau0): res00
       contains tau0_00, res11 contains tau0_11, res01 contains tau0_01^2. The
       square is exact at 512 bits because tau0_01 is a dyadic point, so the
       target carries radius 0 and acb_contains is the right predicate. A code-1
       ball that misses the truth is a wrong tight answer, the worst outcome, and
       must fail loudly.
     - code 2 is an honest "insufficient precision" and is acceptable at low
       precision; it is recorded and the ladder continues.

   Two end conditions. EVENTUAL SUCCESS: the top precision (384) must return 1,
   so a perpetual-2 implementation fails this test. CONVERGENCE: at the top
   precision the two diagonal entries and the off-diagonal res01 = tau12^2 each
   keep at least prec - 30 relative bits (tau0_12 != 0, so res01 is nonzero).
   No per-step radius monotonicity is asserted across the ladder: step count,
   beta, and the internal ellipsoid can change with precision, so a wider ball at
   a higher precision is not a bug, and only the top-precision accuracy is
   pinned. Per-precision outcomes are recorded and printed only on a failure
   path. */
TEST_FUNCTION_START(acb_ppav_periods_from_theta2_escalation, state)
{
    const slong nprecs = 5;
    const slong precs[5] = {24, 48, 96, 192, 384};
    const slong toprec = 384;
    int codes[5];
    slong acc00[5], acc11[5], acc01[5];
    slong pi;
    acb_mat_t tau0;
    acb_t truth01;

    acb_mat_init(tau0, 2, 2);
    acb_init(truth01);

    for (pi = 0; pi < nprecs; pi++)
    {
        codes[pi] = -1;
        acc00[pi] = 0;
        acc11[pi] = 0;
        acc01[pi] = 0;
    }

    _ppav_test_tau0(tau0);

    /* precondition: the fixture is -20-reduced. tau0 is exact, so certifying it
       once at high precision covers every rung of the ladder. */
    if (!acb_siegel_is_reduced(tau0, -20, 512))
        TEST_FUNCTION_FAIL("fixture tau0 is not -20-reduced\n");

    /* res01 target: tau0_01 is a dyadic point, so its square is exact at 512 bits */
    acb_sqr(truth01, acb_mat_entry(tau0, 0, 1), 512);

    for (pi = 0; pi < nprecs; pi++)
    {
        slong prec = precs[pi];
        acb_mat_t res;
        acb_ptr th2, z;
        int code;

        acb_mat_init(res, 2, 2);
        th2 = _acb_vec_init(16);
        z = _acb_vec_init(2);

        /* run the whole pipeline at this precision */
        acb_theta_all(th2, z, tau0, 1, prec);
        code = acb_ppav_periods_from_theta2(res, th2, 2, prec);

        codes[pi] = code;
        acc00[pi] = acb_rel_accuracy_bits(acb_mat_entry(res, 0, 0));
        acc11[pi] = acb_rel_accuracy_bits(acb_mat_entry(res, 1, 1));
        acc01[pi] = acb_rel_accuracy_bits(acb_mat_entry(res, 0, 1));

        /* the fixture is valid, so a certainly-invalid verdict is a bug */
        if (code == 0)
        {
            _escalation_report(precs, codes, acc00, acc11, nprecs);
            TEST_FUNCTION_FAIL("prec %wd: valid fixture got code 0\n", prec);
        }

        /* soundness: a code-1 ball must contain the truth in the res encoding */
        if (code == 1)
        {
            if (!acb_contains(acb_mat_entry(res, 0, 0), acb_mat_entry(tau0, 0, 0)))
            {
                _escalation_report(precs, codes, acc00, acc11, nprecs);
                TEST_FUNCTION_FAIL("prec %wd: res00 does not contain tau0_00\n", prec);
            }
            if (!acb_contains(acb_mat_entry(res, 1, 1), acb_mat_entry(tau0, 1, 1)))
            {
                _escalation_report(precs, codes, acc00, acc11, nprecs);
                TEST_FUNCTION_FAIL("prec %wd: res11 does not contain tau0_11\n", prec);
            }
            if (!acb_contains(acb_mat_entry(res, 0, 1), truth01))
            {
                _escalation_report(precs, codes, acc00, acc11, nprecs);
                TEST_FUNCTION_FAIL("prec %wd: res01 does not contain tau0_01^2\n", prec);
            }
        }
        /* code == 2: honest insufficient precision, escalate */

        _acb_vec_clear(th2, 16);
        _acb_vec_clear(z, 2);
        acb_mat_clear(res);
    }

    /* eventual success: the ladder must terminate in a certified 1 at the top */
    if (codes[nprecs - 1] != 1)
    {
        _escalation_report(precs, codes, acc00, acc11, nprecs);
        TEST_FUNCTION_FAIL("top prec %wd: got code %d, expected eventual 1\n",
            toprec, codes[nprecs - 1]);
    }

    /* convergence: the top-precision diagonal keeps at least prec - 30 rel bits */
    if (acc00[nprecs - 1] < toprec - 30)
    {
        _escalation_report(precs, codes, acc00, acc11, nprecs);
        TEST_FUNCTION_FAIL("top prec %wd: res00 rel accuracy %wd < %wd\n",
            toprec, acc00[nprecs - 1], toprec - 30);
    }
    if (acc11[nprecs - 1] < toprec - 30)
    {
        _escalation_report(precs, codes, acc00, acc11, nprecs);
        TEST_FUNCTION_FAIL("top prec %wd: res11 rel accuracy %wd < %wd\n",
            toprec, acc11[nprecs - 1], toprec - 30);
    }
    if (acc01[nprecs - 1] < toprec - 30)
    {
        _escalation_report(precs, codes, acc00, acc11, nprecs);
        TEST_FUNCTION_FAIL("top prec %wd: res01 rel accuracy %wd < %wd\n",
            toprec, acc01[nprecs - 1], toprec - 30);
    }

    acb_clear(truth01);
    acb_mat_clear(tau0);

    TEST_FUNCTION_END(state);
}
