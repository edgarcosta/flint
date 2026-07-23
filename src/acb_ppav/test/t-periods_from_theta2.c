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

    /* unconditional relative accuracy: the fixture is well conditioned */
    if (acb_rel_accuracy_bits(acb_mat_entry(res, 0, 0)) < prec - 30)
        TEST_FUNCTION_FAIL("res00 lost bits (%wd)\n",
            acb_rel_accuracy_bits(acb_mat_entry(res, 0, 0)));
    if (acb_rel_accuracy_bits(acb_mat_entry(res, 1, 1)) < prec - 30)
        TEST_FUNCTION_FAIL("res11 lost bits (%wd)\n",
            acb_rel_accuracy_bits(acb_mat_entry(res, 1, 1)));

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
