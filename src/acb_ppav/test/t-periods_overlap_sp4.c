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
#include "fmpz_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"
#include "acb_ppav.h"

/* Exercise the certified up-to-Sp(4,Z) comparison of two 2x2 period matrices.

   The three cases are checked without aborting on the first failure so that a
   single run reports the status of every case; this is what makes the
   identity-witness teeth check legible (removing the empty-product witness must
   knock out case (C) specifically). */

TEST_FUNCTION_START(acb_ppav_periods_overlap_sp4, state)
{
    slong prec = 256;
    slong j0 = 0;                        /* the J matrix in sp2gz_fundamental */
    acb_mat_t tau1, tau2, tau3, tauc;
    fmpz_mat_t g0;
    int ra, rb, rc, fail = 0;

    acb_mat_init(tau1, 2, 2);
    acb_mat_init(tau2, 2, 2);
    acb_mat_init(tau3, 2, 2);
    acb_mat_init(tauc, 2, 2);
    fmpz_mat_init(g0, 4, 4);

    acb_siegel_randtest_reduced(tau1, state, prec, 4);

    /* (A) tau2 = g0 . tau1 for a known fundamental g0: the two are Sp(4,Z)
       equivalent by construction, so the comparison must certify (return 1). */
    sp2gz_fundamental(g0, j0);
    acb_siegel_transform(tau2, g0, tau1, prec);
    ra = _acb_ppav_periods_overlap_sp4(tau1, tau2, prec);
    if (ra != 1)
    {
        flint_printf("(A) FAIL: g0.tau1 not certified equivalent (got %d)\n", ra);
        fail = 1;
    }

    /* (B) tau3 = tau1 with 0.3 added to the real part of the (0,1) and (1,0)
       entries (kept symmetric): a genuinely different reduced point that is not
       an Sp(4,Z) neighbour of tau1, so the comparison must return 0. */
    acb_mat_set(tau3, tau1);
    {
        arb_t d;
        arb_init(d);
        arb_set_d(d, 0.3);
        arb_add(acb_realref(acb_mat_entry(tau3, 0, 1)),
                acb_realref(acb_mat_entry(tau3, 0, 1)), d, prec);
        arb_add(acb_realref(acb_mat_entry(tau3, 1, 0)),
                acb_realref(acb_mat_entry(tau3, 1, 0)), d, prec);
        arb_clear(d);
    }
    rb = _acb_ppav_periods_overlap_sp4(tau1, tau3, prec);
    if (rb != 0)
    {
        flint_printf("(B) FAIL: perturbed tau3 certified equivalent (got %d)\n", rb);
        fail = 1;
    }

    /* (C) BLOCKER regression: tau2 = tau1 exactly. Both reduce to the same
       representative, witnessed by the identity (empty product); no non-trivial
       fundamental matrix fixes a reduced tau. This fails on the un-corrected
       enumeration (19 fundamental only) and passes once the identity witness is
       included. */
    acb_mat_set(tauc, tau1);
    rc = _acb_ppav_periods_overlap_sp4(tau1, tauc, prec);
    if (rc != 1)
    {
        flint_printf("(C) FAIL: identical tau not certified equivalent (got %d)\n", rc);
        fail = 1;
    }

    if (fail)
        TEST_FUNCTION_FAIL("periods_overlap_sp4: see per-case diagnostics above\n");

    acb_mat_clear(tau1);
    acb_mat_clear(tau2);
    acb_mat_clear(tau3);
    acb_mat_clear(tauc);
    fmpz_mat_clear(g0);

    TEST_FUNCTION_END(state);
}
