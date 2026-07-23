/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb.h"
#include "fmpz_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"
#include "acb_ppav.h"

/* Sanity filter on a candidate reduced representative r: r must be symmetric
   and have positive-definite imaginary part, i.e. r must look like a genuine
   period matrix. Positive-definiteness of the 2x2 imaginary part Y is certified
   by Sylvester's criterion, requiring certified-positive lower bounds on the
   two leading principal minors Y[0][0] and det(Y). This is deliberately
   redundant with acb_siegel_is_reduced: it is the belt-and-braces guard that a
   witnessed overlap is being read off a valid point, not numerical debris. */
static int
acb_ppav_mat_is_period_like(const acb_mat_t r, slong prec)
{
    arb_t y00, y11, y01, y10, det;
    int res;

    /* symmetric: (0,1) and (1,0) must overlap (real and imaginary parts) */
    if (!acb_overlaps(acb_mat_entry(r, 0, 1), acb_mat_entry(r, 1, 0)))
        return 0;

    arb_init(y00);
    arb_init(y11);
    arb_init(y01);
    arb_init(y10);
    arb_init(det);

    arb_set(y00, acb_imagref(acb_mat_entry(r, 0, 0)));
    arb_set(y11, acb_imagref(acb_mat_entry(r, 1, 1)));
    arb_set(y01, acb_imagref(acb_mat_entry(r, 0, 1)));
    arb_set(y10, acb_imagref(acb_mat_entry(r, 1, 0)));

    /* leading minors of Y: y00 and det(Y) = y00*y11 - y01*y10 */
    arb_mul(det, y00, y11, prec);
    arb_submul(det, y01, y10, prec);

    res = arb_is_positive(y00) && arb_is_positive(det);

    arb_clear(y00);
    arb_clear(y11);
    arb_clear(y01);
    arb_clear(y10);
    arb_clear(det);
    return res;
}

/* Return 1 iff tau1 and tau2 are certified equivalent under Sp(4,Z).

   Both inputs are reduced with acb_siegel_reduce and mapped to reduced
   representatives r1, r2. If either fails to reduce (acb_siegel_reduce silently
   returns the identity transform near the boundary, so an unreduced output is
   meaningless), we cannot certify anything and return 0. Otherwise we enumerate
   the candidate identifications of the reduced domain: the identity (empty
   product), which witnesses the generic interior case where both inputs reduce
   to the SAME representative, followed by the fundamental neighbours from
   sp2gz_fundamental. A candidate g is accepted when g.r1 overlaps r2 entrywise
   (overlap, never containment: r1 and r2 carry comparable radii) and r2 is
   period-like.

   The certificate is "consistent with Sp(4,Z)-equivalence": corner cases whose
   two reduced representatives are related only by a PRODUCT of fundamental
   matrices are out of scope and reported as 0. */
int
_acb_ppav_periods_overlap_sp4(const acb_mat_t tau1, const acb_mat_t tau2,
    slong prec)
{
    slong g = acb_mat_nrows(tau1);
    fmpz_mat_t m, gg;
    acb_mat_t r1, r2, w;
    slong j;
    int res = 0;

    fmpz_mat_init(m, 2 * g, 2 * g);
    fmpz_mat_init(gg, 2 * g, 2 * g);
    acb_mat_init(r1, g, g);
    acb_mat_init(r2, g, g);
    acb_mat_init(w, g, g);

    acb_siegel_reduce(m, tau1, prec);
    acb_siegel_transform(r1, m, tau1, prec);
    acb_siegel_reduce(m, tau2, prec);
    acb_siegel_transform(r2, m, tau2, prec);

    if (!acb_siegel_is_reduced(r1, -20, prec)
        || !acb_siegel_is_reduced(r2, -20, prec)
        || !acb_ppav_mat_is_period_like(r2, prec))
    {
        goto cleanup;
    }

    /* identity witness first: the generic interior case */
    if (acb_mat_overlaps(r1, r2))
    {
        res = 1;
        goto cleanup;
    }

    for (j = 0; j < sp2gz_nb_fundamental(g); j++)
    {
        sp2gz_fundamental(gg, j);
        acb_siegel_transform(w, gg, r1, prec);
        if (acb_mat_overlaps(w, r2))
        {
            res = 1;
            goto cleanup;
        }
    }

cleanup:
    fmpz_mat_clear(m);
    fmpz_mat_clear(gg);
    acb_mat_clear(r1);
    acb_mat_clear(r2);
    acb_mat_clear(w);
    return res;
}
