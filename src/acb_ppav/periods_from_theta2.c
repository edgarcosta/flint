/*
    Copyright (C) 2026 Jean Kieffer
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Ported from HDME theta/theta2_inverse.c (Dupont's PhD thesis, Algorithme 13)
   and theta/theta2_invalid.c (the certain-0 path), both (C) Jean Kieffer.

   Structure. From the 16 squared theta constants th (acb_theta ordering, up to
   a common nonzero scalar) we run four Borchardt (AGM) means m0,...,m3 and
   assemble the reduced period data

       tau11    = i m0 / m1,
       tau22    = i m0 / m2,
       tau12^2  = m0 / m3 + tau11 tau22.

   HDME takes a further Borchardt square root of tau12^2 to get tau12; we stop
   at tau12^2, which is exactly the (0,1) output entry, so that step is omitted.

   Index mapping (the single highest-risk detail). HDME's theta2_inverse selects
   the j-th input of sequence i as th[theta_char_set_label_g2(L[i][j])], where L
   is HDME's label matrix (below) and

       theta_char_set_label_g2(label)      [hdme theta_char_set_label_g2.c]
         = theta_char_set_ab(revbin2(label>>2), revbin2(label & 3), 2)
         = 8 a0 + 4 a1 + 2 b0 + b1,   with label = (a1 a0 b1 b0) in binary.

   That storage index equals acb_theta's k-index: acb_theta packs k = 8 a[0] +
   4 a[1] + 2 b[0] + b[1] with a = ((k>>3)&1,(k>>2)&1), b = ((k>>1)&1,k&1), and
   HDME and acb_theta share the theta_{a,b} definition. Concretely both agree on
     - the naive sum theta_{0,b}(0,tau) = sum q1^{m^2} q2^{n^2} exp(2pi i mn t01)
       (-1)^{m b0 + n b1} (hdme theta_0123_naive.c: th[2] carries the (b0,b1) =
       (1,0) sign pattern, so b0 is the high bit and multiplies the first
       coordinate, matching acb_theta's b[0]); and
     - the duplication theta^2_{a,b}(0,2tau) = 2^{-g} sum_beta (-1)^{<a,beta>}
       theta_{0,beta} theta_{0,b+beta} (hdme theta_duplication.c), which fixes
       the a-convention to agree as well.
   So the acb_theta index of HDME label L is theta_char_set_label_g2(L):

       label : 0 1 2 3 4 6 8 9 12
       k     : 0 2 1 3 8 9 4 6 12

   Applying that map entrywise to HDME's L yields SEQ below. The order within
   each length-4 sequence is significant (the Borchardt step, acb_theta_agm_mul,
   is not symmetric in its inputs), so SEQ preserves HDME's column order.

   HDME label matrix L (for the derivation):
       { {0,1,2,3}, {4,0,6,2}, {8,9,0,1}, {0,8,4,12} }

   Scale freedom. th is defined up to a common nonzero scalar lambda. Each mi is
   a Borchardt mean, homogeneous of degree 1, so mi -> lambda mi; every assembled
   entry is a ratio of means (or a product of two such ratios) and is therefore
   invariant under th -> lambda th. The final recheck compares projectively, so
   it is scale free too.

   Recheck. After assembling tau we recompute the squared theta constants by
   naive summation (acb_theta_sum: independent of both the AGM path here and the
   fast squared algorithm) at low precision and compare projectively with th.
   A certainly-disjoint entry proves th is not of the required form (return 0);
   otherwise, if tau is certifiably reduced and every entry overlaps, return 1;
   anything undecided returns 2. */

#include "arf.h"
#include "arb.h"
#include "arb_mat.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"
#include "acb_ppav.h"

/* acb_theta k-indices of the four Borchardt sequences; see the header comment
   for the derivation from HDME's label matrix. */
static const slong SEQ[4][4] = {
    {0, 2, 1, 3},
    {8, 0, 9, 1},
    {4, 6, 0, 2},
    {0, 4, 8, 12}
};

/* out[k] = theta_{a,b}(0,tau)^2 for the 16 characteristics, computed by naive
   summation. acb_theta_sum returns theta (not its square), so we square; at
   z = 0 the tilde normalization factor exp(-pi y^T Yinv y) = exp(0) = 1, hence
   tilde = 1 needs no unscaling (acb_theta.rst: widetilde{theta} =
   e^{-pi y^T Yinv y} theta). tau must have positive definite imaginary part. */
static void
_theta2_via_sum(acb_ptr out, const acb_mat_t tau, slong prec)
{
    slong g = 2;
    slong n = (WORD(1) << g);
    slong k;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_t ctx_z;
    arb_ptr ds;
    acb_ptr z;

    acb_theta_ctx_tau_init(ctx_tau, 1, g);
    acb_theta_ctx_z_init(ctx_z, g);
    ds = _arb_vec_init(n);
    z = _acb_vec_init(g);

    acb_theta_ctx_tau_set(ctx_tau, tau, prec);
    acb_theta_ctx_z_set(ctx_z, z, ctx_tau, prec);
    acb_theta_eld_distances(ds, z, 1, tau, prec);
    acb_theta_sum(out, ctx_z, 1, ctx_tau, ds, 1, 1, 1, prec);

    for (k = 0; k < n * n; k++)
        acb_sqr(&out[k], &out[k], prec);

    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_clear(ctx_z);
    _arb_vec_clear(ds, n);
    _acb_vec_clear(z, g);
}

/* Index of the entry with the largest certified lower bound on |th[.]|, or -1
   if no entry is certifiably nonzero. Used as the projective normalizer. */
static slong
_pick_nonzero(acb_srcptr th, slong len, slong prec)
{
    slong k, best = -1;
    arf_t lo, hi;

    arf_init(lo);
    arf_init(hi);
    arf_zero(hi);
    for (k = 0; k < len; k++)
    {
        acb_get_abs_lbound_arf(lo, &th[k], prec);
        if (arf_cmp(lo, hi) > 0)
        {
            arf_set(hi, lo);
            best = k;
        }
    }
    arf_clear(lo);
    arf_clear(hi);
    return best;
}

int
acb_ppav_periods_from_theta2(acb_mat_t res, acb_srcptr th, slong g, slong prec)
{
    slong i, j, k, norm, wp, rprec;
    arf_t eps;
    acb_ptr means, seq, rec;
    acb_mat_t tau;
    acb_t w, lhs, rhs;
    arb_mat_t Y, L;
    int ret, code, compat, red, pd;

    FLINT_ASSERT(g == 2);

    wp = prec + 32;
    rprec = 64;

    arf_init(eps);
    arf_set_d(eps, 1e-10);
    means = _acb_vec_init(4);
    seq = _acb_vec_init(4);
    rec = _acb_vec_init(16);
    acb_mat_init(tau, 2, 2);
    acb_init(w);
    acb_init(lhs);
    acb_init(rhs);
    arb_mat_init(Y, 2, 2);
    arb_mat_init(L, 2, 2);

    ret = 2;

    /* Four Borchardt means. acb_ppav_agm's 0 (certainly invalid start) and 2
       (undecidable) propagate unchanged. */
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
            acb_set(&seq[j], &th[SEQ[i][j]]);
        code = acb_ppav_agm(&means[i], seq, eps, 2, wp);
        if (code != 1)
        {
            ret = code;
            goto cleanup;
        }
    }

    /* Assemble r = [[i m0/m1, m0/m3 + tau11 tau22], [., i m0/m2]]. */
    acb_div(acb_mat_entry(res, 0, 0), &means[0], &means[1], wp);
    acb_mul_onei(acb_mat_entry(res, 0, 0), acb_mat_entry(res, 0, 0));
    acb_div(acb_mat_entry(res, 1, 1), &means[0], &means[2], wp);
    acb_mul_onei(acb_mat_entry(res, 1, 1), acb_mat_entry(res, 1, 1));
    acb_div(acb_mat_entry(res, 0, 1), &means[0], &means[3], wp);
    acb_addmul(acb_mat_entry(res, 0, 1), acb_mat_entry(res, 0, 0),
        acb_mat_entry(res, 1, 1), wp);
    acb_set(acb_mat_entry(res, 1, 0), acb_mat_entry(res, 0, 1));

    if (!acb_mat_is_finite(res))
    {
        ret = 2;
        goto cleanup;
    }

    /* tau for the recheck: tau12 is a square root of res01 = tau12^2. theta^2
       is invariant under tau12 -> -tau12 (each theta_{a,b} only picks up
       (-1)^{a1 b1}), so the sign is immaterial to the projective comparison; we
       take Im(tau12) >= 0 so a genuinely reduced input stays reduced. Choosing a
       representative by the midpoint is not a certified 0/1/2 decision. */
    acb_set(acb_mat_entry(tau, 0, 0), acb_mat_entry(res, 0, 0));
    acb_set(acb_mat_entry(tau, 1, 1), acb_mat_entry(res, 1, 1));
    acb_sqrt(w, acb_mat_entry(res, 0, 1), wp);
    if (arf_sgn(arb_midref(acb_imagref(w))) < 0)
        acb_neg(w, w);
    acb_set(acb_mat_entry(tau, 0, 1), w);
    acb_set(acb_mat_entry(tau, 1, 0), w);

    /* Independent low-precision recheck. Only summation over a positive definite
       imaginary part is meaningful, so gate on a certified Cholesky. */
    acb_mat_get_imag(Y, tau);
    pd = arb_mat_cho(L, Y, wp);
    compat = -1;
    if (pd)
    {
        _theta2_via_sum(rec, tau, rprec);
        norm = _pick_nonzero(th, 16, wp);
        if (norm >= 0)
        {
            compat = 1;
            for (k = 0; k < 16; k++)
            {
                /* projective compare: th[k] rec[norm] =?= th[norm] rec[k] */
                acb_mul(lhs, &th[k], &rec[norm], rprec);
                acb_mul(rhs, &th[norm], &rec[k], rprec);
                if (!acb_overlaps(lhs, rhs))
                {
                    compat = 0;
                    break;
                }
            }
        }
    }

    if (compat == 0)
    {
        ret = 0;
        goto cleanup;
    }

    red = acb_siegel_is_reduced(tau, -20, prec);
    if (red && compat == 1)
        ret = 1;
    else
        ret = 2;

cleanup:
    arf_clear(eps);
    _acb_vec_clear(means, 4);
    _acb_vec_clear(seq, 4);
    _acb_vec_clear(rec, 16);
    acb_mat_clear(tau);
    acb_clear(w);
    acb_clear(lhs);
    acb_clear(rhs);
    arb_mat_clear(Y);
    arb_mat_clear(L);
    return ret;
}
