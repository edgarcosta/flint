/*
    Copyright (C) 2026 Jean Kieffer
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Ported from HDME theta/borchardt_mean.c and its helpers (all (C) Jean
   Kieffer): borchardt_step.c, borchardt_mean_m0.c, borchardt_mean_M0_.c,
   borchardt_mean_Delta0.c, borchardt_mean_nb_steps_before_quad_conv.c,
   borchardt_mean_nb_steps_after_quad_conv.c, borchardt_mean_quad_conv_is_reached.c.

   Two substantive adaptations to match the acb_ppav.rst contract:

   1. The good half plane is chosen up front by acb_ppav_agm_half_plane, which
      returns an angle beta such that every rescaled value exp(-2 pi i beta) a_j
      lies in the right half plane. HDME instead divides by a_0 (assuming 1 is
      one of the values) and rotates by the average argument inside
      borchardt_mean_m0. Because our vector is already aligned by beta, m0 is
      just min_j lbound(Re(b_j)); no per-call argument rotation is needed. This
      lower bound is <= HDME's rotated m0, so the derived step cap is only ever
      more conservative.

   2. borchardt_step is replaced by acb_theta_agm_mul (all = 0) followed by a
      division by 2^g, which reproduces the classical ((x+y)/2, sqrt(xy)) step
      including the 2^g normalization (acb_theta.rst duplication formula).

   Ported constants, with provenance:
     - quad-convergence threshold delta < 1/7
       (borchardt_mean_quad_conv_is_reached.c: |a_i - a_0| < |a_0|/7);
     - before-quad-conv step cap
       ceil( log(m0/(7 Delta0)) / log(1 - 2^{-g}) )
       (borchardt_mean_nb_steps_before_quad_conv.c; the 2^{-g} generalizes the
       g = 2 code's 2^{-2});
     - after-quad-conv step count ceil( log2( log2(M0/7) + prec + 1 ) )
       (borchardt_mean_nb_steps_after_quad_conv.c);
     - the residual mathematical error after those steps is <= 2^{-prec}
       (borchardt_mean.c adds it with arb_add_error_2exp_si).

   The certified-decision rule of the module binds this port: every branch is a
   certified comparison on a rigorous bound, never a midpoint test. If the
   spread cannot be certified below 1/7 within the step cap, or a square root
   cannot be certified to have positive real part, we return 2 (precision
   insufficient) rather than a guessed 1. HDME's borchardt_mean returns a bare
   failure flag in these cases; the acb_ppav.rst contract also mentions setting
   r to a coarse disc of radius ||a||_infty, but a coarse enclosure returned as
   1 would violate the certified-decision rule, so we return 2 and leave r
   undefined, matching acb_ppav_agm_half_plane's own return-2 semantics. */

#include "fmpz.h"
#include "arb.h"
#include "acb.h"
#include "acb_theta.h"
#include "acb_ppav.h"

/* One Borchardt step in place: extract square roots with positive real part,
   apply acb_theta_agm_mul, then divide by 2^g. Returns 0 if some square root
   cannot be certified to have positive real part at this precision. */
static int
_agm_step(acb_ptr b, slong g, slong prec)
{
    slong n = (WORD(1) << g);
    slong i;
    int ok = 1;
    acb_ptr roots = _acb_vec_init(n);

    for (i = 0; i < n; i++)
    {
        acb_sqrt(&roots[i], &b[i], prec);
        ok = ok && arb_is_positive(acb_realref(&roots[i]));
    }
    if (ok)
    {
        acb_theta_agm_mul(b, roots, roots, g, 0, prec);
        _acb_vec_scalar_mul_2exp_si(b, b, n, -g);
    }

    _acb_vec_clear(roots, n);
    return ok;
}

/* m0 = min_j lbound Re(b_j); M0 = max_j ubound |b_j|;
   Delta0 = sum_{i != j} ubound |b_i - b_j|. */
static void
_agm_stats(arb_t m0, arb_t M0, arb_t Delta0, acb_srcptr b, slong n, slong prec)
{
    slong i, j;
    arf_t lo, hi, acc, t;
    acb_t d;

    arf_init(lo);
    arf_init(hi);
    arf_init(acc);
    arf_init(t);
    acb_init(d);

    arb_get_lbound_arf(lo, acb_realref(&b[0]), prec);
    acb_get_abs_ubound_arf(hi, &b[0], prec);
    for (i = 1; i < n; i++)
    {
        arb_get_lbound_arf(t, acb_realref(&b[i]), prec);
        arf_min(lo, lo, t);
        acb_get_abs_ubound_arf(t, &b[i], prec);
        arf_max(hi, hi, t);
    }
    arb_set_arf(m0, lo);
    arb_set_arf(M0, hi);

    arf_zero(acc);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                acb_sub(d, &b[i], &b[j], prec);
                acb_get_abs_ubound_arf(t, d, prec);
                arf_add(acc, acc, t, prec, ARF_RND_CEIL);
            }
        }
    }
    arb_set_arf(Delta0, acc);

    arf_clear(lo);
    arf_clear(hi);
    arf_clear(acc);
    arf_clear(t);
    acb_clear(d);
}

/* Number of Borchardt steps guaranteed to reach quadratic convergence,
   ceil( log(m0/(7 Delta0)) / log(1 - 2^{-g}) ). Returns 0 if the vector is not
   certifiably a valid AGM start (m0 not positive, or any bound non-finite). */
static int
_agm_nb_before(fmpz_t nb, acb_srcptr b, slong g, slong prec)
{
    slong n = (WORD(1) << g);
    arb_t m0, M0, Delta0, num, den;
    arf_t sup;
    int res = 1;

    arb_init(m0);
    arb_init(M0);
    arb_init(Delta0);
    arb_init(num);
    arb_init(den);
    arf_init(sup);

    _agm_stats(m0, M0, Delta0, b, n, prec);

    if (!arb_is_positive(m0) || !arb_is_finite(m0) || !arb_is_finite(M0)
        || !arb_is_finite(Delta0) || arb_contains_negative(Delta0))
        res = 0;

    if (res)
    {
        if (arb_is_zero(Delta0))
        {
            fmpz_zero(nb);
        }
        else
        {
            arb_div(num, m0, Delta0, prec);
            arb_div_si(num, num, 7, prec);
            arb_log(num, num, prec);

            arb_one(den);
            arb_mul_2exp_si(den, den, -g);
            arb_neg(den, den);
            arb_add_si(den, den, 1, prec);
            arb_log(den, den, prec);

            arb_div(num, num, den, prec);
            if (arb_is_negative(num))
                arb_zero(num);

            arb_get_ubound_arf(sup, num, prec);
            arf_ceil(sup, sup);
            arf_get_fmpz(nb, sup, ARF_RND_NEAR);
        }
    }

    arb_clear(m0);
    arb_clear(M0);
    arb_clear(Delta0);
    arb_clear(num);
    arb_clear(den);
    arf_clear(sup);
    return res;
}

/* Number of further steps to bring the residual error to <= 2^{-prec} once in
   the quadratic regime, ceil( log2( log2(M0/7) + prec + 1 ) ). */
static void
_agm_nb_after(fmpz_t nb, acb_srcptr b, slong g, slong prec)
{
    slong n = (WORD(1) << g);
    arb_t m0, M0, Delta0, num;
    arf_t sup;

    arb_init(m0);
    arb_init(M0);
    arb_init(Delta0);
    arb_init(num);
    arf_init(sup);

    _agm_stats(m0, M0, Delta0, b, n, prec);

    arb_div_si(num, M0, 7, prec);
    arb_log_base_ui(num, num, 2, prec);
    arb_add_si(num, num, prec, prec);
    arb_add_si(num, num, 1, prec);
    arb_log_base_ui(num, num, 2, prec);

    arb_get_ubound_arf(sup, num, prec);
    arf_ceil(sup, sup);
    arf_get_fmpz(nb, sup, ARF_RND_NEAR);

    arb_clear(m0);
    arb_clear(M0);
    arb_clear(Delta0);
    arb_clear(num);
    arf_clear(sup);
}

/* b = exp(-2 pi i beta) * a, with tib = 2 pi i beta stored for the inverse
   rescaling at the end. */
static void
_agm_rescale(acb_ptr b, acb_t tib, acb_srcptr a, const arf_t beta, slong n,
    slong prec)
{
    acb_t scal;
    acb_init(scal);

    acb_zero(tib);
    arb_const_pi(acb_imagref(tib), prec);
    arb_mul_2exp_si(acb_imagref(tib), acb_imagref(tib), 1);
    arb_mul_arf(acb_imagref(tib), acb_imagref(tib), beta, prec);

    acb_neg(scal, tib);
    acb_exp(scal, scal, prec);
    _acb_vec_scalar_mul(b, a, n, scal, prec);

    acb_clear(scal);
}

int
acb_ppav_agm(acb_ptr r, acb_srcptr a, const arf_t eps, slong g, slong prec)
{
    slong n = (WORD(1) << g);
    arf_t beta, delta, thr;
    acb_t tib, scal;
    acb_ptr b;
    fmpz_t nb_before, nb_after;
    slong wp, pp, i, cap, na;
    int hp, reached, ret;

    arf_init(beta);

    hp = acb_ppav_agm_half_plane(beta, a, eps, g, prec);
    if (hp != 1)
    {
        arf_clear(beta);
        return hp;   /* 0 or 2 */
    }

    arf_init(delta);
    arf_init(thr);
    acb_init(tib);
    acb_init(scal);
    b = _acb_vec_init(n);
    fmpz_init(nb_before);
    fmpz_init(nb_after);

    /* Sizing pass: derive the step counts at a modest precision, then choose a
       working precision whose guard covers the (empirically ~2 bits/step)
       rounding loss so the final enclosure keeps ~prec relative bits. */
    pp = prec + 32;
    _agm_rescale(b, tib, a, beta, n, pp);
    if (!_agm_nb_before(nb_before, b, g, pp))
    {
        ret = 2;
        goto cleanup;
    }
    _agm_nb_after(nb_after, b, g, pp);
    wp = prec + 32 + 4 * (fmpz_get_si(nb_before) + fmpz_get_si(nb_after));

    /* Working pass. */
    _agm_rescale(b, tib, a, beta, n, wp);

    arf_set_si(thr, 1);
    arf_div_si(thr, thr, 7, wp, ARF_RND_DOWN);   /* thr <= 1/7 */

    cap = fmpz_get_si(nb_before);
    reached = 0;
    for (i = 0; i < cap; i++)
    {
        acb_ppav_agm_max_diff(delta, b, g, 30);
        if (arf_cmp(delta, thr) < 0)
        {
            reached = 1;
            break;
        }
        if (!_agm_step(b, g, wp))
        {
            ret = 2;
            goto cleanup;
        }
    }
    if (!reached)
    {
        acb_ppav_agm_max_diff(delta, b, g, 30);
        reached = (arf_cmp(delta, thr) < 0);
    }
    if (!reached)
    {
        ret = 2;   /* spread not certifiably below 1/7 at this precision */
        goto cleanup;
    }

    _agm_nb_after(nb_after, b, g, wp);
    na = fmpz_get_si(nb_after);
    for (i = 0; i < na; i++)
    {
        if (!_agm_step(b, g, wp))
        {
            ret = 2;
            goto cleanup;
        }
    }

    /* Certified enclosure of the limit, then rescale back by exp(+2 pi i beta). */
    acb_set(r, &b[0]);
    arb_add_error_2exp_si(acb_realref(r), -prec);
    arb_add_error_2exp_si(acb_imagref(r), -prec);
    acb_exp(scal, tib, wp);
    acb_mul(r, r, scal, wp);
    ret = 1;

cleanup:
    arf_clear(beta);
    arf_clear(delta);
    arf_clear(thr);
    acb_clear(tib);
    acb_clear(scal);
    _acb_vec_clear(b, n);
    fmpz_clear(nb_before);
    fmpz_clear(nb_after);
    return ret;
}
