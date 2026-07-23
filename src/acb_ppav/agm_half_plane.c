/*
    Copyright (C) 2026 Jean Kieffer
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_ppav.h"

/* Ported and adapted from HDME theta/borchardt_excl_half_planes.c and
   theta/borchardt_mean_invalid.c (both (C) Jean Kieffer). HDME parametrizes
   half planes by arg/pi in [-1,1] and computes only the excluded (certain-fail)
   region with eps = 0. Here we parametrize by beta in R/Z = [0,1), carry a
   strictly positive margin eps, and compute BOTH the certain-yes and
   certain-no regions per ball, as required by the acb_ppav.rst contract. The
   substance carried over is the per-ball endpoint arithmetic: acb_arg with the
   negative-real-axis branch-cut special case, directed rounding via
   arb_get_{l,u}bound_arf, and the treatment of a zero-containing ball as
   contributing to neither region. */

/* Fractional part in [0, 1): r = x - floor(x), computed exactly. */
static void
_agm_hp_frac(arf_t r, const arf_t x)
{
    arf_t f;
    arf_init(f);
    arf_floor(f, x);
    arf_sub(r, x, f, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_clear(f);
}

/* For a single ball z, using cc = arccos(eps)/(2 pi) and pi2 = 2 pi (both arb
   enclosures), compute:
     - the closed certain-yes arc [yesL, yesR] (lifted reals) and set *has_yes
       to 1 iff it is nonempty (yesL <= yesR); the set of beta for which
       Re(e^{-2 pi i beta} w) >= eps |w| certainly holds for every w in z;
     - the open certain-no arc (noS, noS + noLen) with noS in [0,1) and
       noLen in (0,1), and set *has_no accordingly; the set of beta for which
       that inequality certainly fails for every w in z.
   A ball containing zero contributes to neither region. Both arcs are certified
   subsets of the true regions by the directed rounding. */
static void
_agm_hp_ball(int *has_yes, arf_t yesL, arf_t yesR,
             int *has_no, arf_t noS, arf_t noLen,
             const acb_t z, const arb_t cc, const arb_t pi2, slong prec)
{
    arb_t base, bpc, bmc, arg;
    arf_t P, Q, half;

    *has_yes = 0;
    *has_no = 0;

    if (acb_contains_zero(z))
        return;

    arb_init(base);
    arb_init(bpc);
    arb_init(bmc);
    arb_init(arg);
    arf_init(P);
    arf_init(Q);
    arf_init(half);

    if (arb_contains_zero(acb_imagref(z)) && arb_is_negative(acb_realref(z)))
    {
        /* z straddles the branch cut of acb_arg (negative real axis). Negate to
           land near the positive real axis where acb_arg is tight, then shift
           the beta center by 1/2 (since arg(z) = arg(-z) + pi). */
        acb_t nz;
        acb_init(nz);
        acb_neg(nz, z);
        acb_arg(arg, nz, prec);
        acb_clear(nz);
        arb_div(base, arg, pi2, prec);
        arf_set_si(half, 1);
        arf_mul_2exp_si(half, half, -1);
        arb_add_arf(base, base, half, prec);
    }
    else
    {
        acb_arg(arg, z, prec);
        arb_div(base, arg, pi2, prec);
    }

    arb_add(bpc, base, cc, prec);   /* beta center + half cone width */
    arb_sub(bmc, base, cc, prec);   /* beta center - half cone width */

    /* Certain-yes [yesL, yesR]: widest interval contained in every point's yes
       interval. Round the left endpoint up and the right endpoint down. */
    arb_get_ubound_arf(yesL, bmc, prec);
    arb_get_lbound_arf(yesR, bpc, prec);
    *has_yes = (arf_cmp(yesL, yesR) <= 0);

    /* Certain-no (noS, noS + noLen): complement of the union over the ball of
       the closed yes intervals, i.e. (ubound(base + cc), lbound(base - cc) + 1),
       an open arc. */
    arb_get_ubound_arf(P, bpc, prec);
    arb_get_lbound_arf(Q, bmc, prec);
    arf_add_si(Q, Q, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_sub(noLen, Q, P, ARF_PREC_EXACT, ARF_RND_DOWN);
    if (arf_sgn(noLen) > 0)
    {
        _agm_hp_frac(noS, P);
        *has_no = 1;
    }

    arb_clear(base);
    arb_clear(bpc);
    arb_clear(bmc);
    arb_clear(arg);
    arf_clear(P);
    arf_clear(Q);
    arf_clear(half);
}

/* Is x in [0,1) strictly inside some open arc (noS[j], noS[j] + noLen[j])? */
static int
_agm_hp_point_covered(const arf_t x, const arf_struct * noS,
                      const arf_struct * noLen, const int * has_no, slong n)
{
    arf_t d;
    slong j;
    int res = 0;

    arf_init(d);
    for (j = 0; j < n; j++)
    {
        if (!has_no[j])
            continue;
        arf_sub(d, x, &noS[j], ARF_PREC_EXACT, ARF_RND_DOWN);
        _agm_hp_frac(d, d);
        if (arf_sgn(d) > 0 && arf_cmp(d, &noLen[j]) < 0)
        {
            res = 1;
            break;
        }
    }
    arf_clear(d);
    return res;
}

/* Do the open certain-no arcs cover the whole of R/Z? Coverage is constant on
   each open interval between consecutive endpoints, so it suffices to test every
   endpoint and every midpoint between consecutive endpoints. */
static int
_agm_hp_cover(const arf_struct * noS, const arf_struct * noLen,
              const int * has_no, slong n)
{
    arf_ptr cp;
    slong ncp = 0, i, j, k;
    int covered = 1;

    for (j = 0; j < n; j++)
        if (has_no[j] && arf_cmp_si(&noLen[j], 1) >= 0)
            return 1;   /* a full-circle arc covers by itself */

    cp = _arf_vec_init(2 * n);
    for (j = 0; j < n; j++)
    {
        if (!has_no[j])
            continue;
        arf_set(&cp[ncp], &noS[j]);
        ncp++;
        arf_add(&cp[ncp], &noS[j], &noLen[j], ARF_PREC_EXACT, ARF_RND_DOWN);
        _agm_hp_frac(&cp[ncp], &cp[ncp]);
        ncp++;
    }

    if (ncp == 0)
    {
        _arf_vec_clear(cp, 2 * n);
        return 0;   /* no arcs: nothing covered */
    }

    /* insertion sort the endpoints ascending */
    for (i = 1; i < ncp; i++)
        for (k = i; k > 0 && arf_cmp(&cp[k - 1], &cp[k]) > 0; k--)
            arf_swap(&cp[k - 1], &cp[k]);

    for (k = 0; k < ncp && covered; k++)
        if (!_agm_hp_point_covered(&cp[k], noS, noLen, has_no, n))
            covered = 0;

    for (k = 0; k < ncp && covered; k++)
    {
        arf_t mid, hi;
        arf_init(mid);
        arf_init(hi);
        if (k + 1 < ncp)
            arf_set(hi, &cp[k + 1]);
        else
        {
            arf_set(hi, &cp[0]);
            arf_add_si(hi, hi, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
        }
        arf_add(mid, &cp[k], hi, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul_2exp_si(mid, mid, -1);
        _agm_hp_frac(mid, mid);
        if (!_agm_hp_point_covered(mid, noS, noLen, has_no, n))
            covered = 0;
        arf_clear(mid);
        arf_clear(hi);
    }

    _arf_vec_clear(cp, 2 * n);
    return covered;
}

/* Intersect the closed certain-yes arcs [yesL[j], yesR[j]] on R/Z. Each arc has
   length < 1/2, so the running intersection is a single arc; align each new arc
   by an integer shift in {-1, 0, 1}. Returns 1 and sets beta to the midpoint
   (the point farthest from the complement) if the intersection is nonempty. */
static int
_agm_hp_intersect(arf_t beta, const arf_struct * yesL, const arf_struct * yesR,
                  const int * has_yes, slong n)
{
    arf_t Lcur, Rcur, Lj, Rj, len, Ls, Rs, nl, nr;
    slong j, m;
    int nonempty = 1;

    for (j = 0; j < n; j++)
        if (!has_yes[j])
            return 0;

    arf_init(Lcur);
    arf_init(Rcur);
    arf_init(Lj);
    arf_init(Rj);
    arf_init(len);
    arf_init(Ls);
    arf_init(Rs);
    arf_init(nl);
    arf_init(nr);

    _agm_hp_frac(Lcur, &yesL[0]);
    arf_sub(len, &yesR[0], &yesL[0], ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_add(Rcur, Lcur, len, ARF_PREC_EXACT, ARF_RND_DOWN);

    for (j = 1; j < n && nonempty; j++)
    {
        int found = 0;
        _agm_hp_frac(Lj, &yesL[j]);
        arf_sub(len, &yesR[j], &yesL[j], ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_add(Rj, Lj, len, ARF_PREC_EXACT, ARF_RND_DOWN);

        for (m = -1; m <= 1 && !found; m++)
        {
            arf_set(Ls, Lj);
            arf_add_si(Ls, Ls, m, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_add_si(Rs, Rj, m, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_max(nl, Lcur, Ls);
            arf_min(nr, Rcur, Rs);
            if (arf_cmp(nl, nr) <= 0)
            {
                arf_set(Lcur, nl);
                arf_set(Rcur, nr);
                if (arf_cmp_si(Lcur, 1) >= 0)
                {
                    arf_sub_si(Lcur, Lcur, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
                    arf_sub_si(Rcur, Rcur, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
                }
                found = 1;
            }
        }
        if (!found)
            nonempty = 0;
    }

    if (nonempty)
    {
        arf_add(beta, Lcur, Rcur, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul_2exp_si(beta, beta, -1);
        _agm_hp_frac(beta, beta);
    }

    arf_clear(Lcur);
    arf_clear(Rcur);
    arf_clear(Lj);
    arf_clear(Rj);
    arf_clear(len);
    arf_clear(Ls);
    arf_clear(Rs);
    arf_clear(nl);
    arf_clear(nr);
    return nonempty;
}

int
acb_ppav_agm_half_plane(arf_t beta, acb_srcptr a, const arf_t eps,
    slong g, slong prec)
{
    slong n = (WORD(1) << g);
    slong j;
    arb_t pi2, cc, t;
    arf_ptr yesL, yesR, noS, noLen;
    int * has_yes;
    int * has_no;
    int result;

    arb_init(pi2);
    arb_init(cc);
    arb_init(t);

    arb_const_pi(pi2, prec);
    arb_mul_2exp_si(pi2, pi2, 1);          /* 2 pi */
    arb_set_arf(t, eps);
    arb_acos(cc, t, prec);                 /* arccos(eps) */
    arb_div(cc, cc, pi2, prec);            /* half cone width in beta units */

    yesL = _arf_vec_init(n);
    yesR = _arf_vec_init(n);
    noS = _arf_vec_init(n);
    noLen = _arf_vec_init(n);
    has_yes = flint_malloc(n * sizeof(int));
    has_no = flint_malloc(n * sizeof(int));

    for (j = 0; j < n; j++)
        _agm_hp_ball(&has_yes[j], &yesL[j], &yesR[j], &has_no[j], &noS[j],
                     &noLen[j], &a[j], cc, pi2, prec);

    if (_agm_hp_intersect(beta, yesL, yesR, has_yes, n))
        result = 1;
    else if (_agm_hp_cover(noS, noLen, has_no, n))
        result = 0;
    else
        result = 2;

    _arf_vec_clear(yesL, n);
    _arf_vec_clear(yesR, n);
    _arf_vec_clear(noS, n);
    _arf_vec_clear(noLen, n);
    flint_free(has_yes);
    flint_free(has_no);
    arb_clear(pi2);
    arb_clear(cc);
    arb_clear(t);
    return result;
}
