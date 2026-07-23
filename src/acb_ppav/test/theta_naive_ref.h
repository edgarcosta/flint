/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef THETA_NAIVE_REF_H
#define THETA_NAIVE_REF_H

#include "arb.h"
#include "acb.h"
#include "acb_mat.h"

/* Independent, convention-explicit theta constant; a, b in {0,1}^2.

   Returns a certified enclosure of theta_{a,b}(0,tau)^2, with

       theta_{a,b}(0,tau) = sum_{n in Z^2} exp(pi i (v^T tau v + v . b)),  v = n + a/2,

   which is the acb_theta convention (acb_theta.rst) specialized to z = 0:
   theta_{a,b}(0,tau) = sum_{n in Z^2 + a/2} exp(pi i n^T tau n + pi i n^T b).
   The lattice sum is truncated to |n_i| <= N and a rigorous tail bound is
   folded into the radius before squaring. This is the one convention-free
   reference in the periods vertical, so the tail bound must be sound.

   Tail bound. Write tau = X + iY with Y = Im(tau) positive definite, and let
   lambda > 0 be a rigorous lower bound on the smallest eigenvalue of Y. Since
   |exp(pi i (v^T tau v + v.b))| = exp(-pi v^T Y v) <= exp(-pi lambda |v|^2), and
   with q = exp(-pi lambda) in (0,1) the bound factorizes over coordinates. The
   full 1D sum F = sum_{n in Z} q^{(n+a/2)^2} satisfies F <= 2/(1-q), and the 1D
   tail t_N = sum_{|n| >= N+1} q^{(n+a/2)^2} satisfies t_N <= 2 q^{(N-1)^2}/(1-q).
   The dropped set {n : max(|n_i|) >= N+1} is covered by the two coordinate
   slabs, so

       sum_{dropped} |term| <= 2 F t_N <= 8 q^{(N-1)^2} / (1-q)^2,

   the quantity added to the radius. It grows as lambda -> 0, so replacing
   lambda_min by a lower bound keeps it conservative. */
static void
_acb_ppav_naive_theta2_ab(acb_t out, const int a[2], const int b[2],
                          const acb_mat_t tau, slong prec)
{
    slong N, n0, n1, wp = prec + 32;
    arb_t lambda, tail;
    acb_t th, term, quad, pii;
    arb_init(lambda); arb_init(tail);
    acb_init(th); acb_init(term); acb_init(quad); acb_init(pii);
    acb_const_pi(pii, wp); acb_mul_onei(pii, pii);          /* pii = pi*i */

    /* lambda = rigorous lower bound on the smallest eigenvalue of Y = Im(tau):
       for Y = [[y11,y12],[y12,y22]], lambda_min = (tr - sqrt(tr^2 - 4 det))/2. */
    {
        arb_t y11, y22, y12, tr, det, disc;
        arf_t lb;
        arb_init(y11); arb_init(y22); arb_init(y12);
        arb_init(tr); arb_init(det); arb_init(disc);
        arf_init(lb);
        arb_set(y11, acb_imagref(acb_mat_entry(tau, 0, 0)));
        arb_set(y22, acb_imagref(acb_mat_entry(tau, 1, 1)));
        arb_set(y12, acb_imagref(acb_mat_entry(tau, 0, 1)));
        arb_add(tr, y11, y22, wp);
        arb_mul(det, y11, y22, wp); arb_submul(det, y12, y12, wp);
        arb_mul(disc, tr, tr, wp); arb_submul_si(disc, det, 4, wp);
        arb_sqrt(disc, disc, wp);
        arb_sub(lambda, tr, disc, wp); arb_mul_2exp_si(lambda, lambda, -1);
        arb_get_lbound_arf(lb, lambda, wp);     /* collapse to a lower bound */
        arb_set_arf(lambda, lb);
        arf_clear(lb);
        arb_clear(y11); arb_clear(y22); arb_clear(y12);
        arb_clear(tr); arb_clear(det); arb_clear(disc);
    }

    /* Choose N so that tail < 2^{-prec-8}, doubling with a safety cap. If the
       cap is hit the (possibly large) tail is still added honestly below. */
    N = 4;
    while (1)
    {
        arb_t c, num, den;
        arb_init(c); arb_init(num); arb_init(den);
        arb_const_pi(c, wp); arb_mul(c, c, lambda, wp);     /* c = pi*lambda */
        arb_mul_si(num, c, -(N - 1) * (N - 1), wp);
        arb_exp(num, num, wp); arb_mul_si(num, num, 8, wp); /* 8 exp(-c (N-1)^2) */
        arb_neg(den, c); arb_exp(den, den, wp);             /* exp(-c) */
        arb_sub_si(den, den, 1, wp); arb_neg(den, den);     /* 1 - exp(-c) */
        arb_sqr(den, den, wp);                              /* (1 - exp(-c))^2 */
        arb_div(tail, num, den, wp);
        arb_clear(c); arb_clear(num); arb_clear(den);
        if (arf_cmpabs_2exp_si(arb_midref(tail), -(prec + 8)) < 0) break;
        N *= 2;
        if (N > 4096) break;                                /* safety cap */
    }

    /* Main sum over the box |n_i| <= N. */
    acb_zero(th);
    for (n0 = -N; n0 <= N; n0++)
      for (n1 = -N; n1 <= N; n1++)
      {
        arb_t v0, v1;
        acb_t t;
        arb_init(v0); arb_init(v1); acb_init(t);
        arb_set_si(v0, 2 * n0 + a[0]); arb_mul_2exp_si(v0, v0, -1);   /* v = n + a/2 */
        arb_set_si(v1, 2 * n1 + a[1]); arb_mul_2exp_si(v1, v1, -1);
        /* quad = v^T tau v */
        acb_mul_arb(quad, acb_mat_entry(tau, 0, 0), v0, wp); acb_mul_arb(quad, quad, v0, wp);
        acb_mul_arb(t, acb_mat_entry(tau, 1, 1), v1, wp); acb_mul_arb(t, t, v1, wp);
        acb_add(quad, quad, t, wp);
        acb_mul_arb(t, acb_mat_entry(tau, 0, 1), v0, wp); acb_mul_arb(t, t, v1, wp);
        acb_mul_2exp_si(t, t, 1); acb_add(quad, quad, t, wp);        /* + 2 tau01 v0 v1 */
        /* quad += v . b   (v0 reused as scratch to hold v . b) */
        arb_mul_si(v0, v0, b[0], wp); arb_addmul_si(v0, v1, b[1], wp);
        acb_add_arb(quad, quad, v0, wp);
        acb_mul(term, pii, quad, wp); acb_exp(term, term, wp);
        acb_add(th, th, term, wp);
        acb_clear(t); arb_clear(v0); arb_clear(v1);
      }

    /* Fold the tail into the radius of |th|, then square. */
    acb_add_error_arb(th, tail);
    acb_sqr(out, th, prec);

    arb_clear(lambda); arb_clear(tail);
    acb_clear(th); acb_clear(term); acb_clear(quad); acb_clear(pii);
}

#endif
