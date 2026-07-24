/*
    Copyright (C) 2026 Jean Kieffer
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Ported from hdme igusa/thomae_theta2.c and theta/borchardt_root_ui.c
   (commit 309bfd25), both (C) Jean Kieffer.

   Given the 16 fourth-power theta constants th4 (acb_theta ordering, normalized
   so th4[0] = 1) and the Rosenhain triple ros = (l, m, n), recover the squared
   theta constants th2 with th2[k]^2 = th4[k].  The th4[0] = 1 normalization is
   HDME's contract: the derived Thomae quotient degrees are inconsistent under
   arbitrary projective scaling, so a merely-nonzero th4[0] is not sufficient.
   Following HDME thomae_theta2.c we choose only the four "fundamental" square
   roots at the even characteristics kf = 1, 2, 4, 8 and assemble the remaining
   even entries {3, 6, 9, 12, 15} by Thomae's formulae (transcribed verbatim
   below).  th2[0] = 1 and the six odd characteristics {5, 7, 10, 11, 13, 14}
   are identically zero.

   Deviations from HDME.
     - Root SELECTION (no analogue in thomae_theta2.c, which always takes the
       principal root).  Per fundamental we choose a branch and a sign:
         neg, PER-CHARACTERISTIC, is read at bit kf (the acb_theta position
           itself): bit kf = 1 selects the ROTATED root i*sqrt(-z) in place of
           acb_sqrt(z).  This matches the 16-bit per-characteristic mask that
           theta4 returns and the header's negative_th4s array.
         signs, COMPACT, is read at bit j (j = 0..3 over kf = 1, 2, 4, 8 in
           HDME's iteration order): bit j = 1 negates the chosen root.  This is
           identical to HDME's s in [0, 16) driver in thomae_correct_signs.c
           (its `if (signs % 2 == 1) neg; signs /= 2` reads bit j).
       The rotated root is borchardt_root_ui.c's e = 2 case specialized: with
       e = 2 the scalar exp(i pi / e) = exp(i pi / 2) = i, so the whole path is
       acb_neg, acb_sqrt, acb_mul_onei.  (i sqrt(-z))^2 = z is a genuine square
       root; its imaginary part is certainly >= 0, and it avoids acb_sqrt's
       branch cut (the negative real axis) when z lies near it, preserving
       precision there.  The per-characteristic wiring of neg/signs is the
       reduce/split-and-twist design pinned by the reviewed plan, not present in
       HDME.
     - A[] hardcoding.  HDME resolves each HDME label L to its storage index by
       calling theta_char_set_label_g2(L) at runtime; that map equals the
       acb_theta k-index, so the labels appearing in thomae_theta2.c are inlined
       here as the constant table A[] below. */

#include "acb.h"
#include "acb_ppav.h"

/* HDME label L -> acb_theta k-index, i.e. theta_char_set_label_g2(L).  Only the
   ten even labels are referenced by Thomae's formulae; the six odd labels are
   never read and hold -1 as a poison value.  A 10-entry ordinal array would be
   a bug: this must be indexed by the HDME label, not by an even ordinal. */
static const slong A[16] = {
    0, 2, 1, 3,  8, -1, 9, -1,  4, 6, -1, -1,  12, -1, -1, 15
};

/* One fundamental square root of z at even characteristic kf (compact index j),
   with the neg/signs selection documented in the file header. */
static void
_theta2_fundamental_root(acb_t root, const acb_t z, slong kf, slong j,
    slong neg, slong signs, slong prec)
{
    if ((neg >> kf) & 1)
    {
        /* rotated root i*sqrt(-z): away from acb_sqrt's negative-real cut */
        acb_neg(root, z);
        acb_sqrt(root, root, prec);
        acb_mul_onei(root, root);
    }
    else
    {
        acb_sqrt(root, z, prec);
    }

    if ((signs >> j) & 1)
        acb_neg(root, root);
}

/* th2 (output) must not alias th4 or ros. */
void
acb_ppav_g2_theta2(acb_ptr th2, acb_srcptr th4, acb_srcptr ros,
    slong neg, slong signs, slong prec)
{
    slong j, kf;
    acb_t l, m, n;
    acb_t num, den, aux;
    acb_ptr res;

    acb_init(l);
    acb_init(m);
    acb_init(n);
    acb_init(num);
    acb_init(den);
    acb_init(aux);
    res = _acb_vec_init(16);

    acb_set(l, &ros[0]);
    acb_set(m, &ros[1]);
    acb_set(n, &ros[2]);

    acb_one(&res[0]);

    /* The four fundamental even characteristics kf = 1, 2, 4, 8 (compact j). */
    kf = 1;
    for (j = 0; j < 4; j++)
    {
        _theta2_fundamental_root(&res[kf], &th4[kf], kf, j, neg, signs, prec);
        kf *= 2;
    }

    /* Remaining even entries by Thomae's formulae, transcribed from HDME
       thomae_theta2.c with theta_char_set_label_g2(L) -> A[L]. */

    /* res[A[6]] = res[A[2]] / (n res[A[4]]) */
    acb_set(num, &res[A[2]]);
    acb_mul(den, n, &res[A[4]], prec);
    acb_div(&res[A[6]], num, den, prec);

    /* res[A[12]] = res[A[8]] / (l res[A[4]]) */
    acb_set(num, &res[A[8]]);
    acb_mul(den, l, &res[A[4]], prec);
    acb_div(&res[A[12]], num, den, prec);

    /* res[A[3]] = res[A[4]] res[A[6]] (n - 1) / res[A[1]] */
    acb_mul(num, &res[A[4]], &res[A[6]], prec);
    acb_sub_si(aux, n, 1, prec);
    acb_mul(num, num, aux, prec);
    acb_div(&res[A[3]], num, &res[A[1]], prec);

    /* res[A[9]] = res[A[4]] res[A[12]] (l - 1) / res[A[1]] */
    acb_mul(num, &res[A[4]], &res[A[12]], prec);
    acb_sub_si(aux, l, 1, prec);
    acb_mul(num, num, aux, prec);
    acb_div(&res[A[9]], num, &res[A[1]], prec);

    /* res[A[15]] = res[A[1]] res[A[12]] (n - m) / (res[A[2]] (n - 1)) */
    acb_mul(num, &res[A[1]], &res[A[12]], prec);
    acb_sub(aux, n, m, prec);
    acb_mul(num, num, aux, prec);
    acb_set(den, &res[A[2]]);
    acb_sub_si(aux, n, 1, prec);
    acb_mul(den, den, aux, prec);
    acb_div(&res[A[15]], num, den, prec);

    _acb_vec_set(th2, res, 16);

    acb_clear(l);
    acb_clear(m);
    acb_clear(n);
    acb_clear(num);
    acb_clear(den);
    acb_clear(aux);
    _acb_vec_clear(res, 16);
}
