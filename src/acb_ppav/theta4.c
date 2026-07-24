/*
    Copyright (C) 2026 Jean Kieffer
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Ported from hdme igusa/thomae_theta4.c (commit 309bfd25), (C) Jean Kieffer.

   Given the Rosenhain invariants (l, m, n) = (ros[0], ros[1], ros[2]), Thomae's
   formulas give the ten even fourth-power theta constants (normalized so the
   label-0 value is 1); the six odd entries are 0.

   Two deviations from HDME:
     (1) HDME stores each even value at theta_char_set_label_g2(L) = 8a0+4a1+2b0+b1,
         which is already acb_theta's k-index.  A[L] is that label-to-index map
         precomputed as a constant table, purely to avoid porting n_revbin /
         theta_char_set_ab; -1 marks the six odd labels and is never dereferenced.
         This is a numerical no-op relative to HDME.
     (2) The slong return has no HDME analogue; it is the rst's per-characteristic
         16-bit flag (below).

   Precondition: th4 (length 16) must not alias ros (length 3). */

#include "arb.h"
#include "acb.h"
#include "acb_ppav.h"

slong
acb_ppav_g2_theta4(acb_ptr th4, acb_srcptr ros, slong prec)
{
    static const slong A[16] = { 0, 2, 1, 3,  8, -1, 9, -1,  4, 6, -1, -1,  12, -1, -1, 15 };
    acb_t l, m, n;
    acb_t num, den, aux;
    slong flag, k;

    acb_init(l);
    acb_init(m);
    acb_init(n);
    acb_init(num);
    acb_init(den);
    acb_init(aux);

    acb_set(l, &ros[0]);
    acb_set(m, &ros[1]);
    acb_set(n, &ros[2]);

    for (k = 0; k < 16; k++)
        acb_zero(&th4[k]);

    /* Thomae's formulas (hdme thomae_theta4.c), labels remapped through A[]. */

    /* label 0 */
    acb_one(&th4[A[0]]);

    /* label 1: m (l-1) (n-1) / (l n (m-1)) */
    acb_set(num, m);
    acb_sub_si(aux, l, 1, prec);
    acb_mul(num, num, aux, prec);
    acb_sub_si(aux, n, 1, prec);
    acb_mul(num, num, aux, prec);
    acb_mul(den, l, n, prec);
    acb_sub_si(aux, m, 1, prec);
    acb_mul(den, den, aux, prec);
    acb_div(&th4[A[1]], num, den, prec);

    /* label 2: m (l-1) (n-m) / (l (m-1) (n-l)) */
    acb_set(num, m);
    acb_sub_si(aux, l, 1, prec);
    acb_mul(num, num, aux, prec);
    acb_sub(aux, n, m, prec);
    acb_mul(num, num, aux, prec);
    acb_set(den, l);
    acb_sub_si(aux, m, 1, prec);
    acb_mul(den, den, aux, prec);
    acb_sub(aux, n, l, prec);
    acb_mul(den, den, aux, prec);
    acb_div(&th4[A[2]], num, den, prec);

    /* label 4: m / (l n) */
    acb_set(num, m);
    acb_mul(den, l, n, prec);
    acb_div(&th4[A[4]], num, den, prec);

    /* label 8: m (l-m) (n-1) / (n (m-1) (l-n)) */
    acb_set(num, m);
    acb_sub(aux, l, m, prec);
    acb_mul(num, num, aux, prec);
    acb_sub_si(aux, n, 1, prec);
    acb_mul(num, num, aux, prec);
    acb_set(den, n);
    acb_sub_si(aux, m, 1, prec);
    acb_mul(den, den, aux, prec);
    acb_sub(aux, l, n, prec);
    acb_mul(den, den, aux, prec);
    acb_div(&th4[A[8]], num, den, prec);

    /* label 6: th4[A[2]] / (n^2 th4[A[4]]) */
    acb_set(num, &th4[A[2]]);
    acb_sqr(den, n, prec);
    acb_mul(den, den, &th4[A[4]], prec);
    acb_div(&th4[A[6]], num, den, prec);

    /* label 12: th4[A[8]] / (l^2 th4[A[4]]) */
    acb_set(num, &th4[A[8]]);
    acb_sqr(den, l, prec);
    acb_mul(den, den, &th4[A[4]], prec);
    acb_div(&th4[A[12]], num, den, prec);

    /* label 3: th4[A[4]] th4[A[6]] (n-1)^2 / th4[A[1]] */
    acb_mul(num, &th4[A[4]], &th4[A[6]], prec);
    acb_sub_si(aux, n, 1, prec);
    acb_sqr(aux, aux, prec);
    acb_mul(num, num, aux, prec);
    acb_set(den, &th4[A[1]]);
    acb_div(&th4[A[3]], num, den, prec);

    /* label 9: th4[A[4]] th4[A[12]] (l-1)^2 / th4[A[1]] */
    acb_mul(num, &th4[A[4]], &th4[A[12]], prec);
    acb_sub_si(aux, l, 1, prec);
    acb_sqr(aux, aux, prec);
    acb_mul(num, num, aux, prec);
    acb_set(den, &th4[A[1]]);
    acb_div(&th4[A[9]], num, den, prec);

    /* label 15: th4[A[1]] th4[A[12]] (n-m)^2 / (th4[A[2]] (n-1)^2) */
    acb_mul(num, &th4[A[1]], &th4[A[12]], prec);
    acb_sub(aux, n, m, prec);
    acb_sqr(aux, aux, prec);
    acb_mul(num, num, aux, prec);
    acb_set(den, &th4[A[2]]);
    acb_sub_si(aux, n, 1, prec);
    acb_sqr(aux, aux, prec);
    acb_mul(den, den, aux, prec);
    acb_div(&th4[A[15]], num, den, prec);

    /* Branch-cut flag: bit k is set iff th4[k] certainly lies on acb_sqrt's
       branch cut, the negative real axis (imag straddles 0 and real is entirely
       < 0).  This is the exact rotation trigger of hdme theta/borchardt_root_ui.c.
       The rst (theta4 section) says "negative imaginary axis"; the branch cut of
       the principal square root is the negative REAL axis, so this follows the
       branch-cut-consistent reading (flagged to Jean as a proposed rst fix). */
    flag = 0;
    for (k = 0; k < 16; k++)
        if (arb_contains_zero(acb_imagref(&th4[k]))
            && arb_is_negative(acb_realref(&th4[k])))
            flag |= (WORD(1) << k);

    acb_clear(l);
    acb_clear(m);
    acb_clear(n);
    acb_clear(num);
    acb_clear(den);
    acb_clear(aux);

    return flag;
}
