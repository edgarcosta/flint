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

/* Rosenhain invariants of a genus 2 curve from its six ordered Weierstrass
   points w, port of hdme igusa/thomae_reorder.c and igusa/thomae_rosenhain.c
   (commit 309bfd25).

   perm in [0, 720) (PRECONDITION) is decoded by the factorial mixed-radix
   scheme into an image list; new[images[k]] = w[k] reorders w; then the
   Mobius map sending (new[5], new[0], new[1]) -> (infinity, 0, 1) is applied
   and ros = (new[2], new[3], new[4]).  An out-of-range perm goes out of
   bounds either way: perm >= 720 makes images[0] = perm/120 >= 6 an
   out-of-bounds write, perm < 0 makes some image index negative.  Callers
   pass only [0, 720), the range HDME's thomae_correct_signs iterates.

   Orientation: new[images[k]] = w[k] places w[k] at position images[k], so
   new[j] = w[images^{-1}(j)].  The sigma of the "bring w_{sigma(5)},
   w_{sigma(0)}, w_{sigma(1)} to infinity, 0, 1" phrasing is therefore
   images^{-1}, not images.  HDME's assignment is kept byte-compatibly so a
   future search stays compatible with HDME's driver.

   No return value and no branching on ball midpoints: division by a ball
   containing zero yields an indeterminate ball that propagates, the correct
   contract here.  ros (length 3) and w (length 6) may overlap: w is fully
   copied into an internal buffer before any entry of ros is written. */

/* Factorial mixed-radix permutation decode (HDME perm_nb_to_images,
   thomae_reorder.c:4).  DEVIATION from HDME (documented, safe): HDME's base
   case is n==0, which writes images[0]; called with n=6 the recursion
   descends to an n==0 frame at offset +6 and writes images[6], one past the
   length-6 buffer.  We stop at n==1 (whose sole Lehmer digit is 0) and never
   reach the n==0 write.  Numerically identical: at n==1 the incoming nb is
   always 0. */
static void
_perm_images(slong * images, slong nb, slong n)
{
    slong fact = 1, k, q, r;

    if (n <= 1)
    {
        if (n == 1)
            images[0] = 0;
        return;
    }
    for (k = 2; k < n; k++)
        fact *= k;                  /* (n-1)! */
    q = nb / fact;
    r = nb % fact;
    images[0] = q;
    _perm_images(images + 1, r, n - 1);
    for (k = 1; k < n; k++)
        if (images[k] >= q)
            images[k] += 1;
}

void
acb_ppav_g2_rosenhain(acb_ptr ros, acb_srcptr w, slong perm, slong prec)
{
    slong images[6];
    acb_ptr new_w;
    acb_ptr proj;
    acb_t den;
    slong k;

    _perm_images(images, perm, 6);

    new_w = _acb_vec_init(6);
    proj = _acb_vec_init(5);
    acb_init(den);

    for (k = 0; k < 6; k++)
        acb_set(&new_w[images[k]], &w[k]);

    /* bring new_w[5] to infinity */
    for (k = 0; k < 5; k++)
    {
        acb_sub(&proj[k], &new_w[k], &new_w[5], prec);
        acb_inv(&proj[k], &proj[k], prec);
    }
    /* bring new_w[0], new_w[1] to 0, 1 */
    acb_sub(den, &proj[1], &proj[0], prec);
    acb_inv(den, den, prec);
    for (k = 2; k < 5; k++)
    {
        acb_sub(&proj[k], &proj[k], &proj[0], prec);
        acb_mul(&proj[k], &proj[k], den, prec);
    }

    acb_set(&ros[0], &proj[2]);
    acb_set(&ros[1], &proj[3]);
    acb_set(&ros[2], &proj[4]);

    _acb_vec_clear(new_w, 6);
    _acb_vec_clear(proj, 5);
    acb_clear(den);
}
