/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_PPAV_THOMAE_REF_H
#define ACB_PPAV_THOMAE_REF_H

#include "acb.h"

/* Recover the Rosenhain triple (l, m, n) from squared theta constants.
   th2 has length 16 in acb_theta ordering (k = 8 a0 + 4 a1 + 2 b0 + b1);
   only the even indices {0,1,4,8,9,12} are read, at HARD-CODED positions
   (deliberately not routed through any label map, so this helper is an
   independent witness of the ordering).  Scale-invariant: with
   t_k = th2[k]/th2[0],
       l = t_4 / (t_8 t_12)   = th2[4] th2[0] / (th2[8] th2[12])
       m = t_1 t_4 / (t_9 t_12) = th2[1] th2[4] / (th2[9] th2[12])
       n = t_1 / (t_8 t_9)    = th2[1] th2[0] / (th2[8] th2[9])
   ros = (l, m, n), length 3.  ros must not alias th2.  Division by a ball
   containing zero yields an indeterminate ball, which propagates to the
   caller's overlap assertions; no branching here. */
static void
_ros_from_theta2(acb_ptr ros, acb_srcptr th2, slong prec)
{
    acb_t num, den;

    acb_init(num);
    acb_init(den);

    /* l = th2[4] th2[0] / (th2[8] th2[12]) */
    acb_mul(num, &th2[4], &th2[0], prec);
    acb_mul(den, &th2[8], &th2[12], prec);
    acb_div(&ros[0], num, den, prec);

    /* m = th2[1] th2[4] / (th2[9] th2[12]) */
    acb_mul(num, &th2[1], &th2[4], prec);
    acb_mul(den, &th2[9], &th2[12], prec);
    acb_div(&ros[1], num, den, prec);

    /* n = th2[1] th2[0] / (th2[8] th2[9]) */
    acb_mul(num, &th2[1], &th2[0], prec);
    acb_mul(den, &th2[8], &th2[9], prec);
    acb_div(&ros[2], num, den, prec);

    acb_clear(num);
    acb_clear(den);
}

#endif
