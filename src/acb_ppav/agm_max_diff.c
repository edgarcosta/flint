/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_ppav.h"

void
acb_ppav_agm_max_diff(arf_t delta, acb_srcptr a, slong g, slong prec)
{
    slong n = (WORD(1) << g);
    slong j;
    arf_t num, den, q;
    arb_t t;

    arf_init(num);
    arf_init(den);
    arf_init(q);
    arb_init(t);

    /* Certified upper bound on the ratio: bound the numerator from above, the
       denominator from below, and round the quotient up. */
    acb_abs(t, &a[0], prec);
    arb_get_lbound_arf(den, t, prec);

    arf_zero(num);
    for (j = 1; j < n; j++)
    {
        acb_t d;
        acb_init(d);
        acb_sub(d, &a[j], &a[0], prec);
        acb_abs(t, d, prec);
        arb_get_ubound_arf(q, t, prec);
        arf_max(num, num, q);
        acb_clear(d);
    }

    if (arf_sgn(den) <= 0)   /* a0 ball touches zero: ratio is unbounded */
        arf_pos_inf(delta);
    else
        arf_div(delta, num, den, prec, ARF_RND_UP);

    arf_clear(num);
    arf_clear(den);
    arf_clear(q);
    arb_clear(t);
}
