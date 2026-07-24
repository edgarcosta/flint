/*
    Copyright (C) 2026 Jean Kieffer
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_ppav.h"

/* Port of hdme igusa/thomae_roots.c (commit 309bfd25), generalized to
   genus 1 and 2.

   The isolation count returned by acb_poly_find_roots is the whole
   certificate: it reports how many of the roots it could enclose in
   pairwise-disjoint balls, so demanding that count equal the expected
   degree d certifies that all d Weierstrass points are separated. An
   inseparable input (e.g. a repeated root) can never reach d, so it
   returns 0 with no branching on ball midpoints. maxiter = prec matches
   HDME. The roots come out in find_roots order with no canonicalization;
   ordering is delegated to rosenhain's perm argument, so callers must
   treat w as a set.

   Deviation from HDME: thomae_roots trusts acb_poly_degree(crv) as the
   root count. Here w is a fixed-length-d buffer, so we pin d by g,
   reject g outside {1,2}, and reject deg(f) != d up front; a
   higher-degree f would make find_roots write deg(f) entries and overrun
   w. */

int
acb_ppav_weierstrass(acb_ptr w, const acb_poly_t f, slong g, slong prec)
{
    slong d;

    if (g == 1)
        d = 3;
    else if (g == 2)
        d = 6;
    else
        return 0;

    if (acb_poly_degree(f) != d)
        return 0;

    return acb_poly_find_roots(w, f, NULL, prec, prec) == d;
}
