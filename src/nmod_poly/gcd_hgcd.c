/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "gr_poly.h"

slong _nmod_poly_gcd_hgcd(mp_ptr G, mp_srcptr A, slong lenA,
                                   mp_srcptr B, slong lenB, nmod_t mod)
{
    slong cutoff = NMOD_BITS(mod) <= 8 ? NMOD_POLY_SMALL_GCD_CUTOFF : NMOD_POLY_GCD_CUTOFF;
    slong lenG = 0;
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_gcd_hgcd(G, &lenG, A, lenA, B, lenB, NMOD_POLY_HGCD_CUTOFF, cutoff, ctx));
    return lenG;
}

void nmod_poly_gcd_hgcd(nmod_poly_t G,
                             const nmod_poly_t A, const nmod_poly_t B)
{
    if (A->length < B->length)
    {
        nmod_poly_gcd_hgcd(G, B, A);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        nmod_poly_t tG;
        mp_ptr g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            nmod_poly_zero(G);
        }
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            nmod_poly_make_monic(G, A);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                nmod_poly_init2(tG, A->mod.n, FLINT_MIN(lenA, lenB));
                g = tG->coeffs;
            }
            else
            {
                nmod_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }

            lenG = _nmod_poly_gcd_hgcd(g, A->coeffs, lenA,
                                               B->coeffs, lenB, A->mod);

            if (G == A || G == B)
            {
                nmod_poly_swap(tG, G);
                nmod_poly_clear(tG);
            }
            G->length = lenG;

            if (G->length == 1)
                G->coeffs[0] = 1;
            else
                nmod_poly_make_monic(G, G);
        }
    }
}
