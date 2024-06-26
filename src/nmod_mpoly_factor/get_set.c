/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "nmod_mpoly_factor.h"

slong nmod_mpoly_factor_get_exp_si(nmod_mpoly_factor_t f,
                                           slong i, const nmod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    FLINT_ASSERT(i < (ulong) f->num);
    return fmpz_get_si(f->exp + i);
}

void nmod_mpoly_factor_set(nmod_mpoly_factor_t res,
                     const nmod_mpoly_factor_t fac, const nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (res == fac)
        return;

    nmod_mpoly_factor_fit_length(res, fac->num, ctx);
    res->constant = fac->constant;
    for (i = 0; i < fac->num; i++)
    {
        nmod_mpoly_set(res->poly + i, fac->poly + i, ctx);
        fmpz_set(res->exp + i, fac->exp + i);
    }
    res->num = fac->num;
}

void _nmod_mpoly_get_lead0(
    nmod_mpoly_t c,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyl_lead_coeff(c, A, 1, ctx);
}

void _nmod_mpoly_set_lead0(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_t c,
    const nmod_mpoly_ctx_t ctx)
{
    slong deg;
    nmod_mpoly_t t, g;

    nmod_mpoly_init(t, ctx);
    nmod_mpoly_init(g, ctx);

    deg = nmod_mpoly_degree_si(B, 0, ctx);
    FLINT_ASSERT(deg >= 0);
    nmod_mpoly_gen(g, 0, ctx);
    nmod_mpoly_pow_ui(g, g, deg, ctx);
    _nmod_mpoly_get_lead0(t, B, ctx);
    nmod_mpoly_sub(t, c, t, ctx);
    nmod_mpoly_mul(t, t, g, ctx);
    nmod_mpoly_add(A, B, t, ctx);

    nmod_mpoly_clear(t, ctx);
    nmod_mpoly_clear(g, ctx);
}
