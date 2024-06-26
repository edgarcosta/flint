/*
    Copyright (C) 2017, 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "mpoly.h"
#include "nmod_mpoly.h"

void nmod_mpoly_randtest_bits(nmod_mpoly_t A, flint_rand_t state,
                slong length, flint_bitcnt_t exp_bits, const nmod_mpoly_ctx_t ctx)
{
    ulong p = ctx->mod.n;
    slong i, j, nvars = ctx->minfo->nvars;
    fmpz * exp;
    TMP_INIT;

    TMP_START;

    exp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (j = 0; j < nvars; j++)
        fmpz_init(exp + j);

    nmod_mpoly_zero(A, ctx);
    for (i = 0; i < length; i++)
    {
        mpoly_monomial_randbits_fmpz(exp, state, exp_bits, ctx->minfo);
        _nmod_mpoly_push_exp_ffmpz(A, exp, ctx);
        A->coeffs[A->length - 1] = p > 1 ? 1 + n_randint(state, p - 1) : 0;
    }
    nmod_mpoly_sort_terms(A, ctx);
    nmod_mpoly_combine_like_terms(A, ctx);

    for (j = 0; j < nvars; j++)
        fmpz_clear(exp + j);

    TMP_END;
}

void nmod_mpoly_randtest_bound(nmod_mpoly_t A, flint_rand_t state,
                     slong length, ulong exp_bound, const nmod_mpoly_ctx_t ctx)
{
    ulong p = ctx->mod.n;
    slong i, j, nvars = ctx->minfo->nvars;
    ulong * exp;
    TMP_INIT;

    TMP_START;
    exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    nmod_mpoly_zero(A, ctx);
    for (i = 0; i < length; i++)
    {
        for (j = 0; j < nvars; j++)
            exp[j] = n_randint(state, exp_bound);
        _nmod_mpoly_push_exp_ui(A, exp, ctx);
        A->coeffs[A->length - 1] = p > 1 ? 1 + n_randint(state, p - 1) : 0;
    }
    nmod_mpoly_sort_terms(A, ctx);
    nmod_mpoly_combine_like_terms(A, ctx);

    TMP_END;
}

void nmod_mpoly_randtest_bounds(nmod_mpoly_t A, flint_rand_t state,
                  slong length, ulong * exp_bounds, const nmod_mpoly_ctx_t ctx)
{
    ulong p = ctx->mod.n;
    slong i, j, nvars = ctx->minfo->nvars;
    ulong * exp;
    TMP_INIT;

    TMP_START;
    exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    nmod_mpoly_zero(A, ctx);
    for (i = 0; i < length; i++)
    {
        for (j = 0; j < nvars; j++)
            exp[j] = n_randint(state, exp_bounds[j]);
        _nmod_mpoly_push_exp_ui(A, exp, ctx);
        A->coeffs[A->length - 1] = p > 1 ? 1 + n_randint(state, p - 1) : 0;
    }
    nmod_mpoly_sort_terms(A, ctx);
    nmod_mpoly_combine_like_terms(A, ctx);

    TMP_END;
}
