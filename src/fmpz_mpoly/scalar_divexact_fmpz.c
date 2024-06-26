/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpoly.h"
#include "fmpz_mpoly.h"

void fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                    const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if (A != B)
    {
        N = mpoly_words_per_exp(B->bits, ctx->minfo);
        fmpz_mpoly_fit_length(A, B->length, ctx);
        fmpz_mpoly_fit_bits(A, B->bits, ctx);
        A->bits = B->bits;
        mpn_copyi(A->exps, B->exps, N*B->length);
    }
    _fmpz_vec_scalar_divexact_fmpz(A->coeffs, B->coeffs, B->length, c);
    _fmpz_mpoly_set_length(A, B->length, ctx);
}

void fmpz_mpoly_scalar_divexact_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                           ulong c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_ui(t, c);
    fmpz_mpoly_scalar_divexact_fmpz(A, B, t, ctx);
    fmpz_clear(t);
}

void fmpz_mpoly_scalar_divexact_si(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                           slong c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_si(t, c);
    fmpz_mpoly_scalar_divexact_fmpz(A, B, t, ctx);
    fmpz_clear(t);
}
