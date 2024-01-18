/*
    Copyright (C) 2017, 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "n_poly.h"
#include "fq_nmod_mpoly.h"

int fq_nmod_mpoly_equal_fq_nmod(
    const fq_nmod_mpoly_t A,
    const fq_nmod_t c,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;

    if (fq_nmod_is_zero(c, ctx->fqctx))
        return A->length == 0;

    if (A->length != 1)
        return 0;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (!mpoly_monomial_is_zero(A->exps + N*0, N))
        return 0;

    return n_fq_equal_fq_nmod(A->coeffs, c, ctx->fqctx);
}
