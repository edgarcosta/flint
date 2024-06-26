/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_set(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)
{
    if (A != B)
    {
        slong i, j;

        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                fmpz_poly_set(fmpz_poly_mat_entry(B, i, j),
                              fmpz_poly_mat_entry(A, i, j));
    }
}
