/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_set(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (res != x)
    {
        ulong xfield;
        slong field_index;
        ca_field_srcptr res_field;

        xfield = x->field;
        field_index = xfield & ~CA_SPECIAL;

        _ca_make_field_element(res, field_index, ctx);
        res_field = CA_FIELD(res, ctx);
        res->field = xfield;  /* set special flags */

        if (field_index == CA_FIELD_ID_QQ)
        {
            fmpq_set(CA_FMPQ(res), CA_FMPQ(x));
        }
        else if (CA_FIELD_IS_NF(res_field))
        {
            nf_elem_set(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_FIELD_NF(res_field));
        }
        else
        {
            fmpz_mpoly_q_set(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_FIELD_MCTX(res_field, ctx));
        }
    }
}

