/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_all_sqr(acb_ptr th2, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    fmpz_mat_t mat;
    acb_mat_t w;
    acb_ptr x, aux;
    slong k2;

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(w, g, g);
    x = _acb_vec_init(g);
    aux = _acb_vec_init(n * n);

    acb_siegel_reduce(w, mat, tau, prec);
    acb_siegel_transform_ext(x, w, mat, z, tau, prec);
    acb_theta_ql_all_sqr(aux, x, w, prec);

    sp2gz_inv(mat, mat);
    k2 = acb_theta_transform_k2(mat);
    acb_theta_transform_sqr(th2, aux, x, w, mat, k2, prec);

    fmpz_mat_clear(mat);
    acb_mat_clear(w);
    _acb_vec_clear(x, g);
    _acb_vec_clear(aux, n * n);
}