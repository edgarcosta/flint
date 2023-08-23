/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
worker_dim0(acb_ptr dth, const acb_t term, slong* coords, slong g,
    slong ord, slong prec, slong fullprec)
{
    slong n = 1 << g;
    slong nb_max = acb_theta_deriv_nb(ord, g);
    acb_t x;
    fmpz_t m;
    acb_ptr f;
    ulong b;
    slong k, j, i, nb, ind;
    slong* orders;

    acb_init(x);
    fmpz_init(m);
    orders = flint_malloc(g * nb_max * sizeof(slong));
    f = _acb_vec_init(nb_max);

    ind = 0;
    for (k = 0; k <= ord; k++)
    {
        /* Get list of orders */
        nb = acb_theta_deriv_nb(k, g);
        acb_theta_deriv_orders(orders, k, g);

        /* Compute factor for each tuple */
        for (j = 0; j < nb; j++)
        {
            acb_one(&f[j]);
            for (i = 0; i < g; i++)
            {
                fmpz_set_si(m, coords[i]);
                fmpz_pow_ui(m, m, orders[j * g + i]);
                acb_mul_fmpz(&f[j], &f[j], m, prec);
            }
        }
        acb_const_pi(x, prec);
        acb_mul_onei(x, x);
        acb_mul_2exp_si(x, x, 1);
        acb_pow_ui(x, x, ord, prec);
        _acb_vec_scalar_mul(f, f, nb, x, prec);

        /* Loop over b */
        for (b = 0; b < n; b++)
        {
            acb_set(x, term);
            if (acb_theta_char_dot_slong(b, coords, g) % 2 == 1)
            {
                acb_neg(x, x);
            }
            for (j = 0; j < nb; j++)
            {
                acb_addmul(&dth[n * ind + n * j + b], x, &f[j], fullprec);
            }
        }

        ind += nb;
    }

    acb_clear(x);
    fmpz_clear(m);
    flint_free(orders);
    _acb_vec_clear(f, nb_max);
}

static void
acb_theta_naive_0b_jet_gen(acb_ptr dth, slong ord, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    arf_t eps;
    acb_ptr c;
    arb_ptr u;
    acb_ptr new_z;
    slong nb = n * acb_theta_deriv_nb(ord, g + 1);
    slong k;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, nb_z, g);
    arf_init(eps);
    c = _acb_vec_init(nb_z);
    u = _arb_vec_init(nb_z);
    new_z = _acb_vec_init(nb_z * g);

    arf_one(eps);
    arf_mul_2exp_si(eps, eps, -prec);
    acb_theta_naive_ellipsoid(E, new_z, c, u, ord, z, nb_z, tau, eps, prec);
    prec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, new_z, tau, E, prec);

    for (k = 0; k < nb_z; k++)
    {
        acb_theta_naive_worker(&dth[k * nb], nb, &c[k], &u[k], E, D, k,
            ord, prec, worker_dim0);
    }

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    arf_clear(eps);
    _acb_vec_clear(c, nb_z);
    _arb_vec_clear(u, nb_z);
    _acb_vec_clear(new_z, nb_z * g);
}

void
acb_theta_naive_0b_jet(acb_ptr dth, slong ord, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nb = acb_theta_deriv_nb(ord, g + 1);
    acb_ptr res;
    slong k;

    if (g == 1)
    {
        res = _acb_vec_init(4 * nb);
        for (k = 0; k < nb_z; k++)
        {
            acb_modular_theta_jet(res, res + nb, res + 2 * nb, res + 3 * nb,
                z + k * g, acb_mat_entry(tau, 0, 0), nb, prec);
            _acb_vec_set(dth + 2 * k * nb, res + 2 * nb, nb);
            _acb_vec_set(dth + (2 * k + 1) * nb, res + 3 * nb, nb);
        }
        _acb_vec_clear(res, 4 * nb);
    }
    else
    {
        acb_theta_naive_0b_jet_gen(dth, ord, z, nb_z, tau, prec);
    }
}