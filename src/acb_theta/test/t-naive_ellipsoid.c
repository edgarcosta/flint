/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("naive_ellipsoid....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: sum of terms on border of ellipsoid must be less than eps */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1; + n_randint(state, 4);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        acb_theta_eld_t E;
        acb_mat_t tau;
        acb_ptr c, z, new_z;
        arb_ptr u;
        arf_t eps;
        acb_t term;
        arb_t abs, sum;
        slong nb_z = 1 + n_randint(state, 4);
        ulong ab = n_randint(state, n) << g;
        slong ord = 0;
        int all = 0;
        slong nb_pts;
        slong* pts;
        slong k, j;

        acb_mat_init(tau, g, g);
        acb_theta_eld_init(E, g, g);
        z = _acb_vec_init(g * nb_z);
        new_z = _acb_vec_init(g * nb_z);
        c = _acb_vec_init(nb_z);
        u = _arb_vec_init(nb_z);
        arf_init(eps);
        acb_init(term);
        arb_init(abs);
        arb_init(sum);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        for (k = 0; k < g * nb_z; k++)
        {
            acb_randtest_precise(&z[k], state, prec, bits);
        }
        arf_one(eps);
        arf_mul_2exp_si(eps, eps, -prec);

        /* Test: sum of terms on the border is less than u */
        acb_theta_naive_ellipsoid(E, new_z, c, u, ab, all, ord, z, nb_z, tau, eps, prec);
        nb_pts = acb_theta_eld_nb_border(E);
        pts = flint_malloc(g * nb_pts * sizeof(slong));
        acb_theta_eld_border(pts, E);

        arb_zero(sum);
        for (j = 0; j < nb_z; j++)
        {
            arb_zero(sum);
            for (k = 0; k < nb_pts; k++)
            {
                acb_theta_naive_term(term, new_z + j * g, tau, pts + k * g, 2 * prec);
                acb_abs(abs, term, 2 * prec);
                arb_add(sum, sum, abs, 2 * prec);
            }
        
            arb_sub(abs, sum, &u[j], 2 * prec);
            
            if (!arb_is_negative(abs))
            {
                flint_printf("FAIL\n");
                flint_printf("sum, eps, bound:\n");
                arb_printd(sum, 10);
                flint_printf("\n");
                arf_printd(eps, 10);
                flint_printf("\n");
                arb_printd(&u[j], 10);
                flint_printf("\ntau:\n");
                acb_mat_printd(tau, 5);
                flint_printf("new_z:\n");
                _acb_vec_printd(new_z + j * g, g, 10);
                flint_printf("\n");
                acb_theta_eld_print(E);
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        acb_theta_eld_clear(E);
        _acb_vec_clear(z, g * nb_z);
        _acb_vec_clear(new_z, g * nb_z);
        _acb_vec_clear(c, nb_z);
        _arb_vec_clear(u, nb_z);
        arf_clear(eps);
        acb_clear(term);
        arb_clear(abs);
        arb_clear(sum);
        flint_free(pts);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}