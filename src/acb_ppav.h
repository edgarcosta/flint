/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_THETA_H
#define ACB_THETA_H

#ifdef __cplusplus
extern "C" {
#endif

#include "fmpz_types.h"
#include "acb_types.h"

/* Periods from thetas in genus 1 and 2 by the AGM method */

int acb_ppav_agm_half_plane(arf_t angle, acb_srcptr a, const arf_t eps,
    slong g, slong prec);
void acb_ppav_agm_max_diff(arf_t delta, acb_srcptr a, slong g, slong prec);

int acb_ppav_agm(acb_ptr r, acb_srcptr a, const arf_t eps,
    slong g, slong prec);

int acb_ppav_periods_from_theta2(acb_mat_t res, acb_srcptr th,
    slong g, slong prec);

/* Periods of elliptic curves */

/* Periods of genus 2 curves */

int acb_ppav_g2_weierstrass(acb_ptr w, const acb_poly_t crv, slong prec);
void acb_ppav_g2_rosenhain(acb_ptr ros, acb_srcptr w, slong perm, slong prec);
slong acb_ppav_g2_theta4(acb_ptr th4, acb_srcptr ros, slong prec);
void acb_ppav_g2_theta2(acb_ptr th2, acb_srcptr th4, acb_srcptr ros,
    slong signs, slong prec);

int acb_ppav_g2_periods_certify(const acb_mat_t res, const ca_poly_t crv, slong prec);

int acb_ppav_g2_periods(acb_mat_t tau, slong * perm, slong * signs,
    acb_ptr w, const acb_poly_t crv, slong prec);
void acb_ppav_g2_periods_with_info(acb_mat_t tau, const slong * perm,
    const slong * signs, acb_srcptr w_low, const acb_poly_t crv, slong prec);

/* Curve equations from invariants */

slong acb_ppav_g2_curve_type(const ca_struct * j);
slong acb_ppav_g2_curve_type(const fmpq * j);

int acb_ppav_curve_from_j(acb_poly_t crv, slong * crv_data,
    acb_srcptr j, slong crv_type, slong g, slong prec);
void acb_ppav_curve_from_j_with_info(acb_poly_t crv, const slong * crv_data,
    acb_srcptr j, slong crv_type, slong g, slong prec);

/* Periods from invariants */

int acb_ppav_periods_from_j(acb_mat_t tau, slong * perm, slong * signs,
    acb_ptr w, slong * crv_data, acb_srcptr j, slong crv_type, slong g, slong prec);
void acb_ppav_periods_from_j_with_info(acb_mat_t tau, const slong * perm,
    const slong * signs, acb_srcptr w_low, const slong * crv_data,
    acb_srcptr j, slong crv_type, slong g, slong prec);

/* Big period matrices and periods in Hilbert space */

int acb_ppav_curve_isomorphism(acb_mat_t r, const acb_poly_t crv1,
    const acb_poly_t crv2, slong g, slong crv_type, slong prec);

int acb_ppav_g2_periods_hilbert(acb_ptr t, const acb_mat_t tau,
    const acb_mat_t curve_isom, const fmpq_mat_t rm_action, slong prec);

/* Hecke operators */

void acb_ppav_g2_siegel_coset();
void acb_ppag_g2_siegel_2step_coset();
void acb_ppav_g2_hilbert_coset();

void acb_ppav_g2_hecke_images(acb_ptr hecke, const acb_mat_t tau,
    const acb_t cofactor, slong ell, slong prec);
slong acb_ppav_g2_hecke_images_integral(fmpz_vec_t * j, slong ** cosets,
    acb_srcptr hecke, )

slong acb_ppav_g2_hecke_int_images_poss(fmpz_vec_t * j, slong ** cosets,
    const acb_mat_t tau, const acb_t cofactor, slong ell, slong prec);
slong acb_ppav_g2_hecke_hilbert_int_images_poss(fmpz_vec_t * j, slong ** cosets,
    acb_srcptr t, const acb_t cofactor, slong prec); /* ideal missing */

#ifdef __cplusplus
}
#endif

#endif
