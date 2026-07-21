/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_PPAV_H
#define ACB_PPAV_H

#ifdef __cplusplus
extern "C" {
#endif

#include "fmpz_types.h"
#include "fmpq_types.h"
#include "acb_types.h"
#include "ca_types.h"

/* Periods from thetas in genus 1 and 2 by the AGM method */

int acb_ppav_agm_half_plane(arf_t beta, acb_srcptr a, const arf_t eps,
    slong g, slong prec);
void acb_ppav_agm_max_diff(arf_t delta, acb_srcptr a, slong g, slong prec);

int acb_ppav_agm(acb_ptr r, acb_srcptr a, const arf_t eps,
    slong g, slong prec);

int acb_ppav_periods_from_theta2(acb_mat_t res, acb_srcptr th,
    slong g, slong prec);

/* Periods of elliptic curves */

/* Periods of genus 2 curves */

void acb_ppav_g2_igusa(ca_vec_t j, const ca_poly_t f);

int acb_ppav_g2_weierstrass(acb_ptr w, const acb_poly_t f, slong prec);
void acb_ppav_g2_rosenhain(acb_ptr ros, acb_srcptr w, slong perm, slong prec);
slong acb_ppav_g2_theta4(acb_ptr th4, acb_srcptr ros, slong prec);
void acb_ppav_g2_theta2(acb_ptr th2, acb_srcptr th4, acb_srcptr ros,
    slong neg, slong signs, slong prec);

int acb_ppav_g2_periods_gather(fmpz_mat_struct * mats,
    const acb_mat_struct * tau_list, slong nb, slong prec);
int acb_ppav_g2_periods_certify(const acb_mat_t r, const acb_poly_t f,
    const ca_vec_t j, slong prec);

typedef struct
{
    ca_struct * j;
    acb_ptr w;
    slong nb_valid_signs;
    slong * perms;
    slong * signs;
    slong * negative_th4s;
    fmpz_mat_struct * mats;
}
acb_ppav_g2_periods_info_struct;

typedef acb_ppav_g2_periods_info_struct acb_ppav_g2_periods_info_t[1];

void acb_ppav_g2_periods_info_init(acb_ppav_g2_periods_info_t info);
void acb_ppav_g2_periods_info_set_nb(acb_ppav_g2_periods_info_t info, slong nb);
void acb_ppav_g2_periods_info_clear(acb_ppav_g2_periods_info_t info);

int acb_ppav_g2_periods_lowprec(acb_mat_t tau, acb_ppav_g2_periods_info_t info,
    const acb_poly_t f, const ca_vec_struct * j, slong prec);
void acb_ppav_g2_periods_highprec(acb_mat_t tau, const acb_ppav_g2_periods_t info,
    const acb_poly_t f, slong prec);

void acb_ppav_g2_periods_acb(acb_mat_t tau, acb_poly_t f, slong prec);
void acb_ppav_g2_periods_ca_j(acb_mat_t tau, acb_poly_t f,
    const ca_vec_struct * j, slong prec);
void acb_ppav_g2_periods_ca(acb_mat_t tau, ca_poly_t f, slong prec);

/* Big period matrices and periods in Hilbert space */

slong acb_ppav_g2_aut(const ca_poly_t crv);

int acb_ppav_g2_isom(acb_mat_t r, const acb_poly_t f1,
    const acb_poly_t f2, slong aut, slong prec);

int acb_ppav_g2_periods_hilbert(acb_ptr t, const acb_mat_t tau,
    const acb_mat_t curve_isom, const ca_mat_t rm_action, slong prec);

/* Hecke operators */

void acb_ppav_g2_siegel_coset();
void acb_ppav_g2_siegel_2step_coset();
void acb_ppav_g2_hilbert_coset();

void acb_ppav_g2_hecke_images(acb_ptr hecke, const acb_mat_t tau,
    const acb_t cofactor, slong ell, slong prec);
slong acb_ppav_g2_hecke_images_integral(fmpz_vec_t * j, slong ** cosets,
    acb_srcptr hecke, slong prec);

slong acb_ppav_g2_hecke_int_images_poss(fmpz_vec_t * j, slong ** cosets,
    const acb_mat_t tau, const acb_t cofactor, slong ell, slong prec);
slong acb_ppav_g2_hecke_hilbert_int_images_poss(fmpz_vec_t * j, slong ** cosets,
    acb_srcptr t, const acb_t cofactor, slong prec); /* ideal missing */

#ifdef __cplusplus
}
#endif

#endif
