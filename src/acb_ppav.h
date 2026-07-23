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

#ifdef ACB_PPAV_INLINES_C
#define ACB_PPAV_INLINE
#else
#define ACB_PPAV_INLINE static inline
#endif

#include "ulong_extras.h"
#include "fmpz_types.h"
#include "fmpq_types.h"
#include "acb_types.h"
#include "nf_elem.h"
#include "ca_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Periods from thetas in genus 1 and 2 by the AGM method */

int acb_ppav_agm_half_plane(arf_t beta, acb_srcptr a, const arf_t eps,
    slong g, slong prec);
void acb_ppav_agm_max_diff(arf_t delta, acb_srcptr a, slong g, slong prec);

int acb_ppav_agm(acb_ptr r, acb_srcptr a, const arf_t eps,
    slong g, slong prec);

int acb_ppav_periods_from_theta2(acb_mat_t res, acb_srcptr th,
    slong g, slong prec);

/* Periods of elliptic curves */

int acb_ppav_weierstrass(acb_ptr w, const acb_poly_t f, slong g, slong prec);

typedef struct
{
    acb_ptr w;
    slong perm;
    int neg;
    fmpz_mat_struct * mat;
}
acb_ppav_g1_periods_info_struct;

typedef acb_ppav_g1_periods_info_struct acb_ppav_g1_periods_info_t[1];

void acb_ppav_g1_periods_info_init(acb_ppav_g1_periods_info_t info);
void acb_ppav_g1_periods_info_clear(acb_ppav_g1_periods_info_t info);

int acb_ppav_g1_periods_lowprec(acb_t tau, acb_ppav_g1_periods_info_t info,
    const acb_poly_t f, slong prec);
void acb_ppav_g1_periods_big(acb_mat_t pi, const acb_ppav_g1_periods_info_t info,
    const acb_poly_t f, slong prec);

/* Periods of genus 2 curves */

void acb_ppav_g2_igusa(ca_vec_t j, const ca_poly_t f);

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
void acb_ppav_g2_periods_highprec(acb_mat_t tau, const acb_ppav_g2_periods_info_t info,
    const acb_poly_t f, slong prec);

void acb_ppav_g2_periods_acb(acb_mat_t tau, acb_poly_t f, slong prec);
void acb_ppav_g2_periods_ca_j(acb_mat_t tau, acb_poly_t f,
    const ca_vec_struct * j, slong prec);
void acb_ppav_g2_periods_ca(acb_mat_t tau, ca_poly_t f, slong prec);

/* Big period matrices and periods in Hilbert space */

slong acb_ppav_g2_aut(const ca_poly_t crv);

void acb_ppav_g2_isom(acb_mat_t r, const acb_poly_t f1,
    const acb_poly_t f2, slong aut, slong prec);

void acb_ppav_g2_periods_big_ca(acb_mat_t pi, ca_poly_t f, slong prec);
void acb_ppav_g2_periods_rm_ca(acb_ptr t, const ca_poly_t f,
    const ca_mat_t rm, slong prec);

/* Principally polarized abelian surfaces over Q */

typedef struct
{
    fmpz_vec_struct m;
    fmpz_poly_struct f;
    ca_poly_struct e1, e2;
    fmpq_mat_struct rm;
    acb_ppav_g1_periods_info_struct info1, info2;
    acb_ppav_g2_periods_info_struct info;
}
acb_ppav_g2_Q_struct;

typedef acb_ppav_g2_Q_struct acb_ppav_g2_Q_t[1];

void acb_ppav_g2_Q_init(acb_ppav_g2_Q_t A);
void acb_ppav_g2_Q_clear(acb_ppav_g2_Q_t A);

void acb_ppav_g2_Q_set_jac(acb_ppav_g2_Q_t A, const fmpz_poly_t f);
void acb_ppav_g2_Q_set_split(acb_ppav_g2_Q_t A, const fmpz_poly_t e1,
    const fmpz_poly_t e2);
void acb_ppav_g2_Q_set_weil(acb_ppav_g2_Q_t A, const nf_elem_t a, const nf_elem_t b);
void acb_ppav_g2_Q_set_rm(acb_ppav_g2_Q_t A, const fmpz_poly_t f,
    const fmpq_mat_t rm);

int acb_ppav_g2_Q_is_jac(const acb_ppav_g2_Q_t A);
int acb_ppav_g2_Q_is_split(const acb_ppav_g2_Q_t A);
int acb_ppav_g2_Q_is_weil(const acb_ppav_g2_Q_t A);
slong acb_ppav_g2_Q_has_rm(const acb_ppav_g2_Q_t A);

void acb_ppav_g2_Q_modular_invariants(fmpz_vec_t m, const acb_ppav_g2_Q_t A);
void acb_ppav_g2_Q_igusa_invariants(fmpq * j, const acb_ppav_g2_Q_t A);
void acb_ppav_g2_Q_curve(fmpz_poly_t f, const acb_ppav_g2_Q_t A);
void acb_ppav_g2_Q_elliptic_factor(ca_poly_t e, const acb_ppav_g2_Q_t A, slong k);

void acb_ppav_g2_Q_periods(acb_mat_t tau, const acb_ppav_g2_Q_t A, slong prec);
void acb_ppav_g2_Q_periods_big(acb_mat_t pi, const acb_ppav_g2_Q_t A, slong prec);
void acb_ppav_g2_Q_periods_rm(acb_ptr t, const acb_ppav_g2_Q_t A, slong prec);

/* Hecke operators */

ACB_PPAV_INLINE slong
acb_ppav_g2_siegel_coset_nb(slong ell)
{
    return (n_pow(ell,4) - 1) / (ell - 1);
}

ACB_PPAV_INLINE slong
acb_ppav_g2_siegel_2step_coset_nb(slong ell)
{
    return ell * acb_ppav_g2_siegel_coset_nb(ell);
}

slong acb_ppav_g2_hilbert_coset_nb(const nf_elem_t beta, slong q);

void acb_ppav_g2_siegel_coset(fmpz_mat_t mat, slong k, slong ell);
void acb_ppav_g2_siegel_2step_coset(fmpz_mat_t mat, slong k, slong ell);
void acb_ppav_g2_hilbert_coset(fmpz_mat_t mat, slong k, const nf_elem_t beta,
    slong q);

/* Isogenous abelian varieties */

slong acb_ppav_g2_Q_siegel_isog(acb_ppav_g2_Q_struct ** B,
    const acb_ppav_g2_Q_t A, slong ell);
slong acb_ppav_g2_Q_siegel_2step_isog(acb_ppav_g2_Q_struct ** B,
    const acb_ppav_g2_Q_t A, slong ell);
slong acb_ppav_g2_Q_hilbert_isog(acb_ppav_g2_Q_struct ** B,
    const acb_ppav_g2_Q_t A, const nf_elem_t beta, slong q,
    slong hmf_cofactor);

#ifdef __cplusplus
}
#endif

#endif
