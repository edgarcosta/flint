/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb.h"
#include "acb_ppav.h"

/* Independent transcription of the factorial mixed-radix permutation
   decode the library performs internally (HDME perm_nb_to_images,
   thomae_reorder.c:4).  Kept separate from the implementation so a decode
   bug there is caught here rather than masked.  Base case n==1 matches the
   safe variant in rosenhain.c (never reaches HDME's out-of-bounds n==0
   write). */
static void
_rosenhain_decode_ref(slong * images, slong nb, slong n)
{
    slong fact = 1, k, q, r;

    if (n <= 1)
    {
        if (n == 1)
            images[0] = 0;
        return;
    }
    for (k = 2; k < n; k++)
        fact *= k;              /* (n-1)! */
    q = nb / fact;
    r = nb % fact;
    images[0] = q;
    _rosenhain_decode_ref(images + 1, r, n - 1);
    for (k = 1; k < n; k++)
        if (images[k] >= q)
            images[k] += 1;
}

static void
_rosenhain_set_vec_si(acb_ptr v, const slong * a, slong n)
{
    slong k;

    for (k = 0; k < n; k++)
        acb_set_si(&v[k], a[k]);
}

/* Tight ball enclosing the exact rational p/q (radius 0 when p/q is
   dyadic, e.g. 5/2, 5, 10, -1, -3; ~1 ulp for 3/5). */
static void
_rosenhain_set_frac_si(acb_t res, slong p, slong q, slong prec)
{
    acb_set_si(res, p);
    acb_div_si(res, res, q, prec);
}

TEST_FUNCTION_START(acb_ppav_g2_rosenhain, state)
{
    const slong precs[3] = {64, 128, 256};
    const slong w_int[6] = {0, 1, 2, 3, 4, 6};
    slong iter;

    /* (B) decode pin, prec-independent: 144 decodes to [1,2,0,3,4,5].
       144 is a non-involution, so an inverse-orientation slip in the
       decode shows up here as a wrong image list. */
    {
        slong images[6];
        const slong expected[6] = {1, 2, 0, 3, 4, 5};
        slong k;

        _rosenhain_decode_ref(images, 144, 6);
        for (k = 0; k < 6; k++)
            if (images[k] != expected[k])
                TEST_FUNCTION_FAIL("(B) decode: perm=144 gives image[%wd]=%wd, "
                    "expected %wd\n", k, images[k], expected[k]);
    }

    for (iter = 0; iter < 3; iter++)
    {
        slong prec = precs[iter];
        slong k;
        acb_ptr w = _acb_vec_init(6);
        acb_ptr ros = _acb_vec_init(3);
        acb_t exact;

        acb_init(exact);
        _rosenhain_set_vec_si(w, w_int, 6);

        /* (A) identity perm=0: images=[0,1,2,3,4,5], new=w; the Mobius
           sending (new[5],new[0],new[1])=(6,0,1) -> (inf,0,1) is
           M(t)=5t/(6-t), so ros=(M(2),M(3),M(4))=(5/2,5,10) exactly.
           These inputs are exact integers, so straight-line field
           arithmetic must keep near-full relative accuracy. */
        acb_ppav_g2_rosenhain(ros, w, 0, prec);
        {
            const slong num[3] = {5, 5, 10};
            const slong den[3] = {2, 1, 1};

            for (k = 0; k < 3; k++)
            {
                _rosenhain_set_frac_si(exact, num[k], den[k], prec);
                if (!acb_contains(&ros[k], exact))
                    TEST_FUNCTION_FAIL("(A) prec=%wd: ros[%wd]=%{acb} does not "
                        "contain %wd/%wd\n", prec, k, &ros[k], num[k], den[k]);
                if (acb_rel_accuracy_bits(&ros[k]) < prec - 30)
                    TEST_FUNCTION_FAIL("(A) prec=%wd: ros[%wd] rel accuracy "
                        "%wd < %wd\n", prec, k,
                        acb_rel_accuracy_bits(&ros[k]), prec - 30);
            }
        }

        /* (B) perm=144 (the decode + orientation case): images=[1,2,0,3,4,5],
           new[images[k]]=w[k] gives new=(2,0,1,3,4,6); the Mobius sending
           (new[5],new[0],new[1])=(6,2,0) -> (inf,0,1) is M(t)=3(t-2)/(t-6),
           so ros=(M(1),M(3),M(4))=(3/5,-1,-3). */
        acb_ppav_g2_rosenhain(ros, w, 144, prec);
        {
            const slong num[3] = {3, -1, -3};
            const slong den[3] = {5, 1, 1};

            /* Floor before the containment: acb_contains passes vacuously on
               an indeterminate or very wide ball, so pin finite + accuracy. */
            for (k = 0; k < 3; k++)
            {
                if (!acb_is_finite(&ros[k]))
                    TEST_FUNCTION_FAIL("(B) prec=%wd: ros[%wd]=%{acb} not "
                        "finite\n", prec, k, &ros[k]);
                if (acb_rel_accuracy_bits(&ros[k]) < prec - 30)
                    TEST_FUNCTION_FAIL("(B) prec=%wd: ros[%wd] rel accuracy "
                        "%wd < %wd\n", prec, k,
                        acb_rel_accuracy_bits(&ros[k]), prec - 30);
            }
            for (k = 0; k < 3; k++)
            {
                _rosenhain_set_frac_si(exact, num[k], den[k], prec);
                if (!acb_contains(&ros[k], exact))
                    TEST_FUNCTION_FAIL("(B) prec=%wd: ros[%wd]=%{acb} does not "
                        "contain %wd/%wd\n", prec, k, &ros[k], num[k], den[k]);
            }
        }

        /* (C) affine invariance: the Rosenhain invariants are unchanged
           under x -> 3x - 1 applied to every Weierstrass point.  Integer
           alpha, beta keep the transformed w exact; the two ball triples
           must overlap entrywise. */
        {
            slong wt_int[6];
            acb_ptr wt = _acb_vec_init(6);
            acb_ptr ros_ref = _acb_vec_init(3);
            acb_ptr ros_aff = _acb_vec_init(3);

            for (k = 0; k < 6; k++)
                wt_int[k] = 3 * w_int[k] - 1;
            _rosenhain_set_vec_si(wt, wt_int, 6);

            acb_ppav_g2_rosenhain(ros_ref, w, 0, prec);
            acb_ppav_g2_rosenhain(ros_aff, wt, 0, prec);
            /* Floor before the overlap: acb_overlaps passes vacuously on an
               indeterminate or very wide ball, so pin finite + accuracy on
               both triples. */
            for (k = 0; k < 3; k++)
            {
                if (!acb_is_finite(&ros_ref[k]))
                    TEST_FUNCTION_FAIL("(C) prec=%wd: ros_ref[%wd]=%{acb} not "
                        "finite\n", prec, k, &ros_ref[k]);
                if (acb_rel_accuracy_bits(&ros_ref[k]) < prec - 30)
                    TEST_FUNCTION_FAIL("(C) prec=%wd: ros_ref[%wd] rel accuracy "
                        "%wd < %wd\n", prec, k,
                        acb_rel_accuracy_bits(&ros_ref[k]), prec - 30);
                if (!acb_is_finite(&ros_aff[k]))
                    TEST_FUNCTION_FAIL("(C) prec=%wd: ros_aff[%wd]=%{acb} not "
                        "finite\n", prec, k, &ros_aff[k]);
                if (acb_rel_accuracy_bits(&ros_aff[k]) < prec - 30)
                    TEST_FUNCTION_FAIL("(C) prec=%wd: ros_aff[%wd] rel accuracy "
                        "%wd < %wd\n", prec, k,
                        acb_rel_accuracy_bits(&ros_aff[k]), prec - 30);
            }
            for (k = 0; k < 3; k++)
                if (!acb_overlaps(&ros_ref[k], &ros_aff[k]))
                    TEST_FUNCTION_FAIL("(C) prec=%wd: affine ros_aff[%wd]=%{acb} "
                        "does not overlap ros[%wd]=%{acb}\n", prec, k,
                        &ros_aff[k], k, &ros_ref[k]);

            _acb_vec_clear(wt, 6);
            _acb_vec_clear(ros_ref, 3);
            _acb_vec_clear(ros_aff, 3);
        }

        acb_clear(exact);
        _acb_vec_clear(w, 6);
        _acb_vec_clear(ros, 3);
    }

    TEST_FUNCTION_END(state);
}
