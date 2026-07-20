.. _acb-ppav:

**acb_ppav.h** -- Complex abelian varieties and their moduli
===============================================================================

This module provides functions related to periods of curves of genus 1 and 2
through the AGM method. Some related computations, such as Hecke operators, are
also included.

The numerical functions in this module always compute certified error bounds
unless otherwise stated.

The AGM method for periods from theta functions in dimensions 1 and 2
-------------------------------------------------------------------------------

.. function:: int acb_ppav_agm_half_plane(arf_t beta, acb_srcptr a, const arf_t eps, slong g, slong prec)

    Given a tuple *a* of `2^g` complex numbers, this function attemps to
    compute an angle `\beta` such that `e^{-2\pi i\beta} a_k` has real part at
    least *eps* for each `0\leq j\leq 2^g - 1`. The real number *eps* should be
    strictly positive.

    The possible outcomes are as follows:

    - 0 means that such an angle `\beta` never exists: in other words, for
      every `\beta`, there exists an index `j` such that `e^{-2\pi i\beta} a_j`
      certainly has real part less than *eps*;
    - 1 means that the computed value of `\beta` certainly works.
    - 2 means that we are unsure about whether `\beta` exists or not.

    In cases 0 and 2, the output value of *beta* is left undefined.

    USAGE: this function is used to rule out bad orderings of the Weierstrass
    points and bad sign choices in period computations. By [Kie22] ("Sign
    choices in the AGM for genus two theta constants"), we know that this
    function will succeed for a very small *eps* (say `10^{-10}`) when picking
    the correct theta values. An answer of 0 indicates that this choice of
    signs is definitely wrong, while a (rare) answer of 1 may indicate that the
    chosen precision is too small.

    ALGORITHM: for each `0\leq j\leq 2^g-1`, compute subsets `I_{j,1}` of
    `I_{j,1}` of `\mathbb{R}/\mathbb{Z}` that contain values of `2\pi i \beta`
    that certainly work, resp.~certainly do not work, for this given `j`. These
    subsets are encoded as a finite union of closed intervals in `[0,1]` by
    listing all their endpoints as :type:`arf_t`'s. Let `I_1` be the
    intersection of all `I_{j,1}`'s and let `I_2` be the reunion of all
    `I_{j,2}`'s. If `I_1` is nonempty, we return a point in it (actually the
    farthest possible point to the complement of `I_1`); otherwise if `I_2` is
    the whole of `[0,1]` we return 0; otherwise we return 1.

    The relevant HDME code is in hdme/theta/borchardt_excl_half_planes.c but
    has to be modified to respect the above specifications.


