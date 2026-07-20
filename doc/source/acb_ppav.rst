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
    compute an angle `\beta` such that `z = e^{-2\pi i\beta} a_k` has real part
    at least *eps* times `|z|` for each `0\leq j\leq 2^g - 1`. The real number
    *eps* should be strictly positive.

    The possible outcomes are as follows:

    - 0 means that such an angle `\beta` never exists: in other words, for
      every `\beta`, there exists an index `j` such that `e^{-2\pi i\beta} a_j`
      certainly has real part less than *eps*;
    - 1 means that the computed value of `\beta` certainly works.
    - 2 means that we are unsure about whether `\beta` exists or not. This
      indicates the working precision is probably insufficient.

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
    has to be modified to respect the above specifications. See also
    hdme/theta/borchardt_mean_invalid.c for case 0.

.. function:: void acb_ppav_agm_max_diff(arf_t delta, acb_srcptr a, slong g, slong prec)

    Given a tuple *a* of `2^g` complex numbers, compute an upper bound *delta*
    on the quantity

    .. math::

        \delta = \frac{1}{|a_0|}\max_{0\leq j\leq 2^g - 1} |a_j - a_0|.

    The relevant HDME code is in hdme/theta/borchardt_mean_Delta0.c

.. function:: int acb_ppav_agm(acb_ptr r, acb_srcptr a, const arf_t eps, slong g, slong prec)

    Given a tuple *a* of `2^g` complex numbers, compute an enclosure of its
    arithmetic-geometric mean *r* (the limit of the AGM sequence starting from
    *a* with good sign choices at each step), assuming that
    :func:`acb_ppav_agm_half_plane` returns 1 for the given value of
    *eps*. Otherwise, *r* is left undefined; the return value is the same as
    :func:`acb_ppav_agm_half_plane`.

    ALGORITHM: first call :`acb_ppav_agm_half_plane`. If this function outputs
    0 or 2, we're done. Otherwise rescale *a* by `\exp(-2\pi i \beta)`, and
    start performing AGM steps from this vector, extracting the square roots
    with positive real parts at each step. At each step, compute a
    low-precision approximation of `\delta` as in
    :func:`acb_ppav_agm_max_diff`. Continue until `\delta` is of the order of
    `2^{-\mathrm{prec}}` or of the order of the current error bounds on the
    vector entries. If `\delta` is small enough, say less than 1/8, we obtain
    an enclosure of the AGM value through Proposition 7.1 in Dupont's PhD
    thesis (with k=0), after rescaling by `e^{2\pi i \beta}` Otherwise, we can
    set *r* to a big complex disc centered in 0, since in any case we know that
    the absolute value of the AGM is bounded above by the `L^\infty` norm of
    *a*. (By Dupont, Lemme 7.2, Lemme 7.3 and Lemme 7.4, we know that there
    exists an absolute upper bound, depending on *eps*, on how many steps we
    need to reach the quadratic convergence regime and the number of steps
    after that, but we do not need to compute these bounds.)

    The relevant HDME code is in hdme/theta/borchardt_mean.c. Instead of using
    the borchard_step function, we should use the FLINT function
    :func:`acb_theta_agm_mul`.

.. function:: int acb_ppav_periods_from_theta2(acb_mat_t res, acb_srcptr th, slong g, slong prec)

    Assume that the vector *th* contains the values `\theta_{a,b}(0,\tau)^2`
    for all characteristics `(a,b)` ordered as in :ref:`acb_theta.h
    <acb-theta>`, maybe up to a common nonzero scaling factor, for some matrix
    `\tau\in \mathcal{H}_2` that is in or very close to the fundamental domain.
    Here we define "very close to" rigorously as "the function
    :func:`acb_siegel_is_reduced` returns true with *tol_exp* equal to -20."
    Then this function computes the matrix

    .. math::

        r = \begin{pmatrix} \tau_{1,1} & \tau_{1,2}^2 \\ \tau_{1,2}^2 & \tau_{2,2} \end{pmatrix}.

    The possible outcomes are as follows:

    - 0 means that the input vector *th* is certainly not of the required form.
    - 1 means that the input vector *th* indeed leads to well-behaved AGM
      sequences as expected, that the output matrix `\tau` is certainly very
      close to the fundamental domain according to the above definition, and
      that its squared theta values computed using :func:`acb_theta_sum` to low
      precision indeed overlap with the input.
    - 2 means we are unsure whether *th* satisfies the specifications or
      not. This indicates the working precision is probably insufficient.

    In cases 0 and 2, the returned value *res* is undefined.

    ALGORITHM: we run the AGM method as in Dupont's PhD thesis,
    Algorithme 13. By [Kie22], we know that Conjecture 9.1 is true, and that
    the AGM sequences can be computed by :func:`acb_ppav_agm` if the tolerance
    parameter *eps* is sufficiently small, say `10^{-10}`.

    The relevant HDME code is in hdme/theta/theta2_inverse.c and
    hdme/theta/theta2_invalid.c in case 0.

Periods of elliptic curves
-------------------------------------------------------------------------------

TBC. It would make sense to have periods of elliptic curves in FLINT too, and
we will probably need them for the gluing procedure.

Periods of genus 2 curves
-------------------------------------------------------------------------------

.. function:: int acb_ppav_g2_weierstrass(acb_ptr w, const acb_poly_t crv, slong prec)

    Computes the Weierstrass points of the curves with equation `y^2 = f(x)`,
    where `f` is the polynomial given by *crv*. It is assumed that *crv* has
    degree exactly 6. A return value of 1 means that the individual roots have
    been successfully separated by :func:`acb_poly_find_roots`. Otherwise the
    return value is 0 and *w* is left undefined.

    The relevant HDME code is in hdme/igusa/thomae_roots.c.

.. function:: void acb_ppav_g2_rosenhain(acb_ptr ros, acb_srcptr w, slong perm, slong prec)

    Given the Weierstrass points *w* of a genus 2 hyperelliptic curve and a
    permutation `\sigma\in S_6` encoded by the integer *perm*, return the
    Rosenhain invariants associated to the ordered Weierstrass points
    `w_{\sigma(0)},\ldots,w_{\sigma(5)}`. In other words, we bring
    `w_{\sigma(5)}`, `w_{\sigma(0)}` and `w_{\sigma(1)}` to `\infty, 0, 1`
    respectively; the output vector *ros* is then `(w_{\sigma(2)},
    w_{\sigma(3)}, w_{\sigma(4)})`.

    The relevant HDME code is in hdme/igusa/thomae_rosenhain.c.

.. function:: slong acb_ppav_g2_theta4(acb_ptr th4, acb_srcptr ros, slong prec)

    Given a tuple of Rosenhain invariants, computes the associated values of
    `\theta_{a,b}^4(0,\tau)` for all characteristics `(a,b)` up to a common
    scalar factors using Thomae's formulas. (Only the even characteristics are
    relevant as odd theta constants are identically zero.)

    The return value is a 16-bit integer indicating, for each characteristic,
    whether the computed approximation of `\theta^4` value intersects the
    negative imaginary axis (1) or does not intersect it (0).

    The relevant HDME code is in hdme/igusa/thomae_theta4. We should reorder
    the values to use the ordering from :ref:`acb_theta.h <acb-theta>`
    throughout.

.. function:: void acb_ppav_g2_theta2(acb_ptr th2, acb_srcptr th4, acb_srcptr ros, slong neg, slong signs, slong prec)

    Given vectors *th4*, *ros* as in :func:`acb_theta_ppav_g2_theta4`, extract
    square roots of `\theta_{a,b}^4(0,\tau)` for the characteristics `(a,b)`
    indexed by 1, 2, 4, 8, then use those to determine the remaining values
    `\theta_{a,b}^2(0,\tau)`. The square roots we compute are determined by
    *signs* and *neg*, the second being as output by
    :func:`acb_ppav_g2_theta4`. If the bit number `j` in *neg* is 0, we compute
    the square root with positive real part if the bit number `j` in *signs* is
    0, and its opposite otherwise. If the bit number `j` in *neg* is 1, we
    compute the square root with positive imaginary part if the bit number `j`
    in *signs* is 0, and its opposite otherwise. (We are careful not to lose
    precisions in the square root in that case.)

    The relevant HDME code is in hdme/igusa/thomae_theta2.c and
    hdme/theta/borchardt_root_ui.c.

.. function:: int acb_ppav_g2_periods_certify(const acb_mat_t r, const ca_poly_t crv, slong prec)

    Given a genus 2 curve `y^2 = f(x)` where `f` is the polynomial *crv* with
    exact coefficients, this function attempts to prove that the given "ball"
    matrix *res* indeed contains a matrix of the form

    .. math::

        r = \begin{pmatrix} \tau_{1,1} & \tau_{1,2}^2 \\ \tau_{1,2}^2 & \tau_{2,2} \end{pmatrix}.

    where `\tau` is a legitimate period matrix of the curve. The return value
    is 1 is this can be established with certainty, and 0 otherwise.

    ALGORITHM: There are two cases: This function

    1. The function which to `(\tau_{1,1}, \tau_{1,2}^2, \tau_{2,2})`
       associates the (fractional) Igusa invariants is locally a
       biholomorphism. This can be computed from the curve in Magma using exact
       computations using [KPR25], "Computing isogenies from modular equations
       in genus two", Thm. 3.10, and taking the determinant of that
       matrix. Then we can apply something like [Kie22], "Certified Newton
       schemes for the evaluation of low-genus theta functions", Theorem 2.3,
       or another kind of explicit open image theorem for analytic functions.
       Upper bounds on higher derivatives of theta functions around `\tau` are
       easy to obtain.

    2. That function is *not* a biholomorphism locally. Then things get more
       complicated, but we can still rely on the fact that the function
       `\Theta` which to `(\tau_{1,1}, \tau_{1,2}^2, \tau_{2,2})` associates
       `(\theta_{0,b}(0,\tau/2))_{1\leq b\leq 3}` is well-defined and a
       biholomorphism locally; see [Kie22]. We could proceed as follows: first
       recognize the theta values as algebraic integers of huge degree
       heuristically; check using exact computations that these theta values
       indeed give the right invariants; then use the biholomorphism strategy
       with the newly found exact theta values instead. This will be very
       expensive but hopefully doable.

    There is no comparable code in HDME (and there should have been). This is
    not a priority target for implementation as I suspect/hope it will not be
    needed in most cases.

.. function:: 
