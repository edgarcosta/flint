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
    compute an angle `\beta` such that `z = e^{-2\pi i\beta} a_j` has real part
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
    signs is definitely wrong, while a (rare) answer of 2 may indicate that the
    chosen precision is too small.

    ALGORITHM: for each `0\leq j\leq 2^g-1`, compute subsets `I_{j,1}` and
    `I_{j,2}` of `\mathbb{R}/\mathbb{Z}` that contain values of `2\pi i \beta`
    that certainly work, resp.~certainly do not work, for this given `j`. These
    subsets are encoded as a finite union of closed intervals in `[0,1]` by
    listing all their endpoints as :type:`arf_t`'s. Let `I_1` be the
    intersection of all `I_{j,1}`'s and let `I_2` be the union of all
    `I_{j,2}`'s. If `I_1` is nonempty, we return 1 and set *beta* to a point in
    it (actually the farthest possible point to the complement of `I_1`);
    otherwise if `I_2` is the whole of `[0,1]` we return 0; otherwise we
    return 2.

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

    ALGORITHM: first call :func:`acb_ppav_agm_half_plane`. If this function
    outputs 0 or 2, we're done. Otherwise rescale *a* by `\exp(-2\pi i \beta)`,
    and start performing AGM steps from this vector, extracting the square
    roots with positive real parts at each step. At each step, compute a
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
      sequences as expected, that the matrix `\tau` corresponding to the output
      *r* is certainly very close to the fundamental domain according to the
      above definition (the choice of sign for `\tau_{1,2}` does not matter
      here), and that its squared theta values computed using
      :func:`acb_theta_sum` to low precision indeed overlap with the input.
    - 2 means we are unsure whether *th* satisfies the specifications or
      not. This indicates the working precision is probably insufficient.

    In cases 0 and 2, the returned value *res* is undefined.

    ALGORITHM: we run the AGM method as in Dupont's PhD thesis,
    Algorithme 13. By [Kie22], we know that Conjecture 9.1 is true, and that
    the AGM sequences can be computed by :func:`acb_ppav_agm` if the tolerance
    parameter *eps* is sufficiently small, say `10^{-10}`.

    The relevant HDME code is in hdme/theta/theta2_inverse.c. In the case where
    the return code is 0, see hdme/theta/theta2_invalid.c and
    hdme/igusa/thomae_discard.c and hdme/igusa/thomae_keep_candidate.c
    (negating the conditions for the latter).

Periods of elliptic curves
-------------------------------------------------------------------------------

TBC. It would make sense to have periods of elliptic curves in FLINT too, and
we will probably need them for the gluing procedure.

Periods of genus 2 curves
-------------------------------------------------------------------------------

.. function:: void acb_ppav_g2_igusa(ca_vec_t j, const ca_poly_t f)

    Computes the fractional Igusa invariants of the Jacobian of `y^2 = f(x)`,
    where *f* has the degree 5 or 6 and is squarefree.

    ALGORITHM: use Ueberschiebungs as in :func:`acb_theta_g2_covariants`.

.. function:: int acb_ppav_g2_weierstrass(acb_ptr w, const acb_poly_t f, slong prec)

    Computes the Weierstrass points of the curves with equation `y^2 =
    f(x)`. It is assumed that *f* has degree exactly 6. A return value of 1
    means that the individual roots have been successfully separated by
    :func:`acb_poly_find_roots`. Otherwise the return value is 0 and *w* is
    left undefined.

    The relevant HDME code is in hdme/igusa/thomae_roots.c.

.. function:: void acb_ppav_g2_rosenhain(acb_ptr ros, acb_srcptr w, slong perm, slong prec)

    Given the Weierstrass points *w* of a genus 2 hyperelliptic curve and a
    permutation `\sigma\in S_6` encoded by the integer *perm*, return the
    Rosenhain invariants associated to the ordered Weierstrass points
    `w_{\sigma(0)},\ldots,w_{\sigma(5)}`. In other words, we bring
    `w_{\sigma(5)}`, `w_{\sigma(0)}` and `w_{\sigma(1)}` to `\infty, 0, 1`
    respectively; the output vector *ros* is then `(w_{\sigma(2)},
    w_{\sigma(3)}, w_{\sigma(4)})`.

    The relevant HDME code is in hdme/igusa/thomae_rosenhain.c and
    hdme/igusa/thomae_reorder.c.

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

    Given vectors *th4*, *ros* as in :func:`acb_ppav_g2_theta4`, extract square
    roots of `\theta_{a,b}^4(0,\tau)` for the characteristics `(a,b)` indexed
    by 1, 2, 4, 8, then use those to determine the remaining values
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

.. function:: int acb_ppav_g2_periods_gather(fmpz_mat_struct * mats, const acb_mat_struct * tau_list, slong nb, slong prec)

    Given a list of *nb* matrices *tau_list* that lie in or very close to the
    fundamental domain, attempt to find symplectic transformations that will
    send them to a common period matrix *tau*. This will be useful in the case
    where the period matrix lies on the boundary of the fundamental domain in
    :func:`acb_ppav_g2_periods_lowprec`; in that case there will be several
    valid sign/permutation choices, but if :func:`acb_ppav_g2_periods_gather`
    succeeds, one can still construct a period matrix that is certainly valid.

    The possible return values are 1 (*tau* was successfully identified) or 0
    (we failed to identify the symplectic transformations and *tau* is left
    undefined).

    ALGORITHM: start from one of the given matrices `\tau_0`, and enumerate the
    symplectic transformations (expressed as a product of transformations
    defining the boundary of the fundamental domain, see
    :func:`sp2gz_fundamental`) that do not certainly decrease `\im(\tau_0)` and
    so that the `L^\infty` norm of the real part remains strictly less
    than 1. There will be only finitely many of those. Then check whether the
    rest of the matrices in *tau_list* overlap with one of the images of
    `\tau_0`. If not, return 0. If yes, say that `\tau_i` overlaps with
    `\gamma_i\tau_0`. Then we let `\tau` be the reunion (in the sense of
    interval arithmetic) of the `\gamma_i^{-1}\tau_i` for each `i`.

    There is no comparable code in HDME (and there should have been).

.. function:: int acb_ppav_g2_periods_certify(const acb_mat_t r, const acb_poly_t f, const ca_vec_t j, slong prec)

    Given a genus 2 curve `y^2 = f(x)` and the Igusa invariants of this curve
    as exact complex numbers, this function attempts to prove that the given
    "ball" matrix *r* indeed contains a matrix of the form

    .. math::

        r = \begin{pmatrix} \tau_{1,1} & \tau_{1,2}^2 \\ \tau_{1,2}^2 & \tau_{2,2} \end{pmatrix}.

    where `\tau` is a legitimate period matrix of the curve. The return value
    is 1 if this can be established with certainty, and 0 otherwise.

    ALGORITHM: There are two cases:

    1. The function which to `(\tau_{1,1}, \tau_{1,2}^2, \tau_{2,2})`
       associates the three (fractional) Igusa invariants is locally a
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
       recognize these theta values as algebraic numbers heuristically; check
       using exact computations that these theta values indeed give the right
       invariants; then use the biholomorphism strategy with the newly found
       exact theta values instead. This will be very expensive but hopefully
       doable.

    There is no comparable code in HDME (and there should have been). This is
    not a priority target for implementation as I suspect/hope it will not be
    needed, as all the special cases where you have several valid
    sign/permutation choices should be covered by
    :func:`acb_ppav_g2_periods_gather`.

.. type:: acb_ppav_g2_periods_info_struct

.. type:: acb_ppav_g2_periods_info_t

    These types are intended to store the information from one run of computing
    periods (typically at low precision) for a given genus 2 curve or from
    given Igusa invariants: choice of curve in Mestre's algorithm, ordering of
    Weierstrass points, possible permutation and sign choices when computing
    squared theta values, and symplectic transformations computed by
    :func:`acb_ppav_g2_periods_gather` if any.

    An :type:`acb_ppav_g2_periods_info_t` is an array of length one of type
    :type:`acb_ppav_g2_periods_info_struct` encoding an ellipsoid as described
    above, alllowing it to be passed by reference.

.. function:: void acb_ppag_g2_periods_info_init(acb_ppav_g2_periods_info_t info, slong nb)

    Initialize *info*.

.. function:: void acb_ppav_g2_periods_info_set_nb(acb_ppav_g2_periods_info_t info, slong nb)

    Sets *info* to contain enough space to contain data on *nb* valid signs and
    permutation choices.

.. function:: void acb_ppav_g2_periods_info_clear(acb_ppav_g2_periods_info_t info)

    Clears *info*.

.. function:: int acb_ppav_g2_periods_lowprec(acb_mat_t tau, acb_ppav_g2_periods_info_t info, const acb_poly_t f, const ca_vec_t j, slong prec)

    Attempts to compute a small period matrix *tau* of the curve `y^2 = f(x)`,
    where *f* has degree 6, and store the relevant information in the *info*
    structure. This function is intended for use at low to moderate
    precisions. If *j* is not NULL, then it should contain the Igusa invariants
    of the genus 2 curve as exact complex numbers; these exact invariants are
    then used to certify the computation if all else fails. The possible return
    values are:

    - 0: The computation failed at some point and return values are undefined,
         most likely due to the working precision being to low.
    - 1: the computed matrix *tau* is certainly correct and *info* contains
         information about all the valid signs and permutation choices.

    ALGORITHM:
    1. Compute Weierstrass points using :func:`acb_ppav_g2_weierstrass`. If
       successful, store them in the *w* field of *info*, otherwise abort. Also
       set the *use_j* field in *info* and, if *use_j* is true, the field *j*.
    2. Loop over permutations, and compute Rosenhain invariants using
       :func:`acb_ppav_g2_rosenhain` and `\theta^4` values for each of them. Store
       the output of :func:`acb_ppav_g2_theta4`.
    3. Loop over sign choices, and compute `\theta^2` values for each of them
       using :func:`acb_ppav_g2_theta2`.
    4. For each candidate tuple `\theta^2`, call
       :func:`acb_ppav_periods_from_theta2` and extract a square root to get
       `\tau_{1,2}` (the sign doesn't matter). If 2 is returned, abort and
       return 0. If 0 is returned, discard this pair of permutation and sign
       choices. If 1 is returned, store this pair as a possible choice, and
       store the output period matrix too.
    5. At this point, there is at least one valid choice of permutation and
       signs. Store all these valid choices in the *info* structure.
    6. If there is exactly 1 valid choice, set *tau* to the computed period
       matrix and exit with return value 1.
    7. If there are at least 2 valid choices, call
       :func:`acb_ppav_g2_periods_gather` on the matrices we obtain (say
       `\tau_0,\ldots,\tau_n`). If 1 is returned, store the symplectic matrices
       (say `\gamma_0,\ldots,\gamma_n`) in *info*, let *tau* be the union (in
       the sense of interval arithmetic) of the matrices `\gamma_j^{-1}\tau_j`
       for `0\leq j\leq n`, and return 1; otherwise continue.
    8. (It would be surprising to enter this step.) If *j* is NULL, return 0,
       otherwise consider each valid pair of permutation and choice of signs in
       order, and call :func:`acb_ppav_g2_periods_certify`. If this function
       returns 1 on some pair, stop, set *info* to contain this unique choice
       of permutation and signs, and return 1; otherwise return 0.

    The relevant HDME code is in hdme/igusa/thomae_correct_signs.c, but that
    code doesn't use exactly the same procedure and was subtly heuristic in
    places. It would probably be better to start from a clean slate.

.. function:: void acb_ppav_g2_periods_highprec(acb_mat_t tau, const acb_ppav_g2_periods_t info, const acb_poly_t f, slong prec)

    Assuming that *info* was successfully set by an earlier call to
    :func:`acb_ppav_g2_periods_lowprec` for (a lower-precision approximation
    of) the exact same curve *f*, run the period computation again using these
    guidelines. This function is intended for use at higher precisions. The
    should not fail on well-formed input. If it does then *tau* is set to an
    infinity value.

    ALGORITHM:
    1. Compute Weierstrass points using :func:`acb_ppav_g2_weierstrass`, and
       match them with the Weierstrass points stored in *info*. Overlaps should
       be 1-to-1; otherwise we abort and return 0. Reorder the Weierstrass
       points to match the ordering from the *info* structure.
    2. Compute Rosenhain invariants, `\theta^4`, then `\theta^2` using the
       values of *perm*, *signs* and *neg* stored in *info*. (There might be
       several valid sign choices, do them all.) This is where the parameter
       *neg* matters: it might be the case that `\theta^4` intersected the
       negative axis at low precision, but doesn't anymore at high
       precision. Still we want to keep the same square root as in the
       low-precision computation.
    3. Call :func:`acb_ppav_periods_from_theta2` on all these tuples and
       extract square roots to get `\tau_{1,2}`. Discard any value for which
       the return code was 0. There will be at least 1 value of *tau* left.
    4. If the number of period matrices we computed is at least 2, use the
       symplectic transformations stored in *info* as in step 7 of
       :func:`acb_ppav_periods_lowprec`. If the results overlap, set *tau* to
       their union and return 1, otherwise continue.
    5. This is the case where we had several valid sign choices at low
       precision, but they somehow don't agree at high precision anymore. This
       is unlikely, but we can still try using
       :func:`acb_ppav_g2_periods_certify` if the parameter *use_j* in info is
       set to true. Otherwise, we output 0.

    In that way we have an algorithm that performs well in practice (as I
    expect step 4 to always succeed) and has theoretical guarantees thanks to
    :func:`acb_ppav_g2_periods_certify`.

.. function:: void acb_ppav_g2_periods_acb(acb_mat_t tau, acb_poly_t f, slong prec)

    Computes the periods of the curve `y^2 = f(x)`.

    ALGORITHM: pick a low starting precision, like 100 bits. Call
   :func:`acb_ppav_g2_periods_lowprec` where *j* is NULL. If 0 is returned,
   increase the precision (up to *prec*) until it works. Once 1 is returned,
   call :func:`acb_ppav_g2_periods_highprec`. If the result is 0, set *tau* to
   an infinite value.

.. function:: void acb_ppav_g2_periods_ca_j(acb_mat_t tau, acb_poly_t f, const ca_vec_t j, slong prec)

    Computes the periods of the curve `y^2 = f(x)`, where we are given in
    addition the Igusa invariants of the curve as exact complex numbers.

    ALGORITHM: Same as :func:`acb_ppav_g2_periods_acb`, except that we use the
    given vector *j* in :func:`acb_ppav_g2_periods_lowprec`.

.. function:: void acb_ppav_g2_periods_ca(acb_mat_t tau, ca_poly_t f, slong prec)

    Computes the periods of the curve `y^2 = f(x)` where *f* is given as an
    exact polynomial.

    ALGORITHM: Compute exact Igusa invariants using :func:`acb_ppav_g2_igusa`
    then call :func:`acb_ppav_g2_periods_ca_j`.

Big period matrices and periods in Hilbert space
-------------------------------------------------------------------------------

.. function:: slong acb_ppav_g2_aut(const ca_poly_t f)

    Compute the automorphism group of the curve `y^2 = f(x)`. There are only
    finitely many possibilities, each corresponding to a possible automorphism
    groups (TBC).

    ALGORITHM: Cardona--Quer.

.. function:: int acb_ppav_g2_isom(acb_mat_t r, const acb_poly_t f1, const acb_poly_t f2, slong aut, slong prec)


