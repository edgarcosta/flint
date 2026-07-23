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

.. function:: int acb_ppav_weierstrass(acb_ptr w, const acb_poly_t f, slong g, slong prec)

    Computes the Weierstrass points of the curves with equation `y^2 =
    f(x)`. It is assumed that *f* has degree exactly 6 when `g=2`, and exactly
    3 when `g=1`. A return value of 1 means that the individual roots have been
    successfully separated by :func:`acb_poly_find_roots`. Otherwise the return
    value is 0 and *w* is left undefined.

    The relevant HDME code is in hdme/igusa/thomae_roots.c.

.. type:: acb_ppav_g1_periods_info_struct

.. type:: acb_ppav_g1_periods_info_t

    These types are intended to store the information from one run of computing
    periods (typically at low precision) for a given elliptic curve: ordering
    of Weierstrass points, choice of ordering for the period lattice basis, and
    reduction matrix.

    An :type:`acb_ppav_g1_periods_info_t` is an array of length one of type
    :type:`acb_ppav_g1_periods_info_struct` encoding an ellipsoid as described
    above, alllowing it to be passed by reference.

.. function:: void acb_ppag_g1_periods_info_init(acb_ppav_g1_periods_info_t info, slong nb)

    Initialize *info*.

.. function:: void acb_ppav_g1_periods_info_clear(acb_ppav_g1_periods_info_t info)

    Clears *info*.

.. function:: int acb_ppav_g1_periods_lowprec(acb_t tau, acb_ppav_g1_periods_info_t info, const acb_poly_t f, slong prec)

    Attempts to compute a reduced period *tau* of the elliptic curve `y^2 =
    f(x)`, where *f* has degree 3, and store the relevant information in the
    *info* structure. This function is intended for use at low to moderate
    precisions. The possible return values are:

    - 0: The computation failed at some point and return values are undefined,
         most likely due to the working precision being to low.
    - 1: the computed *tau* is certainly correct and *info* contains
         information about all the valid signs and permutation choices.

    ALGORITHM: we use Cremona, "The complex AGM, periods of elliptic curves
    over C and complex elliptic logarithms", theorem 19.

    1. Compute Weierstrass points using :func:`acb_ppav_weierstrass`.
    2. Find an ordering `(e_1,e_2,e_3)` of the Weierstrass points such that the
       pairs `(e_1 - e_3, e_1 - e_2)` and `(e_2 - e_1, e_2 - e_3)` certainly
       lie in a common half plane (each). In other words, if `e_1,e_2,e_3` are
       close to being aligned, we want to make sure that `e_3` is the one in
       the middle. Set the *w* field of *info* to the ordered triple of
       Weierstrass points.
    3. Compute the AGM sequences with good sign choices starting from
       `(\sqrt{e_1 - e_3}, \sqrt{e_1 - e_2})` and `(\sqrt{e_2 -
       e_1},\sqrt{e_2 - e_2})` (where the sign choices are good; we can
       determine those because of step 2). Let `\pi_1` and `\pi_2` be the
       limits. Then we know that `(\pi_1, \pi_2)` is a basis of the period
       lattice. Let `\tau = \pi_2/\pi_1`. if the sign of `\mathrm{Im}(\tau)`
       cannot be decided, output 0 and abort.
    4. If `\tau` has negative imaginary part, set the *neg* bit of *info* to 1
       and negate `\tau`, otherwise set the *neg* bit to 0.
    5. Reduce `\tau` to the fundamental domain using :func:`acb_siegel_reduce`
       and set the *mat* field of *info* to the corresponding matrix. Set *tau*
       to the resulting period point and return 1.

.. function:: void acb_ppav_g1_periods_big(acb_mat_t pi, const acb_ppav_g1_periods_info_t info, const acb_poly_t f, slong prec)

    Assuming that *info* was successfully set by an earlier call to
    :func:`acb_ppav_g1_periods_lowprec` for (a lower-precision approximation
    of) the exact same curve *f*, this function runs the period computation
    again using these guidelines. This is intended for use at higher precisions
    and should not fail on well-formed input. If it does then *tau* is set to
    an infinity value.

    ALGORITHM:

    1. Compute Weierstrass points using :func:`acb_ppav_weierstrass`, and
       reorder them according to the field *w* in *info*.
    2. Compute the AGM sequences, negate `\tau` if necessary, and apply the
       matrix *mat* as in :func:`acb_ppav_g1_periods_lowprec`.
    3. Call :func:`acb_siegel_reduce` again on the result (at high precision)
       to make sure the result lands close to the fundamental domain.
    4. In addition, keep the original AGM values (before quotienting to get
       `\tau`) and keep track of the action of `\mathrm{SL}_2(\mathbb{Z})` to
       get the big period matrix *pi*.

Periods of genus 2 curves
-------------------------------------------------------------------------------

.. function:: void acb_ppav_g2_igusa(ca_vec_t j, const ca_poly_t f)

    Computes the fractional Igusa invariants of the Jacobian of `y^2 = f(x)`,
    where *f* has the degree 5 or 6 and is squarefree.

    ALGORITHM: use Ueberschiebungs as in :func:`acb_theta_g2_covariants`.

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
    periods (typically at low precision) for a given genus 2 curve: ordering of
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
    1. Compute Weierstrass points using :func:`acb_ppav_weierstrass`. If
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
    of) the exact same curve *f*, this function runs the period computation
    again using these guidelines. This function is intended for use at higher
    precisions and should not fail on well-formed input. If it does then *tau*
    is set to an infinity value.

    ALGORITHM:
    1. Compute Weierstrass points using :func:`acb_ppav_weierstrass`, and
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

    ALGORITHM: Cardona--Quer, "Field of moduli and field of definition for
    curves of genus 2". The invariants `R`, etc. should be computed by
    Ueberschiebungs as in :func:`acb_theta_g2_covariants`.

.. function:: void acb_ppav_g2_isom(acb_mat_t r, const acb_poly_t f1, const acb_poly_t f2, slong aut, slong prec)

    Given two degree 6 polynomials `f_1` and `f_2` such that the genus 2 curves
    `y^2 = f_1(x)` and `y^2 = f_2(x)` are isomorphic, compute an automorphism
    between them in the form of a matrix `r\in \mathrm{GL}_2(\mathbb{C})`.

    ALGORITHM: Follow the hyperellisom function in Pari/GP, which is itself
    following
    https://github.com/JRSijsling/hyperelliptic/blob/main/magma/toolbox/isgl2equiv.m
    We should figure out how this algorithm behaves when the curves are inexact
    approximations of curves with extra automorphisms. I suspect we will need
    to know what the automorphism group is a priori, hence the *aut* argument.

.. function:: void acb_ppav_g2_periods_big_ca(acb_mat_t pi, ca_poly_t f, slong prec)

    Given a squarefree polynomial `f` of degree 5 or 6, compute a big period
    matrix `\Pi` of the genus 2 curve `y^2 = f(x)`, with respect to the
    canonical basis of differential forms `(x dx/y, dx/y)` and some homology
    basis.This matrix will be `2\times 4`, of the form `(r\ r\tau)` where
    `\tau` is a small period matrix very close to the Siegel fundamental
    domain.

    ALGORITHM: first compute a small period matrix `\tau` using
    :func:`acb_ppav_g2_periods_ca`. Then compute a second curve equation `y^2 =
    f_\tau(x)` whose canonical basis of differential forms corresponds to the
    basis `(2\pi i dz_1, 2\pi i dz_2)` on the complex torus attached to
    `\tau`. This curve can be computed using
    :func:`acb_theta_g2_sextic_chi5`. Then compute a matrix `r` giving an
    isomorphism between `y^2 = f_\tau(x)` and `y^2 = f(x)` using
    :func:`acb_ppav_g2_isom`, after computing the automorphism group using
    :func:`acb_ppav_g2_aut` and casting the polynomial to
    :type:`acb_poly_t`. This `r` should be exactly the one appearing in the big
    period matrix (or maybe its inverse/transpose/multiplication by `2\pi
    i`... need to sort out the details here).

.. function:: void acb_ppav_g2_periods_rm_ca(acb_ptr t, const ca_poly_t f, const ca_mat_t rm, slong prec)

    Given a squarefree polynomial `f` of degree 5 or 6 such that the Jacobian
    of the genus 2 curve `y^2 = f(x)` has real multiplication by the full ring
    of integers `\mathbb{Z}_K` in a real quadratic field `K =
    \mathbb{Q}(\sqrt{\Delta})`, compute a pair of complex numbers `t =
    (t_1,t_2)` giving the corresponding moduli point on the Hilbert surface
    `\mathrm{SL}_2(\mathbb{Z}_K)\backslash \mathcal{H}_1^2`. The input matrix
    *rm* should encode the action of `\sqrt{\Delta}` on the canonical basis of
    differential forms `(x dx/y, dx/y)` of the curve.

    ALGORITHM: This uses the description in Birkenhake-Wilhelm, "Humbert
    surfaces and the Kummer plane", 4.1 to 4.5. The first step will be to
    compute the rational representation of `\sqrt{\Delta}` acting on the
    lattice `\tau\mathbb{Z}^g + \mathbb{Z}^g`; this will be possible after we
    compute the big period matrix using :func:`acb_ppav_g2_periods_big_ca`, and
    provides the integral coefficients `a,b,c,d,e` of the Humbert singular
    relation satisfied by `\tau`. Then Proposition 4.5 provides the successive
    steps to move `\tau` to the linear image of `\mathcal{H}_1^2`.

    The relevant HDME code is in hdme/hilbert/hilbert_inverse.c (after the
    coefficients `a,b,c,d,e` are computed).

Principally polarized abelian surfaces over Q
-------------------------------------------------------------------------------

.. type:: acb_ppav_g2_Q_struct

.. type:: acb_ppav_g2_Q_t

    These types are intended to store information on a principally polarized
    abelian surface over Q. This abelian surface could be the Jacobian of a
    genus 2 curve, or a product of two elliptic curves, or the Weil restriction
    of an elliptic curve over a quadratic field. We store information related
    to modular invariants, curve equations, real endomorphisms, and complex periods.

    An :type:`acb_ppav_g2_Q_t` is an array of length one of type
    :type:`acb_ppav_g2_Q_struct` encoding an ellipsoid as described above,
    alllowing it to be passed by reference.

.. function:: void acb_ppav_g2_Q_init(acb_ppav_g2_Q_t A)

    Initializes *A*.

.. function:: void acb_ppav_g2_Q_clear(acb_ppav_g2_Q_t A)

    Clears *A*.

.. function:: void acb_ppav_g2_Q_set_jac(acb_ppav_g2_Q_t A, const fmpz_poly_t f)

    Sets *A* to the Jacobian of the genus 2 curve `y^2 = f(x)`. The polynomial
    *f* must be squarefree of degree 5 or 6.

    ALGORITHM: Compute the modular invariants using Ueberschiebungs; compute
    the Igusa invariants; call :func:`acb_ppav_g2_periods_lowprec` at
    increasing precisions until it succeeds to fill in the periods information,
    and set *prec* to the precision where this computation succeeded.

.. function:: void acb_ppav_g2_Q_set_split(acb_ppav_g2_Q_t A, const fmpz_poly_t e1, const fmpz_poly_t e2)

    Sets *A* to the product `E_1\times E_2`, where `E_i` is the elliptic curves
    `y^2 = e_i(x)` for `i=1,2`. The polynomials *e1* and *e2* should be in
    short Weierstrass form `x^3 + ax + b` with nonzero discriminant.

    ALGORITHM: Compute the modular invariants from the Weierstrass
    coefficients; convert *e1* and *e2* to :type:`ca_poly_t`; call
    :func:`acb_ppav_g1_periods_lowprec` to set the *info1* and *info2* fields.

.. function:: void acb_ppav_g2_Q_set_weil(acb_ppav_g2_Q_t A, const nf_elem_t a, const nf_elem_t b)

    Sets *A* to the Weil restriction of the elliptic curve `y^2 = x^3 + ax +
    b`. The elements `a,b` must be integers in a common quadratic field and the
    discriminant must be nonzero.

    ALGORITHM: Compute the polynomial `x^3 + ax + b` and its conjugate as
    :type:`ca_poly_t`'s, then proceed as in :func:`acb_ppav_g2_Q_set_ell`. We
    could check that it isn't actually isomorphic to a product of elliptic
    curves over `\mathbb{Q}` as well. This should be easily seen from the
    invariants.

.. function:: void acb_ppav_g2_Q_set_rm(acb_ppav_g2_Q_t A, const fmpz_poly_t f, const fmpq_mat_t rm)

    Sets *A* to the Jacobian of the genus 2 curve `y^2 = f(x)`. The polynomial
    *f* must be squarefree of degree 5 or 6. The additional assumption in this
    function is that the curve should have real multiplication by the full ring
    of integers in a real quadratic field `\mathbb{Q}(\sqrt{\Delta})` defined
    over the rationals, and that *rm* must be the matrix of the action of
    `\sqrt{\Delta}` on the canonical basis of differential forms `(x dx/y,
    dx/y)` on the curve.

    ALGORITHM: proceed as in :func:`acb_ppav_g2_Q_set_curve` and fill in the
    *rm* field.

.. function:: int acb_ppav_g2_Q_is_jac(const acb_ppav_g2_Q_t A)

    Returns true iff *A* is a Jacobian.

.. function:: int acb_ppav_g2_Q_is_split(const acb_ppav_g2_Q_t A)

    Returns true iff *A* is the product of two elliptic curves over
    `\mathbb{Q}`.

.. function:: int acb_ppav_g2_Q_is_weil(const acb_ppav_g2_Q_t A)

    Returns true iff *A* is the Weil restriction of an elliptic curve over a
    quadratic extension.

.. function:: slong acb_ppav_g2_Q_has_rm(const acb_ppav_g2_Q_t A)

    Returns positive iff *A* has (known) real multiplication by the full ring
    of integers in a real quadratic field. In that case the return value is the
    discriminant of the RM field. Otherwise returns 0.

.. function:: void acb_ppav_g2_Q_modular_invariants(const acb_ppav_g2_Q_t A)

    Sets *m* to the modular invariants of *A* (a vector of length 4), as
    defined in our previous isogeny classes paper.

.. function:: void acb_ppav_g2_Q_igusa_invariants(const acb_ppav_g2_Q_t A)

    Sets *j* to the Igusa invariants of *A* (a vector of length 3). This
    function throws if *A* is not a Jacobian.

.. function:: void acb_ppav_g2_Q_curve(fmpz_poly_t f, const acb_ppav_g2_Q_t A)

    Sets *f* to an equation of the genus 2 curve whose Jacobian is *A*. This
    function throws if *A* is not a Jacobian.

.. function:: void acb_ppav_g2_Q_elliptic_factor(ca_poly_t e, const acb_ppav_g2_Q_t A, slong k)

    Sets *e* to one of the two elliptic factors of *A*. If *A* is split, this
    will be an integral polynomial, otherwise it will be a polynomial defined
    over a quadratic extension. This function throws if *A* is a Jacobian
    (i.e. is not geometrically split).

.. function:: void acb_ppav_g2_Q_periods(acb_mat_t tau, const acb_ppav_g2_Q_t A, slong prec)

    Computes a small period matrix of *A* that lies very close to the Siegel
    fundamental domain. If *A* is a Jacobian, this calls
    :func:`acb_ppav_g2_periods_highprec` using the precomputed information;
    otherwise, we compute elliptic curve periods using
    :func:`acb_ppav_g1_periods_big`.

.. function:: void acb_ppav_g2_Q_periods_big(acb_mat_t pi, const acb_ppav_g2_Q_t A, slong prec)

    Computes a big period matrix of *A* corresponding to its canonical basis of
    differential forms; the homology basis is chosen so that the attached small
    period matrix lies in or very close to the Siegel fundamental domain.

    In the genus 2 case, this will mostly be a repetition of
    :func:`acb_ppav_g2_periods_big_ca`. Maybe keep just one of them? In genus
    1, this just calls :func:`acb_ppav_g1_periods_big`.

.. function:: void acb_ppav_g2_Q_periods_rm(acb_ptr t, const acb_ppav_g2_Q_t A, slong prec)

    In the genus 2 case, this will mostly be a repetition of
    :func:`acb_ppav_g2_periods_rm_ca`. Maybe keep just one of them?

Hecke operators
-------------------------------------------------------------------------------

.. function:: acb_ppav_g2_siegel_coset_nb(slong ell)

    Returns the number of cosets for the Hecke operator `T(\ell)` for PPAV's of
    dimension 2, which is `(\ell^4 - 1)/(\ell - 1)`.

.. function:: acb_ppav_g2_siegel_2step_coset_nb(slong ell)

    Returns the number of cosets for the Hecke operator `T_1(\ell^2)` for
    PPAV's of dimension 2, which is `\ell (\ell^4 - 1)/(\ell - 1)`.

.. function:: acb_ppav_g2_hilbert_coset_nb(const nf_elem_t beta, slong q)

    In this function, `\beta` denotes a totally positive algebraic integer in a
    real quadratic field `K` such that the ideal `(\beta)` decomposes as
    `\mathfrak{q}\mathfrak{c}^2`, where `\mathfrak{q}` is a product of split
    primes in `K` of total norm *q* that is not divisible by any integer `N >
    1`, and `\mathfrak{c}` has norm coprime to *q*. Note that the ideals
    `\mathfrak{q}` and `\mathfrak{c}` are determined uniquely from the input
    data. We are then interested in the Hecke operator which, to a PPAV `A` of
    dimension 2 with RM by `\mathbb{Z}_K`, associates the quotients `A/G`,
    where `G = A[\mathfrak{c}]\oplus H` and `H` is a maximal isotropic subgroup
    in `A[\mathfrak{q}]` (we could take it to be cyclic, but that would be a
    slightly different Hecke operator). Call this Hecke operator `T(\beta,
    q)`. These Hecke operators will allow us to span the isogeny class of those
    abelian surfaces; see the Arxiv note "Spanning isogeny classes of
    principally polarized abelian surfaces with RM".

    Then this function returns the number of cosets for the Hecke operator
    `T(\beta, q)`.

    ALGORITHM: count the number of possible elements that
    :func:`acb_ppav_g2_hilbert_coset` can return!

.. function:: void acb_ppav_g2_siegel_coset(fmpz_mat_t mat, slong k, slong ell)

    Sets *mat* to the matrix encoding the coset number *k* for the Hecke
    operator `T(\ell)` for PPAV's of dimension 2.

    The relevant HDME code is in hdme/hecke/siegel_coset.c.

.. function:: void acb_ppav_g2_siegel_2step_coset(fmpz_mat_t mat, slong k, slong ell)

    Sets *mat* to the matrix encoding the coset number *k* for the Hecke
    operator `T_1(\ell^2)` for PPAV's of dimension 2.

    The relevant HDME code is in hdme/hecke/siegel_T1_coset.c.

.. function:: void acb_ppav_g2_hilbert_coset(fmpz_mat_t mat, slong k, const nf_elem_t beta, slong q)

    Sets *mat* to the matrix encoding the coset number *k* for the Hecke
    operator `T(\beta, q)` defined as in
    :func:`acb_ppav_g2_hilbert_coset_nb`.

    ALGORITHM: Factor *q* as a product of primes `\ell_1\cdots \ell_r`. This
    corresponds to the decomposition of `\mathfrak{q}` as a product of prime
    ideals. Then the isogeny `A\to A/H` (ignoring polarizations) can be written
    as a composition of cyclic isogenies of degrees `\ell_i`. So we should
    obtain all the possible cosets by multiplying the "usual" cosets for the
    Hecke operator `T(\ell_i)` in the elliptic curves case. (Figure out the
    math in the paper...)

    The relevant HDME code (in the prime case) is in
    hdme/hecke/hilbert_coset.c.

Isogenous abelian varieties
-------------------------------------------------------------------------------

.. function:: slong acb_ppav_g2_Q_siegel_isog(acb_ppav_g2_Q_struct ** B, const acb_ppav_g2_Q_t A, slong ell)

    Sets **B* to the list of all PPAV's of genus 2 over Q that are 1-step
    `\ell`-isogenous to *A*, and return the number of such abelian
    varieties. The return value is the number of such abelian varieties. The
    vector **B* will have to be freed by the user.

    ALGORITHM: this is as in our previous paper.
    1. Compute a period matrix at low precision (this should already have been
       done when setting the *A* structure),
    2. Compute the correct cofactor to use in Hecke enumeration (say from the
       determinant of the `r` part of the big period matrix),
    3. Enumerate the modular invariants of all Hecke images using
       :func:`acb_ppav_g2_siegel_coset` and :func:`acb_theta_g2_even_weight`
       and normalize them, compute an upper bound on the absolute value of
       their product. Figure out which of them could possibly be integral after
       increasing the working precision if necessary.
    4. If none are integral, exit with an empty B;
    5. Pick a higher precision that depends on the product of norms computed in
       step 3; refine the period matrix and recompute the remaining modular
       invariants in step 3 to this much higher precision; conclude that the
       corresponding abelian varieties are indeed defined over Q;
    6. (This is a new step compared to our previous paper) Compute the
       normalized curve equations associated to the isogenous periods. They
       should be rational at the very least, and perhaps we can show they have
       integral coefficients too if properly normalized. If that doesn't work,
       we'll still have reconstructed rational curves with the correct
       invariants, so we can skip Mestre, but we'll need to check that we
       indeed have the correct twist.  Also take into account the fact that we
       might find products of elliptic curves and/or Weil restrictions.
    7. Set the corresponding entries of **B*. Here it's a bit of a waste to
       recompute the periods information since we already have them, but it's
       going to be much less expensive than the rest anyway.

    The relevant HDME code can be found in hdme/hecke/hecke_collect_siegel.c,
    hdme/hecke/hecke_normalize_entry.c, hdme/hecke/hecke_make_integral.c,
    hdme/hecke/hecke_has_integral_precision.c,
    hdme/hecke/hecke_integral_highprec.c, hdme/hecke/hecke_all_isog_Q.c, and
    hdme/modular/siegel_direct_isog_Q.c (we can probably avoid the loop there.)

    It might be worth it to mutualize some of these functions with the 2step
    and Hilbert cases, but I'm not sure how exactly at this point. Maybe
    writing the same code three times isn't horrible?

.. function:: slong acb_ppav_g2_Q_siegel_2step_isog(acb_ppav_g2_Q_struct ** B, const acb_ppav_g2_Q_t A, slong ell);

    Same as :func:`acb_ppav_g2_Q_siegel_isog` but in the 2-step case.

    The relevant HDME code can be found in hdme/hecke/hecke_collect_T1.c, other
    hecke files as in :func:`acb_ppav_g2_Q_siegel_isog`, and
    hdme/modular/siegel_2step_direct_isog_Q.c

.. function:: slong acb_ppav_g2_Q_hilbert_isog(acb_ppav_g2_Q_struct ** B, const acb_ppav_g2_Q_t A, const nf_elem_t beta, slong q, slong hmf_cofactor)

    Same as :func:`acb_ppav_g2_Q_siegel_isog` but in the RM case. We will need
    some information about the corresponding graded algebra of Hilbert modular
    forms, say an integer `M` with the following property: for any even `2\leq
    k\leq D` where `D` is known (depends on `\beta`), if `g` is a HMF with
    integral Fourier coefficients, then `M^k g` is an integral polynomial in
    terms of the pullbacks of modular invariants to the Hilbert surface. This
    integer *M* will be an input to this function as *hmf_cofactor*. Another
    issue is that the pullbacks of modular invariants will only generate the
    symmetric Hilbert modular forms, so we might have to process a given
    `\beta` and its conjugate simultaneously.

    Some relevant HDME code can be found in hdme/hecke/hecke_collect_hilbert.c
    and hdme/hecke/hecke_collect_hilbert_sym.c but the correct rescaling in the
    Hilbert case wasn't worked out there.

