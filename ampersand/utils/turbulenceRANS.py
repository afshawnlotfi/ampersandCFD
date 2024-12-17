"""
-------------------------------------------------------------------------------
  ***    *     *  ******   *******  ******    *****     ***    *     *  ******   
 *   *   **   **  *     *  *        *     *  *     *   *   *   **    *  *     *  
*     *  * * * *  *     *  *        *     *  *        *     *  * *   *  *     *  
*******  *  *  *  ******   ****     ******    *****   *******  *  *  *  *     *  
*     *  *     *  *        *        *   *          *  *     *  *   * *  *     *  
*     *  *     *  *        *        *    *   *     *  *     *  *    **  *     *  
*     *  *     *  *        *******  *     *   *****   *     *  *     *  ******   
-------------------------------------------------------------------------------
 * AmpersandCFD is a minimalist streamlined OpenFOAM generation tool.
 * Copyright (c) 2024 THAW TAR
 * All rights reserved.
 *
 * This software is licensed under the GNU General Public License version 3 (GPL-3.0).
 * You may obtain a copy of the license at https://www.gnu.org/licenses/gpl-3.0.en.html
 */
"""


# input: U, nu, turbulence intensity, eddy viscosity ratio
# output: k, epsilon, omega
def kEpsilon(U: float, nu: float, I: float = 0.16) -> tuple[float, float, float]:
    k = 1.5 * (U * I) ** 2
    epsilon = 1.5 * k ** 1.5 / (0.09 * nu)
    omega = 0.09 * k / epsilon
    print(f"k: {k}, epsilon: {epsilon}, omega: {omega}")
    return k, epsilon, omega

# calculate turbulent intensity for pipe flow
# input: flow velocity (U), nu, D
def calc_intensity(U: float, nu: float, D: float) -> float:
    # TODO: confirm L is D
    Re = U * D / nu
    I = 0.16 * Re ** (-1. / 8.)
    return I

# calculate turbulent length scale for pipe flow
# input: D
def calc_length_scale(D: float) -> float:
    return 0.07 * D

# calculate turbulent length scale for channel flow
# input: channel width (W), channel depth (H)
def calc_length_scale_channel(W: float, H: float) -> float:
    A = W * H
    P = 2 * (W + H)
    D = 4 * A / P
    l = 0.07 * D
    return l

# calculate turbulent kinetic energy 
# input: U, I
def calc_k(U: float, I: float) -> float:
    return 1.5 * (U * I) ** 2

# calculate turbulent dissipation rate
# input: k, nu
def calc_epsilon(k: float, l: float) -> float:
    Cmu = 0.09
    eps = Cmu ** (3. / 4.) * k ** (3. / 2.) / l
    return eps

# calculate turbulent dissipation rate
# input: k, l
def calc_omega(k: float, l: float) -> float:
    Cmu = 0.09
    omega = Cmu ** (-1. / 4.) * k ** 0.5 / l
    return omega


