import numpy as np
import scipy.constants as sc

from scipy.special import jv, kv
from numpy import ndarray

from .mkid import Mkid


def compute_s21_variation(mkid: Mkid, x, power: ndarray):
    print(f"x {x.shape}")

    # Make sure that x and power are arrayLike for loop even with 1 element
    # assert(isinstance(x, ndarray))
    # assert(isinstance(power, ndarray))

    # dx = np.zeros((len(x), len(power)))
    # dqi = np.copy(dx)

    delta_0 = 1.76 * sc.k * 1.2
    n0 = 1.72e10 / 1.6e-19 * 1e18
    tau0 = 438e-9
    temperature = 0.1

    # tau_qp = tau0 / (
    #     np.sqrt(np.pi * temperature / mkid.t_c)
    #     * (2 * delta_0 / sc.k / mkid.t_c) ** (5 / 2)
    #     * np.exp(-delta_0 / sc.k / temperature)
    # )
    # tau_qp = 600e-6
    tau_qp = 1e-3
    # print(f"tau_qp {tau_qp}")

    dx = (
        mkid.alpha
        * s2(2 * np.pi * (x * mkid.omega_r + mkid.omega_r))
        * 0.7
        * tau_qp
        / 4
        / n0
        / delta_0**2
        / mkid.volume
    )
    dx = dx[:, np.newaxis] @ power[np.newaxis, :]
    print(f"dx {dx}")

    dqi = (
        mkid.alpha
        * s1(2 * np.pi * (x * mkid.omega_r + mkid.omega_r))
        * 0.7
        * tau_qp
        / 2
        / n0
        / delta_0**2
        / mkid.volume
    )
    print(f"dqi {dqi}")
    dqi = dqi[:, np.newaxis] @ power[np.newaxis, :]
    return np.array([dx, dqi])


def s2(omega, temperature=0.1):
    delta_0 = 1.76 * sc.k * 1.2
    xsi = sc.hbar * omega / (2 * sc.k * temperature)
    var = 1 + np.sqrt((2 * delta_0) / (np.pi * sc.k * temperature)) * np.exp(-xsi) * jv(
        0, xsi
    )

    return var


def s1(omega, temperature=0.1):
    delta_0 = 1.76 * sc.k * 1.2
    xsi = sc.hbar * omega / (2 * sc.k * temperature)
    var = (
        2
        / np.pi
        * np.sqrt((2 * delta_0) / (np.pi * sc.k * temperature))
        * np.sinh(xsi)
        * kv(0, xsi)
    )

    return var


def compute_new_s21(mkid, evaluated_x, power):

    variations = compute_s21_variation(mkid, evaluated_x, power)
    new_s21 = []
    print(f" var {variations.shape}")

    new_qi = 1 / mkid.q_i
    new_x = mkid.x
    for i, var in enumerate(variations.T):
        # new_qi = mkid.q_i
        # new_qi = 1 / mkid.q_i
        print(f"qi {var[:, 1]}")

        # new_x += var[:, 0] 
        # new_qi += var[:, 1]
        new_qi = 1 / mkid.q_i + var[:, 1]
        new_x = mkid.x + var[:, 0]

        new_qr = 1 / (1 / mkid.q_c + new_qi)
        
        pre_poly = np.sum(
            mkid.baseline_poly
            * np.array(
                # [np.ones(len(mkid.x)), mkid.x + var[:, 0], (mkid.x + var[:, 0]) ** 2]
                [np.ones(len(mkid.x)), new_x, new_x ** 2]
            ).T,
            axis=1,
        )
        new_s21.append(
            pre_poly
            * (
                1
                - (
                    new_qr
                    * np.exp(mkid.init_phase)
                    / mkid.q_c
                    / (1 + 2 * 1j * new_qr * (new_x))
                )
            )
        )
    # print("*")
    # print(len(new_s21))
    # print(new_s21[0].shape)
    all_s21 = np.vstack((mkid.s21, new_s21))
    return all_s21
