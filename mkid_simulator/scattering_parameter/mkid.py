import numpy as np
from numpy import ndarray

from typing import Union


class Mkid:
    def __init__(
        self,
        q_i,
        q_c,
        alpha,
        volume,
        x=np.linspace(-1, 1, 1001),
        delta=0,
        omega_r=0.0,
        t_c=1.2,
        init_phase=0.0,
        baseline_poly=np.array([1, 0, 0]),
    ) -> None:
        # Factory params
        self.q_i = q_i
        self.q_c = q_c
        self.q_r = q_i * q_c / (q_i + q_c)
        self.center = 1 - self.q_r / 2 / q_c
        self.radius = self.q_r / 2 / q_c

        # Frequency to evaluate
        self.omega_r = omega_r
        self.delta = delta
        self.x = x

        self.init_phase = init_phase
        self.baseline_poly = baseline_poly
        self.s21 = self.evaluate_s21()

        if omega_r != 0:
            self.omega = 2 * np.pi * ((self.x * omega_r) + omega_r)

        self.phi = np.arctan(2 * self.q_r * x)

        # Circuit relatives
        self.alpha = alpha
        self.volume = volume
        self.t_c = t_c

    def evaluate_s21(self) -> Union[complex, ndarray]:
        pre_poly = np.sum(
            self.baseline_poly * np.array([np.ones(len(self.x)), self.x, self.x**2]).T,
            axis=1,
        )
        print(self.baseline_poly)
        print(pre_poly)
        s21 = pre_poly * (
            1
            - (
                np.exp(1j * self.init_phase)
                * self.q_r
                / self.q_c
                / (1 + 2 * 1j * self.q_r * (self.x + self.delta))
            )
        )

        return s21
