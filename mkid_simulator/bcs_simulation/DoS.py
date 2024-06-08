import numpy as np

import scipy.constants as sc
from scipy.integrate import quad

from bcs_gap_energy import BCSGapEnergy
from material import Material


class DensityOfState:
    def __init__(self, temperature: float, m: Material) -> None:
        self.temperature = temperature
        self.m = m

        gap_calculator = BCSGapEnergy(self.m)
        tau = self.temperature / self.m.critical_temperature
        self.delta = gap_calculator.evaluate(tau) * self.m.delta_0 / sc.e
        print(self.delta)

    # def dos(self, energy: float) -> float:
    def dos(self, eps: float) -> float:

        # ns_n0: float = energy / np.sqrt(energy**2 - self.delta**2)
        # print(f"{ns_n0} \t {energy} \t {self.delta}")

        ns_n0: float = eps / np.sqrt(eps**2 - 1)
        # print(f"{ns_n0} \t {1/eps}")

        return ns_n0

    def fermi(self, energy: float):
        if self.temperature == 0:
            return 0
        else:
            print(1 / (1 + np.exp(energy / sc.k / self.temperature)))
            return 1 / (1 + np.exp(energy / sc.k / self.temperature))

    # def particules_density(self, energy: float) -> float:
    def particules_density(self, eps: float) -> float:
        print(eps * self.delta / sc.k * sc.e / 0.1)
        return self.dos(eps) * self.fermi(eps * self.delta * sc.e)
        # return self.dos(energy) * self.fermi(energy)

    def compute_n_qp(self):

        # n_qp = quad(self.particules_density, self.m.fermi_energy + self.delta, np.inf)
        n_qp = quad(
            self.particules_density, self.m.fermi_energy / self.delta + 1, np.inf
        )

        return n_qp[0]
