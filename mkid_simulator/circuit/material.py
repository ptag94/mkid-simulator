from dataclasses import dataclass

import scipy.constants as sc


@dataclass
class Material:
    critical_temperature: float
    debye_temperature: float
    delta_0: float

    # In eV
    fermi_energy: float

    # conductivity in normal state (S/m)
    sigma_n: float


def load_materials():

    data = {
        "niobium": Material(
            critical_temperature=9.5,
            debye_temperature=275,
            delta_0=1.76 * sc.k * 9.5,
            fermi_energy=5.32,
            sigma_n=6.93e6,
        ),
        "aluminium": Material(
            critical_temperature=1.2,
            debye_temperature=420,
            delta_0=1.76 * sc.k * 1.2,
            fermi_energy=11.7,
            sigma_n=37.7e6,
        ),
    }

    return data
