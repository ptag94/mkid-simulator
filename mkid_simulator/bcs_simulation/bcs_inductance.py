from typing import Union

import numpy as np
from numpy import ndarray

import scipy.constants as sc
from scipy.integrate import quad
import scipy.special as ss

from .bcs_gap_energy import BCSGapEnergy
from ..circuit import Material


class BCSInductance:
    """
    Computing the inductance for one or a set of frequency(s) and temperature(s)
    from the Mattis Bardeen integrals. Note that the return conductivity is
    normalized by sigma_n (which can be found in the Material object)
    For dimensionless equations we set:
        - theta = Delta(T) / (h_bar omega)
        - z = E / (h_bar omega)

    Methods
    -------
    fermi: Evaluate the Fermi distribution
    evaluate: Evaluate global conductivity (simga_1 + 1J * sigma_2) / sigma_n
    evaluate_approx: Evaluate the conductivity from the approximate formulas
    evaluate_real_conductivity: Evaluate sigma_1
    evaluate_im_conductivity: Evaluate sigma_2
    """

    def fermi(
        self, energy: float, temperature: Union[float, ndarray]
    ) -> Union[float, ndarray]:
        """
        Computes Fermi distribution

        Params
        ------
        energy: float
            Energy of the system

        temperature: Union[float, ndarray]
            Temperature of the system

        Return
        ------
        Union[float, ndarray]
            Value of the distribution
        """
        return 1 / (1 + np.exp(energy / sc.k / temperature))

    def evaluate_approx(
        self, temperature: Union[float, ndarray], f: Union[float, ndarray], m: Material
    ) -> Union[complex, ndarray]:
        """
        Evaluate the complex conductivity from the approximated Mattis Bardeen
        equations.

        If <temperature> and <f> are ndarray, the shape of the returned array is:
            - first dimension corresponding to the <temperature> dimension
            - second dimension corresponding to the <f> dimension 

        Params
        ------
        temperature: Union[float, ndarray]
            Temperature of the system

        f: Union[float, ndarray]
            Frequency(s) to evaluate

        m: Material
            Material of the superconductor

        Return
        ------
        Union[complex, ndarray]
            The complex conductivity (sigma_1 + 1J * sigma_2)
        """

        # if both are arrays
        if isinstance(temperature, ndarray) and isinstance(f, ndarray):
            assert isinstance(f, ndarray)
            assert isinstance(temperature, ndarray)

            delta_energy: Union[float, ndarray] = sc.hbar * f * 2 * np.pi
            temp: Union[float, ndarray] = (
                0.5 * delta_energy / sc.k / temperature[:, np.newaxis]
            )
        else:
            delta_energy = sc.hbar * f * 2 * np.pi
            temp = 0.5 * delta_energy / sc.k / temperature

        # Compute the gap for each temperature
        gap_calculator = BCSGapEnergy(m)
        if isinstance(temperature, float):
            delta: Union[float, ndarray] = (
                gap_calculator.evaluate(temperature / m.critical_temperature)
                * m.delta_0
            )
        else:
            assert isinstance(temperature, ndarray)

            delta: Union[float, ndarray] = np.zeros(len(temperature))
            for i, t in enumerate(temperature):
                delta[i] = (
                    gap_calculator.evaluate(t / m.critical_temperature) * m.delta_0
                )

        # If both ndarray
        if isinstance(temperature, ndarray) and isinstance(f, ndarray):
            assert isinstance(temperature, ndarray)
            assert isinstance(f, ndarray)
            assert isinstance(delta, ndarray)

            real_conductivity = (
                2
                * delta[:, np.newaxis]
                / delta_energy
                * np.exp(-m.delta_0 / sc.k / temperature[:, np.newaxis])
                * ss.kv(0, temp)
                * 2
                * np.sinh(temp)
            )
            im_conductivity = (
                np.pi
                * delta[:, np.newaxis]
                / delta_energy
                * (
                    1
                    - 2
                    * np.exp(-m.delta_0 / sc.k / temperature[:, np.newaxis])
                    * np.exp(-temp)
                    * ss.jv(0, temp)
                )
            )
        else:
            real_conductivity = (
                2
                * delta
                / delta_energy
                * np.exp(-m.delta_0 / sc.k / temperature)
                * ss.kv(0, temp)
                * 2
                * np.sinh(temp)
            )

            im_conductivity = (
                np.pi
                * delta
                / delta_energy
                * (
                    1
                    - 2
                    * np.exp(-m.delta_0 / sc.k / temperature)
                    * np.exp(-temp)
                    * ss.jv(0, temp)
                )
            )

        return real_conductivity + 1j * im_conductivity

    def evaluate(
        self,
        temperature: Union[float, ndarray],
        frequency: Union[float, ndarray],
        m: Material,
    ) -> Union[complex, ndarray]:
        """
        Evaluate the complexe conductivity from Mathis and Bardeen equations

        Params
        ------
        temperature: Union[float, ndarray]
            Temperature of the superconductor

        frequency: Union[float, ndarray]
            Frequency of the incidence wave

        m: Material
            Material of the superconductor

        Return
        ------
        sigma: Union[complex, ndarray]
            Complex conductivity (simga_1 + 1J * sigma_2)
            if 2D array, the columns are function of frequency and rows of temperature
        """

        # Test if we must compute list or just float
        test: ndarray = np.array(
            [
                isinstance(temperature, float),
                isinstance(frequency, float),
            ]
        )

        # Both float
        if all(test):
            assert isinstance(temperature, float)
            assert isinstance(frequency, float)

            real_conductivity = self.evaluate_real_conductivity(
                temperature, frequency, m
            )
            im_conductivity = self.evaluate_im_conductivity(temperature, frequency, m)
            sigma: Union[complex, ndarray] = real_conductivity + 1j * im_conductivity
            assert isinstance(sigma, complex)

        # Both ndarray
        elif all(~test):
            assert isinstance(temperature, ndarray)
            assert isinstance(frequency, ndarray)

            sigma: Union[complex, ndarray] = np.zeros(
                (len(temperature), len(frequency)), dtype=complex
            )

            for i, temp in enumerate(temperature):
                for j, f in enumerate(frequency):
                    print(f"Evaluating for {temp} K {f/1e9} GHz")
                    real_conductivity = self.evaluate_real_conductivity(temp, f, m)
                    im_conductivity = self.evaluate_im_conductivity(temp, f, m)

                    sigma[i, j] = real_conductivity + 1j * im_conductivity
            assert isinstance(sigma, ndarray)

        # If only temperature is a float
        elif test[0]:
            assert isinstance(temperature, float)
            assert isinstance(frequency, ndarray)

            sigma: Union[complex, ndarray] = np.zeros(len(frequency), dtype=complex)

            for i, f in enumerate(frequency):
                real_conductivity = self.evaluate_real_conductivity(temperature, f, m)
                im_conductivity = self.evaluate_im_conductivity(temperature, f, m)
                sigma[i] = real_conductivity + 1j * im_conductivity
            assert isinstance(sigma, ndarray)

        # If only frequency is a float
        else:
            assert isinstance(temperature, ndarray)
            assert isinstance(frequency, float)

            sigma: Union[complex, ndarray] = np.zeros(len(temperature), dtype=complex)

            for i, temp in enumerate(temperature):
                real_conductivity = self.evaluate_real_conductivity(temp, frequency, m)
                im_conductivity = self.evaluate_im_conductivity(temp, frequency, m)
                sigma[i] = real_conductivity + 1j * im_conductivity
            assert isinstance(sigma, ndarray)

        return sigma

    def evaluate_im_conductivity(
        self, temperature: float, f: float, m: Material
    ) -> float:
        """
        Evaluate sigma_2 / sigma_n.

        Params
        ------
        temperature: float
            Temperature of the system

        f: float
            Frequency to evaluate

        m: Material
            Material of the superconductor

        Methods
        -------
        sigma2: Function used as integrand

        Return
        ------
        cond: float
            Imaginary part of the complex conductivity
        """
        delta_energy: float = sc.hbar * f * 2 * np.pi

        def sigma2(z: float) -> float:
            """
            Integrand to compute sigma_2 / sigma_n

            Params
            ------
            z: float
                Dimensionless parameter z = E / (h_bar omega)

            Return
            ------
            float
                The value of the integrand
            """
            diff_fermi: Union[float, ndarray] = 1 - 2 * self.fermi(
                delta_energy * (1 + z), temperature
            )
            assert isinstance(diff_fermi, float)

            g: float = (
                (z**2 + z + theta**2)
                / np.sqrt(theta**2 - z**2)
                / np.sqrt((z + 1) ** 2 - theta**2)
            )

            return diff_fermi * g

        gap_calculator = BCSGapEnergy(m)
        delta = (
            gap_calculator.evaluate(temperature / m.critical_temperature) * m.delta_0
        )
        theta = delta / delta_energy

        cond = quad(sigma2, max(theta - 1, -theta), theta)[0]

        return cond

    def evaluate_real_conductivity(
        self, temperature: float, f: float, m: Material
    ) -> float:
        """
        Evaluate the real part sigma_1 / sigma_n of the complex conductivity

        Params
        ------
        temperature: float
            Temperature of the system

        f: float
            Frequency to evaluate

        m: Material
            Material of the superconductor

        Methods
        -------
        sigma1_therm: Resistivity from the thermal environment
        sigma1_pair: Resistivity from the pair breaking

        Return
        ------
        cond: float
            sigma_1 / sigma_n
        """
        delta_energy: float = sc.hbar * f * 2 * np.pi

        def sigma1_therm(z: float) -> float:
            """
            Thermal part of the resistivity

            Params
            ------
            z: float
                dimensionless parameter z = E / (h_bar omega)

            Return
            ------
            float
                Value of the integrand
            """
            diff_fermi: Union[float, ndarray] = self.fermi(
                z * delta_energy, temperature
            ) - self.fermi(delta_energy * (z + 1), temperature)
            assert isinstance(diff_fermi, float)

            g: float = (
                (z**2 + z + theta**2)
                / np.sqrt(z**2 - theta**2)
                / np.sqrt((z + 1) ** 2 - theta**2)
            )

            return 2 * diff_fermi * g

        def sigma1_pair(z: float) -> float:
            """
            Pair breaking part of the resistivity

            Params
            ------
            z: float
                dimensionless parameter z = E / (h_bar omega)

            Return
            ------
            float
                Value of the integrand
            """
            diff_fermi: Union[float, ndarray] = 1 - self.fermi(
                delta_energy * (z + 1), temperature
            )
            assert isinstance(diff_fermi, float)

            g: float = (
                (z**2 + z + theta**2)
                / np.sqrt(z**2 - theta**2)
                / np.sqrt((z + 1) ** 2 - theta**2)
            )

            return diff_fermi * g

        gap_calculator = BCSGapEnergy(m)
        delta: float = (
            gap_calculator.evaluate(temperature / m.critical_temperature) * m.delta_0
        )
        theta: float = delta / delta_energy

        # If it can't break Cooper pairs, we don't take the second integrand
        if delta_energy < 2 * delta:
            cond: float = quad(sigma1_therm, theta, np.inf)[0]

        else:
            cond: float = (
                quad(sigma1_therm, theta, np.inf)[0]
                - quad(sigma1_pair, theta - 1, -theta)[0]
            )

        return cond
