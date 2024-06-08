from typing import Union
import numpy as np

import scipy.constants as sc
from scipy.integrate import quad
from scipy.optimize import bisect


from numpy import ndarray

from ..circuit import Material


class BCSGapEnergy:
    """
    Computing Gap Energy for a given material. First computes the eta parameter,
    then with a bisection algorithm optimize the delta parameter.
    For dimensionless purpose, we give:
        -tau = T/T_c
        -delta = Delta(T)/Delta_0

    Methods
    -------
    evaluate: Evaluate the gap energy for a given tau
    evaluate_approx: Evaluate the gap energy for a given tau with the approximate equation
    """

    def __init__(self, m: Material) -> None:
        """
        Attributes
        ----------
        m: Material
            The Material object

        omega_d: float
            Debye frequency

        eta_0: float
            Value of the integral independant of delta's value

        Params
        ------
        m: Material
            Material object for critical temperature, delta_0...
        """
        self.m: Material = m
        self.omega_d: float = m.debye_temperature * sc.k / sc.hbar
        self.eta_0: float = np.arcsinh(self.omega_d * sc.hbar / m.delta_0)

    def evaluate(self, tau: float) -> float:
        """
        Evaluate the gap energy for a given tau

        Params
        ------
        tau: float
            Tau parameter

        Methods
        -------
        eta_function: Integrand
        f: Equation to resolve <int(eta_function) - eta = 0>

        Return
        ------
        float
            The evaluate gap energy
        """

        def eta_function(z: float, delta: float) -> float:
            """
            Integrand to fit delta

            Params
            ------
            z: float
                Integrator

            delta: float
                delta = Delta(T) / Delta_0

            Return
            ------
            eta: float
                value of eta
            """
            temp: float = np.sqrt(1 + z**2)
            eta: float = np.tanh(1.76 / 2 * delta / tau * temp) / temp
            return eta

        def f(delta: float) -> float:
            """
            Optimize the delta parameter for:
                <integrand(delta) - eta_0 = 0

            Params
            ------
            delta: float
                Delta(T) / Delta_0

            Return
            ------
            substract: float
                integrand evaluated at delta minus eta_0
            """
            upper: float = np.sinh(self.eta_0) / delta
            integ: float = quad(eta_function, 0, upper, args=(delta))[0]
            substract: float = self.eta_0 - integ
            return substract

        # Manually sets delta when tau = 0 or >= 1
        if 0 < tau < 1:
            delta_temp = bisect(f, 1e-20, 10)
            assert isinstance(delta_temp, float)
        elif tau == 0:
            delta_temp = 1
        else:
            delta_temp = 0

        return delta_temp

    def evaluate_approx(self, tau: Union[float, ndarray]) -> Union[float, ndarray]:
        """
        Evaluate the gap energy with the approximate equation

        Params
        ------
        tau: float
            T / T_c
        
        Raise
        -----
        ValueError: if tau is 1 because of the division by zero

        Return
        ------
        var: float
            Computed value
        """

        # Check if zero is in tau
        if isinstance(tau, ndarray):
            if len(np.where(tau == 0)[0]) != 0:
                raise ValueError(
                    "Presence of 0 detected! If you want to handle the 0 case,"
                    + " you need to add it manually in your data"
                )

        elif isinstance(tau, float):
            if tau == 1:
                raise ValueError(
                    "Presence of 0 detected! If you want to handle the 0 case,"
                    + " you need to add it manually in your data"
                )
        
        return np.tanh(1.74 * np.sqrt(1 / tau - 1))
