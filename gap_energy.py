import numpy as np

import mkid_simulator as ms


def main():
    # Loads usuals materials
    materials = ms.load_materials()

    calculator = ms.BCSGapEnergy(materials["aluminium"])

    tau_list = np.linspace(0.01, 1, 200)
    delta = []
    for tau in tau_list:
        delta.append(calculator.evaluate(tau))

    delta_approx = calculator.evaluate_approx(tau_list)

    data = []
    data.append(
        {
            "x": tau_list,
            "y": delta,
            "title": "Gap energy depending of temperature for Aluminium",
            "xlabel": r"$\frac{T}{T_c}$",
            "ylabel": r"$\frac{\Delta(T)}{\Delta_0}$",
            "plot_options": {"label": "Full"},
            "add_plot": [
                {
                    "x": tau_list,
                    "y": delta_approx,
                    "plot_options": {"label": "Approximated", "linestyle": "dashed"},
                }
            ],
        }
    )

    ms.plots(data)


if __name__ == "__main__":
    main()
