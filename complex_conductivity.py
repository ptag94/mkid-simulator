import numpy as np

import mkid_simulator as ms


def main():
    # Loads usuals materials
    materials = ms.load_materials()

    temp_list = np.array([0.1, 0.5])
    freq_list = np.linspace(1, 300, 100) * 1e9
    ind_calc = ms.BCSInductance()

    conductivity = ind_calc.evaluate(temp_list, freq_list, materials["aluminium"])
    conductivity_approx = ind_calc.evaluate_approx(
        temp_list, freq_list, materials["aluminium"]
    )

    assert isinstance(conductivity, np.ndarray)
    assert isinstance(conductivity_approx, np.ndarray)

    data = []

    # Add real part (sigma_1)
    data.append(
        {
            "x": freq_list / 1e9,
            "y": np.real(conductivity[0, :]),
            "title": r"Normalized $\sigma_1$ for aluminium",
            "xlabel": r"$f$ (GHz)",
            "ylabel": r"$\frac{\sigma_1}{\sigma_n}$",
            "ax_options": {"scale": "y"},
            "plot_options": {"label": f"Full {temp_list[0]} K"},
            "add_plot": [
                {
                    "x": freq_list / 1e9,
                    "y": np.real(conductivity[1, :]),
                    "plot_options": {"label": f"Full {temp_list[1]} K"},
                },
                {
                    "x": freq_list / 1e9,
                    "y": np.real(conductivity_approx[0, :]),
                    "plot_options": {
                        "label": f"Approximated {temp_list[0]} K",
                        "linestyle": "dashed",
                    },
                },
                {
                    "x": freq_list / 1e9,
                    "y": np.real(conductivity_approx[1, :]),
                    "plot_options": {
                        "label": f"Approximated {temp_list[1]} K",
                        "linestyle": "dashed",
                    },
                },
            ],
        }
    )

    # Add imaginary part (sigma_2)
    data.append(
        {
            "x": freq_list / 1e9,
            "y": np.imag(conductivity[0, :]),
            "title": r"Normalized $\sigma_2$ for aluminium",
            "xlabel": r"$f$ (GHz)",
            "ylabel": r"$\frac{\sigma_1}{\sigma_n}$",
            "ax_options": {"scale": "y"},
            "plot_options": {"label": f"Full {temp_list[0]} K"},
            "add_plot": [
                {
                    "x": freq_list / 1e9,
                    "y": np.imag(conductivity[1, :]),
                    "plot_options": {"label": f"Full {temp_list[1]} K"},
                },
                {
                    "x": freq_list / 1e9,
                    "y": np.imag(conductivity_approx[0, :]),
                    "plot_options": {
                        "label": f"Approximated {temp_list[0]} K",
                        "linestyle": "dashed",
                    },
                },
                {
                    "x": freq_list / 1e9,
                    "y": np.imag(conductivity_approx[1, :]),
                    "plot_options": {
                        "label": f"Approximated {temp_list[1]} K",
                        "linestyle": "dashed",
                    },
                },
            ],
        }
    )

    ms.plots(data)


if __name__ == "__main__":
    main()
