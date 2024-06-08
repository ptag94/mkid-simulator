import numpy as np

import mkid_simulator as ms


def main():
    # Loads usuals materials
    materials = ms.load_materials()

    alpha = 0.5
    volume = 4.8e-17
    omega_r = 4e9
    evaluated_x = np.linspace(-1, 1, 10001)

    mkids: list = []
    mkids.append(
        ms.Mkid(
            1e3,
            1e4,
            alpha,
            volume,
            omega_r=omega_r,
            x=evaluated_x,
            init_phase=-0,
            baseline_poly=np.array([1, 0, 0])
        )
    )
    # mkids.append(Mkid(3e1, 1e3, alpha, volume, omega_r=omega_r, delta=0.04))

    # evaluated_x = np.array([-0.1, -0.2])
    # evaluated_photon = np.array([1e-23, 5e-23, 9e-23, 1e-21])
    evaluated_photon = np.array([7.2e-11 / 1e3, 7.2e-11 / 1e3, 7.2e-11 / 1e3])
    # evaluated_photon = np.ones(5) * 1e-12

    variated_s21 = ms.compute_new_s21(mkids[0], evaluated_x, evaluated_photon)

    add_plots = [
        {
            "x": mkids[0].x,
            "y": np.abs(sub),
            "plot_options": {"label": f"+{i + 1}ms d'exposition"},
        }
        for i, sub in enumerate(variated_s21)
    ]

    data = []
    data.append(
        {
            "x": mkids[0].x,
            "y": np.abs(mkids[0].s21),
            "title": r"Scattering parameter $S_{21}$",
            "xlabel": r"$x=\frac{f - f_r}{f_r}$",
            "ylabel": r"$|S_{21}|$",
            "ax_options": {"xlim": [-0.025, 0.025]},
            "plot_options": {},
        }
    )
    data.append(
        {
            "x": mkids[0].x,
            "y": np.abs(mkids[0].s21),
            "title": r"Scattering parameter $S_{21}$",
            "xlabel": r"$x=\frac{f - f_r}{f_r}$",
            "ylabel": r"$|S_{21}|$",
            "ax_options": {"xlim": [-0.05, 0.05]},
            "plot_options": {},
            "add_plot": add_plots,
        }
    )

    ms.plots(data)


if __name__ == "__main__":
    main()
