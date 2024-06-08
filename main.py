import numpy as np

import mkid_simulator as ms


def main():
    # Loads usuals materials
    materials = ms.load_materials()

    alpha = 0.005
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
            baseline_poly=np.array([1, 0, 0]),
        )
    )
    mkids.append(
        ms.Mkid(
            1e3,
            1e4,
            alpha,
            volume,
            omega_r=omega_r,
            x=evaluated_x,
            init_phase=0.1 * np.pi,
            baseline_poly=np.array([0.1, 5e-3, 1e-3]),
        )
    )

    evaluated_photon = np.array(
        [
            1.3278431947791491e-12,
            1.3276746645087918e-12,
            1.327635679723037e-12,
            1.327573871985013e-12,
            1.3276729954703174e-12,
        ]
    )

    variated_s21 = ms.compute_new_s21(mkids[0], evaluated_x, evaluated_photon)
    print(variated_s21.shape)

    add_plots = [
        {
            "x": mkids[0].x,
            "y": np.abs(sub),
        }
        for i, sub in enumerate(variated_s21)
    ]

    data = []
    data.append(
        {
            "x": mkids[0].x,
            "y": np.abs(mkids[0].s21),
            "title": r"$S_{21}$ avec $\phi_0 = 0.1\pi$ et $g(x) = 0.1 + 0.005x + 0.001x^2$",
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

    ms.add_circular_plot_s21(data, mkids[0])
    ms.add_variation_photon_circular_plot(data, variated_s21)

    ms.plots(data)


if __name__ == "__main__":
    main()
