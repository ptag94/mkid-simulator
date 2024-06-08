import numpy as np

import mkid_simulator as ms


def main():
    # Loads usuals materials
    materials = ms.load_materials()

    alpha = 1
    volume = 4.8e-17
    omega_r = 4e9
    # evaluated_x = np.linspace(-0.01, 0.01, 5)
    evaluated_x = np.array([0.0])

    mkids: list = []
    mkids.append(
        ms.Mkid(
            1e3,
            1e4,
            alpha,
            volume,
            omega_r=omega_r,
            x=evaluated_x,
            init_phase=0,
            baseline_poly=np.array([1, 0, 0])
        )
    )

    evaluated_photon = np.array([1e-12])

    variated_s21 = ms.compute_new_s21(mkids[0], evaluated_x, evaluated_photon)
    print(variated_s21.shape)
    


if __name__ == "__main__":
    main()
