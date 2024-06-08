import numpy as np

import mkid_simulator as ms


def main():
    # Loads usuals materials
    materials = ms.load_materials()

    e_max = m.delta_0 * 2

    e = np.linspace(0, e_max, 100)

    temperature_list = [0, 1.1, 0.1]
    dos = []
    delta_energy = []
    part_dens = []

    for temp in temperature_list:
        dos_calculator = DensityOfState(temp, m)
        dos.append(dos_calculator.dos(e) / sc.e - m.fermi_energy)

        gap_calculator = BCSGapEnergy(m)
        delta_energy.append(
            gap_calculator.evaluate(temp / m.critical_temperature) * m.delta_0 / sc.e
        )

        part_dens.append(dos_calculator.particules_density(e))

    print(delta_energy)

    fig, ax = plt.subplots(1, 3)
    for ind, i in enumerate(dos):
        ax[0].plot(e / sc.e, i)
        ax[0].axvline(delta_energy[ind])
    ax[0].set_xlim([0, e_max / sc.e])

    ax[1].plot(e / sc.e, part_dens[ind] / sc.e)
    # ax[1].set_xlim([0, e_max / sc.e])
    ax[1].axvline(delta_energy[ind])

    dos_calc = DensityOfState(0.8, m)
    ax[2].plot(e / sc.e, dos_calc.particules_density(e), label="1")
    ax[2].plot(e / sc.e, dos_calc.dos(e), label="2")
    ax[2].set_xlim([0, e_max / sc.e])
    ax[2].axvline(m.fermi_energy)

    plt.show()


if __name__ == "__main__":
    main()
