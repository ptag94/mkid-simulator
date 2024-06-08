import matplotlib.pyplot as plt
import numpy as np

from ..scattering_parameter import Mkid

# from mkid import Mkid
# from analyse import compute_s21_variation, compute_new_s21, s1, s2


def plots(data_list):
    fig = plt.figure()
    n = len(data_list)  # Nombre de subplots à créer

    for i, data in enumerate(data_list):
        # Calculer le nombre de lignes et de colonnes nécessaires
        cols = int(np.ceil(np.sqrt(n)))
        rows = int(np.ceil(n / cols))

        # Ajouter un subplot
        ax = fig.add_subplot(rows, cols, i + 1)
        match_type(ax, data)

        # Options for the plots
        if "ax_options" in data:
            for option in data["ax_options"]:
                match option:
                    case "scale":
                        if "x" in data["ax_options"][option]:
                            ax.set_xscale("log")
                        if "y" in data["ax_options"][option]:
                            ax.set_yscale("log")

                    case "xlim":
                        ax.set_xlim(data["ax_options"][option])

        if "plot_options" in data:
            if "label" in data["plot_options"]:
                ax.legend()

        if "add_plot" in data:
            for i, sub_data in enumerate(data["add_plot"]):
                match_type(ax, sub_data)

            # Check if any label in add_plots
            check_labels: list[bool] = [
                "label" in sub["plot_options"] if "plot_options" in sub else False
                for sub in data["add_plot"]
            ]
            if any(check_labels):
                if ("label" in data["add_plot"][0]["plot_options"]) or (
                    "label" in data["plot_options"]
                ):
                    ax.legend()

        if "title" in data:
            ax.set_title(data["title"])
        if "xlabel" in data:
            ax.set_xlabel(data["xlabel"])
        if "ylabel" in data:
            ax.set_ylabel(data["ylabel"])

    plt.tight_layout()
    plt.show()


def match_type(ax, data):
    if "type" in data:
        match data["type"]:
            case "plot":
                if "plot_options" in data:
                    ax.plot(data["x"], data["y"], **data["plot_options"])
                else:
                    ax.plot(data["x"], data["y"])
            case "scatter":
                if "plot_options" in data:
                    ax.scatter(data["x"], data["y"], **data["plot_options"])
                else:
                    ax.scatter(data["x"], data["y"])
    else:
        if "plot_options" in data:
            ax.plot(data["x"], data["y"], **data["plot_options"])
        else:
            ax.plot(data["x"], data["y"])


def add_plot_s21(data, mkid):

    data.append(
        {
            "x": mkid.x,
            "y": np.abs(mkid.s21),
            "title": r"$S_{21}$: $Q_i$ = "
            + f"{mkid.q_i}"
            + r"|$Q_c = $"
            + f"{mkid.q_c}",
            "xlabel": r"$x$",
            "ylabel": r"$S_{21}$",
            "type": "plot",
            "add_plot": None,
            "plot_options": {},
        }
    )


def add_circular_plot_s21(data, mkid):
        phi = np.arctan(2 * mkid.q_r * np.linspace(-0.5, 0.5, 10000))

        data.append(
            {
                "x": np.real(mkid.s21),
                "y": np.imag(mkid.s21),
                "title": "Cercle de résonance du mkid ",
                "xlabel": r"$\Re(S_{21})$",
                "ylabel": r"$\Im(S_{21})$",
                "type": "scatter",
                "plot_options": {"s": 0.8},
                "add_plot": [
                    {
                        "x": mkid.center - mkid.radius * np.cos(-2 * phi),
                        "y": -mkid.radius * np.sin(-2 * phi),
                        "type": "plot",
                        "plot_options": {"linewidth": 0.3},
                    },
                    {
                        "x": mkid.center,
                        "y": 0,
                        "type": "scatter",
                        "plot_options": {"s": 0.5, "color": "black"},
                    },
                ],
            }
        )


def add_variation_photon_circular_plot(data, variated_s21):

    add_plot = [
            {
                "x": np.real(s21),
                "y": np.imag(s21),
                "type": "scatter",
                "plot_options": {"s": 0.8},
            }
        for s21 in variated_s21[1:, :]
    ]

    data.append(
        {
            "x": np.real(variated_s21[0]),
            "y": np.imag(variated_s21[0]),
            "title": "Cercle de résonance du mkid ",
            "xlabel": r"$\Re(S_{21})$",
            "ylabel": r"$\Im(S_{21})$",
            "type": "scatter",
            "plot_options": {"s": 0.8},
            "add_plot": add_plot,
        }
    )


# def plot_variation_photon(data, mkid, evaluated_x, new_s21, power=None):

#     add_plot = [
#         (
#             {
#                 "x": evaluated_x,
#                 "y": np.abs(s21),
#                 "type": "plot",
#                 "plot_options": {},
#             }
#             if power is None
#             else {
#                 "x": evaluated_x,
#                 "y": np.abs(s21),
#                 "type": "plot",
#                 "plot_options": {"label": f"{power[i]:.2e} W"},
#             }
#         )
#         for i, s21 in enumerate(new_s21)
#     ]

#     data.append(
#         {
#             "x": mkid.x,
#             "y": np.abs(mkid.s21),
#             "title": r"Variation de $S_{21}$ lors de l'absorption d'un photon",
#             "xlabel": r"$x$",
#             "ylabel": r"$S_{21}$",
#             "type": "plot",
#             "plot_options": {},
#             "add_plot": add_plot,
#         }
#     )


# def plot_test(data, mkids):
#     omega = np.linspace(1e9, 10e9, 100)

#     data.append(
#         {
#             "x": omega,
#             "y": s1(omega),
#             "title": "s",
#             "xlabel": r"$\omega$",
#             "ylabel": r"$s$",
#             "type": "plot",
#             "plot_options": {},
#             "add_plot": [{"x": omega, "y": s2(omega), "type": "plot", "plot_options": {}}],
#         }
#     )
