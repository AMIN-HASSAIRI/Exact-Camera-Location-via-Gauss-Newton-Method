"""Starting Prpject"""

import math
import numpy as np

print("Hello Exact Camera Location via Gauss–Newton Method!")
np.set_printoptions(precision=50, suppress=True)

y_0 = [9660, 13105, 4120, -0.040, -0.170, -0.030, -0.072]
print(f"y_0={y_0}")
M = 6
print(f"M={M}")
u = [-0.0480, -0.0100, 0.0490, -0.0190, 0.0600, 0.0125]
print(f"u={u}")
v = [0.029, 0.0305, 0.0285, 0.0115, -0.0005, -0.0270]
print(f"v={v}")
x1 = [9855, 8170, 2885, 8900, 5700, 8980]
print(f"x1={x1}")
x2 = [5680, 5020, 730, 7530, 7025, 11120]
print(f"x2={x2}")
x3 = [3825, 4013, 4107, 3444, 3008, 3412]
print(f"x3={x3}")
epsilon = 10 ** (-3)
print(f"epsilon={epsilon}")
N = 50
print(f"N={N}")

for i in range(0, M):
    print("For the position", i + 1)
    y = y_0
    n = 0
    e = 1
    while (e > epsilon) and (n <= N):
        fi_1 = (
            y[4]
            - (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * (y[3] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[4] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
        ) * (x3[i] - y[2]) - (
            y[5]
            - (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] ** 2 + y[4] ** 2)
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
        ) * (
            x2[i] - y[1]
        )
        fi_2 = (
            y[5]
            - (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] ** 2 + y[4] ** 2)
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
        ) * (x1[i] - y[0]) - (
            y[3]
            + (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * (y[4] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
        ) * (
            x3[i] - y[2]
        )
        fi_3 = (
            y[3]
            + (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * (y[4] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
        ) * (x2[i] - y[1]) - (
            y[4]
            - (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * (y[3] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[4] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
        ) * (
            x1[i] - y[0]
        )
        # ######################### computation of Df(y) #########################

        Df_1_1 = 0

        Df_1_2 = y[5] - (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i]) * (
            (y[3] ** 2 + y[4] ** 2)
            / math.sqrt((y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2))
        )
        Df_1_3 = (
            -y[4]
            + (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * (y[3] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            - (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[4] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
        )
        Df_1_4 = (
            -(math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * (y[4] ** 2 / (y[3] ** 2 + y[4] ** 2) ** (3 / 2))
            - (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] * y[4] * y[5] * (2 * y[3] ** 2 + 2 * y[4] ** 2 + y[5] ** 2))
                / (
                    (y[3] ** 2 + y[4] ** 2) ** (3 / 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (x3[i] - y[2]) + (
            (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] * y[5] ** 2)
                / (
                    math.sqrt(y[3] ** 2 + y[4] ** 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (
            x2[i] - y[1]
        )
        Df_1_5 = (
            (
                1
                + (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
                * ((y[3] * y[4]) / (y[3] ** 2 + y[4] ** 2) ** (3 / 2))
                + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
                * (
                    (y[5] * (y[3] ** 4 - y[4] ** 4 + y[3] ** 2 * y[5] ** 2))
                    / (y[3] ** 2 + y[4] ** 2) ** (3 / 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (x3[i] - y[2]) + (
            (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[4] * y[5] ** 2)
                / (
                    math.sqrt(y[3] ** 2 + y[4] ** 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (
            x2[i] - y[1]
        )
        Df_1_6 = (
            (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[4] * math.sqrt(y[3] ** 2 + y[4] ** 2))
                / (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
            )
        ) * (x3[i] - y[2]) - (
            1
            + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[5] * math.sqrt(y[3] ** 2 + y[4] ** 2))
                / (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
            )
        ) * (
            x2[i] - y[1]
        )
        Df_1_7 = (
            (y[3] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            * (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            + (
                (y[4] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
            * (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
        ) * (x3[i] - y[2]) + (
            (
                (y[3] ** 2 + y[4] ** 2)
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
            * (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
        ) * (
            x2[i] - y[1]
        )
        Df_2_1 = -y[5] + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i]) * (
            (y[3] ** 2 + y[4] ** 2)
            / math.sqrt((y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2))
        )
        Df_2_2 = 0

        Df_2_3 = (
            y[3]
            + (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * (y[4] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
        )
        Df_2_4 = (
            -(math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] * y[5] ** 2)
                / (
                    math.sqrt(y[3] ** 2 + y[4] ** 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (x1[i] - y[0]) - (
            (
                1
                - (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
                * (y[3] * y[4])
                / (y[3] ** 2 + y[4] ** 2) ** (3 / 2)
            )
            + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[5] * (-y[3] ** 4 + y[4] ** 4 + y[4] ** 2 * y[5] ** 2))
                / (
                    (y[3] ** 2 + y[4] ** 2) ** (3 / 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (
            x3[i] - y[2]
        )
        Df_2_5 = (
            -(math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[4] * y[5] ** 2)
                / (
                    math.sqrt(y[3] ** 2 + y[4] ** 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (x1[i] - y[0]) - (
            (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * (y[3] ** 2 / (y[3] ** 2 + y[4] ** 2) ** (3 / 2))
            - (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] * y[4] * y[5] * (2 * y[3] ** 2 + 2 * y[4] ** 2 + y[5] ** 2))
                / (
                    (y[3] ** 2 + y[4] ** 2) ** (3 / 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (
            x3[i] - y[2]
        )
        Df_2_6 = (
            1
            + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[5] * math.sqrt(y[3] ** 2 + y[4] ** 2))
                / (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
            )
        ) * (x1[i] - y[0]) - (
            (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] * math.sqrt(y[3] ** 2 + y[4] ** 2))
                / (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
            )
        ) * (
            x3[i] - y[2]
        )
        Df_2_7 = (
            -(
                (y[3] ** 2 + y[4] ** 2)
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
            * (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
        ) * (x1[i] - y[0]) - (
            (y[4] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            * (-math.sin(y[6]) * u[i] + math.cos(y[6]) * v[i])
            + (
                (y[3] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
            * (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
        ) * (
            x3[i] - y[2]
        )
        Df_3_1 = (
            y[4]
            - (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * (y[3] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[4] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
        )
        Df_3_2 = (
            -y[3]
            - (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * (y[4] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            - (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
        )
        Df_3_3 = 0

        Df_3_4 = (
            1
            - (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * ((y[3] * y[4]) / (y[3] ** 2 + y[4] ** 2) ** (3 / 2))
            + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[5] * (-y[3] ** 4 + y[4] ** 4 + y[4] ** 2 * y[5] ** 2))
                / (
                    (y[3] ** 2 + y[4] ** 2) ** (3 / 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (x2[i] - y[1]) + (
            (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * ((y[4] ** 2) / (y[3] ** 2 + y[4] ** 2) ** (3 / 2))
            - (-math.sin(y[6]) * u[i] + math.cos(y[6]) * v[i])
            * (
                (y[3] * y[4] * y[5] * (2 * y[3] ** 2 + 2 * y[4] ** 2 + y[5] ** 2))
                / (
                    (y[3] ** 2 + y[4] ** 2) ** (3 / 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (
            x1[i] - y[0]
        )
        Df_3_5 = (
            (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * ((y[3] ** 2) / (y[3] ** 2 + y[4] ** 2) ** (3 / 2))
            - (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] * y[4] * y[5] * (2 * y[3] ** 2 + 2 * y[4] ** 2 + y[5] ** 2))
                / (
                    (y[3] ** 2 + y[4] ** 2) ** (3 / 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (x2[i] - y[1]) - (
            1
            + (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
            * ((y[3] * y[4]) / (y[3] ** 2 + y[4] ** 2) ** (3 / 2))
            + (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[5] * (y[3] ** 4 - y[4] ** 4 + y[3] ** 2 * y[5] ** 2))
                / (
                    (y[3] ** 2 + y[4] ** 2) ** (3 / 2)
                    * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
                )
            )
        ) * (
            x1[i] - y[0]
        )
        Df_3_6 = (
            (math.sin(y[6]) * u[i] - math.cos(y[6]) * v[i])
            * (
                (y[3] * math.sqrt(y[3] ** 2 + y[4] ** 2))
                / (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
            )
        ) * (x2[i] - y[1]) + (
            (-math.sin(y[6]) * u[i] + math.cos(y[6]) * v[i])
            * (
                (y[4] * math.sqrt(y[3] ** 2 + y[4] ** 2))
                / (y[3] ** 2 + y[4] ** 2 + y[5] ** 2) ** (3 / 2)
            )
        ) * (
            x1[i] - y[0]
        )
        Df_3_7 = (
            (y[4] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            * (-math.sin(y[6]) * u[i] + math.cos(y[6]) * v[i])
            + (
                (y[3] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
            * (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
        ) * (x2[i] - y[1]) + (
            (y[3] / math.sqrt(y[3] ** 2 + y[4] ** 2))
            * (-math.sin(y[6]) * u[i] + math.cos(y[6]) * v[i])
            - (
                (y[4] * y[5])
                / math.sqrt(
                    (y[3] ** 2 + y[4] ** 2) * (y[3] ** 2 + y[4] ** 2 + y[5] ** 2)
                )
            )
            * (math.cos(y[6]) * u[i] + math.sin(y[6]) * v[i])
        ) * (
            x1[i] - y[0]
        )
        ############################ computation of h ##########################

        Df = np.array(
            [
                [Df_1_1, Df_1_2, Df_1_3, Df_1_4, Df_1_5, Df_1_6, Df_1_7],
                [Df_2_1, Df_2_2, Df_2_3, Df_2_4, Df_2_5, Df_2_6, Df_2_7],
                [Df_3_1, Df_3_2, Df_3_3, Df_3_4, Df_3_5, Df_3_6, Df_3_7],
            ]
        )

        Df_transpose = np.transpose(Df)  # Df_transpose=Df'
        w1 = Df_transpose.dot(Df)
        r = Df_transpose.dot(-1)

        fi = np.array([fi_1, fi_2, fi_3])  # fi_transpose=fi'
        fi_transpose = np.transpose(fi)
        w2 = r.dot(fi_transpose)

        w1_determinant = np.linalg.det(w1)
        if w1_determinant != 0:  # h_column=w1\w2;
            h_column = np.linalg.solve(w1, w2)
        elif w1_determinant == 0:
            solution, residuals, rank, singular_values = np.linalg.lstsq(
                w1, w2, rcond=None
            )
            h_column = solution

        h = np.transpose(h_column)

        Z = np.array(y) + h
        y = Z
        e = np.linalg.norm(h)
        n = n + 1
    print(f"h= {h}")
    print(f"The number of iterations is n= {n}")
    print(f"e= {e}")
    print(f"The solution is y= {y}")
