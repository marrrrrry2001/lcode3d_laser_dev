import numpy as np
import matplotlib.pyplot as plt
import numba

k_0 = 1
L = 1
t_max = 1
t_step = 1
hx = 0.5
hy = 0.5
hxi = 0.5
Lx = 5
Ly = 5
Lxi = -2.5
x = np.arange(0, Lx, hx)
y = np.arange(0, Ly, hy)
xi = np.arange(Lxi, 0, hxi)
t = np.arange(0, t_max, t_step)
nx = len(x)
ny = len(y)
nxi = len(xi)


# print(nx, ny, nxi)


# @numba.njit
def Tridiag_matrix_alg(a, b, c, d):  # TODO: переписать через готовый метод (какой?)
    N = len(a)
    assert len(a) == len(b) == len(c) == len(d), "Error!"

    #   Прямой ход метода исключения Гаусса:
    for i in range(1, N):
        xi = a[i] / b[i - 1]
        print(xi)
        b[i] -= xi * c[i - 1]
        d[i] -= xi * d[i - 1]

    #   Обратный ход метода исключения Гаусса:
    y = np.zeros(N)
    y[-1] = d[-1] / b[-1]

    for i in range(N - 2, -1, -1):
        y[i] = (d[i] - c[i] * y[i + 1]) / b[i]

    return y


def laser0(xi, x, y):  # начальное приближение лазера, в перетяжке только реальная часть импульса
    global rerhs
    xi_0 = -1.4
    a = 0.1
    r0 = Lx / 2
    sigma_xi = 0.2
    # rerhs = []
    # for i in xi:
    #     for j in x:
    #         for k in y:
    #             if - np.pi * sigma_xi + xi_0 < i < np.pi * sigma_xi + xi_0:
    #                 rerhs.append(a * np.sin((i - xi_0) / sigma_xi) * np.exp(-(j ** 2 + k ** 2) / r0 ** 2))
    #             else:
    #                 rerhs.append(0)
    # if xi > - np.pi * sigma_xi + xi_0:
    print(xi)
    for i in xi:
        print(i.T())
        if i < np.pi * sigma_xi + xi_0:
            # rerhs.append(a * np.sin((xi - xi_0) / sigma_xi) * np.exp(-(x ** 2 + y ** 2) / r0 ** 2))
            rerhs = a * np.sin((xi - xi_0) / sigma_xi) * np.exp(-(x ** 2 + y ** 2) / r0 ** 2)
        else:
            rerhs = 0
    return rerhs


def w(xi, x, y):
    wxi, wx, wy = [], [], []
    j, p = x, y
    # print(xi)
    for i in xi:
        wxi.append(laser0(i, j, p))
        for j in x:
            wx.append(laser0(i, j, p))
            for p in y:
                wy.append(laser0(i, j, p))
    # wxi.append(laser0(xi, x, y))
    # wx.append(laser0(xi, x, y))
    # wy.append(laser0(xi, x, y))
    return [wxi, wx, wy]


# print(w(-3, 4, 5))

def RHS(xi, x, y, t):
    xi, x, y = np.meshgrid(xi, x, y)
    rew = w(xi, x, y)
    # print('rew:', rew[0])
    reRHS = []
    for i in range(len(xi) - 2):
        for j in range(len(x) - 1):
            for p in range(len(y) - 1):
                for tau in range(len(t)):
                    # print(i, j, p, tau)
                    n = i + 2
                    l = j + 1
                    k = p + 1
                    reRHS.append(
                        - rew[n][l][k - 1] / hy ** 2 - rew[n][l][k + 1] / hy ** 2 + 8 / (t_step * hxi) *
                        rew[n - 1][l][k] - 6 / (t_step * hxi) * rew[n][l][k] - 2 / (t_step * hxi) * rew[n - 2][l][
                            k] - 1 / hx ** 2 * (rew[n][l - 1][k] + 2 * tau + rew[n][l + 1][k]))
    return reRHS


# Создадим трехдиагональную матрицу для направления xi

reAn1 = 1. + np.zeros(nxi)
reBn1 = 1. + np.zeros(nxi)
reCn1 = -6. * hx ** 2 / (t_step * hxi) + np.zeros(nxi)
reDn1 = hx ** 2 * np.array(RHS(xi, x, y, t))

# imAn1 = np.zeros(nxi)
# imBn1 = np.zeros(nxi)
# imCn1 = 4. * k_0 * hx ** 2 / t_step + np.zeros(nxi)
# imDn1 = np.zeros(nxi)

# reAn2 = 1. + np.zeros(nxi)
# reBn2 = 1. + np.zeros(nxi)
# reCn2 = -6. * hx ** 2 / (t_step * hxi) + np.zeros(nxi)
# reDn2 = hx ** 2 * Tridiag_matrix_alg(reAn1, reBn1, reCn1, reDn1)
# imAn2 = np.zeros(nx)
# imBn2 = np.zeros(nx)
# imCn2 = 4. * k_0 * hx ** 2 / t_step + np.zeros(nx)
# imDn2 = hx ** 2 * Tridiag_matrix_alg(imAn1, imBn1, imCn1, imDn1)

ReAxi = Tridiag_matrix_alg(reAn1, reBn1, reCn1, reDn1)
# ImAxi = Tridiag_matrix_alg(imAn1, imBn1, imCn1, imDn1)
# ImAxi = Tridiag_matrix_alg(imAn2, imBn2, imCn2, imDn2)
# print(ReAxi)
# print(RHS(xi, 0, 0, 0))
# plt.plot(xi, (np.sqrt(RHS(2, 2, 2, t)[1]) ** 2 + (RHS(xi, 2, 2, t)[0]) ** 2))
plt.plot(xi, ReAxi)
plt.plot(xi, w(xi, x, y)[0])
plt.show()
# w(0, 0, 1)
