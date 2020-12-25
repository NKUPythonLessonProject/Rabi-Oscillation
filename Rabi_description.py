from qutip import *
from scipy import optimize
import numpy as np

N = 50
# other ops
a = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))


def Get_H_operator(wc, wa, wr, isRWA):  # 哈密顿量
    if isRWA:
        H = wc * a.dag() * a + wa * sm.dag() * sm + 0.5 * wr * (a.dag() * sm + a * sm.dag())
    else:
        H = wc * a.dag() * a + wa * sm.dag() * sm + 0.5 * wr * (a.dag() + a) * (sm + sm.dag())
    return H


def Get_initial_state(n, ce, cg):  # 初态
    psi0 = tensor(basis(N, n), cg * basis(2, 0) + ce * basis(2, 1))
    return psi0


def Get_c_ops(kappa, gamma, n_th_a):  # 衰减
    c_op_list = []

    rate = kappa * (1 + n_th_a)
    if rate > 0.0:
        c_op_list.append(np.sqrt(rate) * a)

    rate = kappa * n_th_a
    if rate > 0.0:
        c_op_list.append(np.sqrt(rate) * a.dag())

    rate = gamma
    if rate > 0.0:
        c_op_list.append(np.sqrt(rate) * sm)
    return c_op_list


def Get_omega(x, y, p=None):  # 三角拟合
    test_func = lambda r, a1, a2, a3: a1 * np.cos(a2 * r + a3)
    params, params_covariance = optimize.curve_fit(test_func, x, y, p0=p)
    return params[1]


def Get_N_op_photon():  #
    return a.dag() * a


def Get_N_op_excited():
    return sm.dag() * sm


def Get_N_op_ground():
    return sm * sm.dag()
