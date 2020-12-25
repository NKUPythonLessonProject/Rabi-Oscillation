from Rabi_description import *
import matplotlib.pyplot as plt


def example1():  # 破坏wa=wc
    wc = 1.0 * 2 * np.pi
    wa = 1.5 * 2 * np.pi
    wr = 0.07 * 2 * np.pi
    H = Get_H_operator(wc, wa, wr, True)

    psi0 = Get_initial_state(0, 1, 0)

    tlist = np.linspace(0, 25, 100)

    c_op_list = Get_c_ops(0, 0, 0)

    output = mesolve(H, psi0, tlist, c_op_list, [a.dag() * a, sm.dag() * sm])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(tlist, output.expect[0], label="Cavity")
    ax.plot(tlist, output.expect[1], label="Atom excited state")
    ax.legend()
    ax.set_xlabel('Time')
    ax.set_ylabel('Occupation probability')
    ax.set_title('wa!=wc')
    plt.savefig("1.png")


def example2():  # 无衰减振荡，画eg能量图，画ec能量图，验证wr
    wc = 1.5 * 2 * np.pi
    wa = 1.5 * 2 * np.pi
    wr = 0.05 * 2 * np.pi
    H = Get_H_operator(wc, wa, wr, True)

    psi0 = Get_initial_state(0, 1, 0)

    tlist = np.linspace(0, 25, 100)

    output1 = mesolve(H, psi0, tlist, [], [Get_N_op_ground(), Get_N_op_excited()])

    output2 = mesolve(H, psi0, tlist, [], [Get_N_op_photon(), Get_N_op_excited()])

    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.plot(tlist, output1.expect[0], label="Atom ground state")
    ax1.plot(tlist, output1.expect[1], label="Atom excited state")
    ax1.legend()
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Occupation probability')
    ax1.set_title('Atom ground state and Atom excited state')
    plt.savefig("2.png")

    plt.clf()

    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(tlist, output2.expect[0], label="Cavity")
    ax2.plot(tlist, output2.expect[1], label="Atom excited state")
    ax2.legend()
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Occupation probability')
    ax2.set_title('Cavity and Atom ground state')
    plt.savefig("3.png")

    omega = Get_omega(tlist, output1.expect[1], [1, wr, 0])
    print(omega / wr)
    print(omega)


def example3():  # RWA影响
    wc = 1.5 * 2 * np.pi
    wa = 1.5 * 2 * np.pi
    wr = 2 * 2 * np.pi
    H1 = Get_H_operator(wc, wa, wr, True)
    H2 = Get_H_operator(wc, wa, wr, False)

    psi0 = Get_initial_state(0, 1, 0)

    tlist = np.linspace(0, 3, 100)

    c_op_list = Get_c_ops(0, 0, 0)

    output1 = mesolve(H1, psi0, tlist, c_op_list, [Get_N_op_photon(), Get_N_op_excited()])
    output2 = mesolve(H2, psi0, tlist, c_op_list, [Get_N_op_photon(), Get_N_op_excited()])

    fig1, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(tlist, output1.expect[0], label="Cavity")
    ax1.plot(tlist, output1.expect[1], label="Atom excited state")
    ax1.legend()
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Occupation probability')
    ax1.set_title('use RWA')
    plt.savefig("4.png")

    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(tlist, output2.expect[0], label="Cavity")
    ax2.plot(tlist, output2.expect[1], label="Atom excited state")
    ax2.legend()
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Occupation probability')
    ax2.set_title('not use RWA')
    plt.savefig("5.png")


def example4():  # 衰减的影响
    wc = 1.5 * 2 * np.pi
    wa = 1.5 * 2 * np.pi
    wr = 0.07 * 2 * np.pi
    kappa = 0.008
    gamma = 0.04
    n_th_a = 0.1
    H = Get_H_operator(wc, wa, wr, True)

    c_op_list = Get_c_ops(kappa, gamma, n_th_a)

    psi0 = Get_initial_state(0, 1, 0)

    tlist = np.linspace(0, 50, 100)

    output = mesolve(H, psi0, tlist, c_op_list, [Get_N_op_photon(), Get_N_op_excited()])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(tlist, output.expect[0], label="Cavity")
    ax.plot(tlist, output.expect[1], label="Atom excited state")
    ax.legend()
    ax.set_xlabel('Time')
    ax.set_ylabel('Occupation probability')
    ax.set_title('consider dissipation')
    plt.savefig("6.png")
