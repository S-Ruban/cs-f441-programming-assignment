import matplotlib.pyplot as plt  # matplotlib for plotting graphs
import numpy as np  # numpy for mathematical functions and handling arrays

A = 1.0  # area of the patch (in cm²)
C_m = 1.0  # Membrance Capacitance (in μF)
G_Na_bar = 120  # maximal sodium conductance (in mS/cm²)
G_K_bar = 36  # maximal potassium conductance (in mS/cm²)
G_m = 0.3  # voltage independent "leak" conductance (in mS/cm²)
E_Na = 115  # sodium reverse potential (in mV)
E_K = -12  # potassium reverse potential (in mV)
V_rest = 10.613  # reverse potential (in mV)
dt = 0.1  # timestep of the simulation
t_start = 0  # start time of the simulation (in ms)
t_end = 100  # end time of the simulation (in ms)
p_start = 1  # start of the pulse current (in ms)
p_dur = 99  # duration of the pulse current (in ms)
p_end = p_start + p_dur  # end time of the pulse current (in ms)
I_inj = ((0.0635 * 1e-3) / (900 * np.pi * (1e-4) ** 2)) * A  # current injected (in μA)


def alpha_n(V):
    return (10 - V) / (100 * (np.exp((10 - V) / 10) - 1))


def beta_n(V):
    return 0.125 * np.exp(-V / 80)


def alpha_m(V):
    return (25 - V) / (10 * (np.exp(-((V - 25) / 10)) - 1))


def beta_m(V):
    return 4 * np.exp(-V / 18)


def alpha_h(V):
    return 0.07 * np.exp(-V / 20)


def beta_h(V):
    return 1 / (np.exp((30 - V) / 10) + 1)


def tau_n(V):
    return 1 / (alpha_n(V) + beta_n(V))


def tau_m(V):
    return 1 / (alpha_m(V) + beta_m(V))


def tau_h(V):
    return 1 / (alpha_h(V) + beta_h(V))


def n_inf(V):
    return alpha_n(V) * tau_n(V)


def m_inf(V):
    return alpha_m(V) * tau_m(V)


def h_inf(V):
    return alpha_h(V) * tau_h(V)


t = np.arange(t_start, t_end, dt)
n = np.zeros(len(t))
m = np.zeros(len(t))
h = np.zeros(len(t))
V = np.zeros(len(t))
I = np.zeros(len(t))

m[0] = m_inf(V[0])
h[0] = h_inf(V[0])
n[0] = n_inf(V[0])

for i in range(int(p_start / dt), int(p_end / dt)):
    I[i] = I_inj  # setting the current

for i in range(len(t) - 1):
    n[i + 1] = n_inf(V[i]) + (n[i] - n_inf(V[i])) * np.exp(-dt / tau_n(V[i]))
    m[i + 1] = m_inf(V[i]) + (m[i] - m_inf(V[i])) * np.exp(-dt / tau_m(V[i]))
    h[i + 1] = h_inf(V[i]) + (h[i] - h_inf(V[i])) * np.exp(-dt / tau_h(V[i]))
    G_Na = G_Na_bar * (m[i + 1] ** 3) * h[i + 1]  # sodium channel conductance
    G_K = G_K_bar * (n[i + 1] ** 4)  # potassium channel conductance
    G_eff = G_Na + G_K + G_m  # 1/R = 1/R1 + 1/R2 + 1/R3 = G1 + G2 + G3 = G
    I_Na = G_Na * E_Na  # I = V/R = VG
    I_K = G_K * E_K
    I_leak = G_m * V_rest
    I_eff = I_Na + I_K + I_leak
    V_inf = (I_eff + (I[i] / A)) / G_eff
    tau_V = C_m / G_eff  # τ = RC = C/G
    V[i + 1] = V_inf + (V[i] - V_inf) * np.exp(
        -dt / tau_V
    )  # Using Exponential Euler rule

plt.plot(t, V)  # Plot voltage varying with time
plt.plot(t, I)  # Plot current varying with time
plt.title("Voltage/Current vs Time")
plt.legend(["Voltage", "Current"])
plt.xlabel("Time (ms)")
plt.ylabel("Voltage (mV) / Current (μA)")
plt.annotate(
    "Current density = {:.2f} μA/cm²".format(round(I_inj, 2)),
    xy=(50, 10),
    xytext=(30, I_inj + 5),
)
plt.show()  # Show plot
