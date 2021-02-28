import numpy as np
import matplotlib.pyplot as plt

tau = 15.5      # TODO [-] Taux de compression == Vmax / Vmin
D = 365e-3      # TODO [m] Alesage == Diamètre de l'alésage piston
C = 420e-3      # TODO [m] Course == 2*Longeur_de_la_manivelle
L = 0.25        # TODO [m] Longueur bielle
mpiston = 5     # TODO [kg] Masse piston
mbielle = 1     # TODO [kg] Masse bielle
Q = 1650e3      # [J/kg_inlet gas] Chaleur emise par fuel par kg de melange admis

Vc = 263.4e-3   # TODO [m^3] Cylindrée == (Vmax - Vmin) * nbre_cylindres
Vmin = Vc / (tau - 1)
Vmax = Vc * tau / (tau - 1)
R = C/2     # 2*Vc / (np.pi * D**2) ou C/2 ?
beta = L/R  # Supposé très grand
gamma = 1.3


def myfunc(rpm, s, thetaC, deltaThetaC):
    """
    :param rpm: la vitesse du moteur [tour/minute] [2500, 4000]
    :param s: le taux de suralimentation [-] (multiplie la pression atmosphérique)
    :param thetaC: l'angle d'allumage (degrés avant PMH) [°] [15°, 40°]
    :param deltaThetaC: la durée de la combustion [°] [40°, 70°]
    :return: t: l'épaisseur critique (section de la bielle)
    """

    # Calcul de la pression
    # h = float(input("Pas de temps [rad] \n>>"))  # pas de temps [rad]
    h = 0.05
    thetaC *= np.pi / 180
    deltaThetaC *= np.pi / 180
    thetaEnd = thetaC + deltaThetaC

    dQ = lambda theta: Q/2 * np.sin(np.pi * (theta - thetaC) / deltaThetaC) * (np.pi / deltaThetaC)
    V  = lambda theta: Vc/2 * (1 - np.cos(theta) + beta - np.sqrt(beta**2 - np.sin(theta)**2)) + Vc / (tau - 1)
    dV = lambda theta: Vc/2 * (np.sin(theta) + np.cos(theta) * np.sin(theta) / np.sqrt(beta**2 - np.sin(theta)**2))
    f  = lambda theta, p: -gamma * p/V(theta) * dV(theta) + (gamma-1) * 1/V(theta) * dQ(theta)

    Theta = np.arange(thetaC, thetaEnd, h)
    P = np.zeros_like(Theta)
    P[0] = 1e5*s
    # Theta, P = Heun(f, thetaC, thetaEnd, 1e5*s, len(Theta))
    for t in range(len(Theta)-1):
        K1 = f(Theta[t], P[t])
        K2 = f(Theta[t] + h / 2, P[t] + h / 2 * K1)
        K3 = f(Theta[t] + h / 2, P[t] + h / 2 * K2)
        K4 = f(Theta[t] + h, P[t] + h * K3)
        P[t+1] = P[t] + h/6 * (K1 + 2*K2 + 2*K3 + K4)

    # Plot
    plt.figure("Pression")
    plt.xlabel("$\\theta [°]$")
    plt.ylabel("$p [bar]$")
    plt.plot(180/np.pi*Theta, P/1e5)
    plt.show()

    # Calcul des efforts
    omega = 2*np.pi*rpm/60  # [rad/s]
    a_bielle = R * omega**2 * np.cos(Theta)
    F_I_pied = mpiston * a_bielle  # Hypothèse: accélération de la bielle == accélération du piston
    F_I_tete = (mpiston + mbielle) * a_bielle
    F_p = np.pi*D**2 / 4 * P
    F_pied = F_p - F_I_pied   # [Connexion avec le piston]
    F_tete = - F_p + F_I_tete # [Connexion avec le vilbrequin]

    Theta *= 180 / np.pi
    # Plot
    plt.figure("Forces")
    plt.xlabel("$\\theta [°]$")
    plt.ylabel("$F [N]$")
    plt.plot(Theta, F_pied, label="$F_{pied} [N]$")
    plt.plot(Theta, F_tete, label="$F_{tete} [N]$")
    plt.plot(Theta, F_p, label="$F_{p} [N]$")
    plt.plot(Theta, F_I_tete, label="$F_{Itete} [N]$")
    plt.plot(Theta, F_I_pied, label="$F_{Ipied} [N]$")
    plt.plot(Theta, a_bielle * mbielle, label="$F_{bielle} [N]$")
    plt.legend()
    plt.show()

    # F_f = [min(F_tete[i], - F_pied[i]) for i in range(len(F_tete))]  # [N]
    # F_tete + F_pied = mbielle * a_bielle
    # acc = np.array([max(abs(F_tete[i]), abs(F_pied[i])) for i in range(len(F_tete))]) - (abs(F_tete) - abs(F_pied))
    F_max = max(abs(F_tete - F_pied))

    # Déterminer la section de la bielle
    # Fc = pi*E / (2*L**2) * S**2
    E = 69e9                              # [Pa]
    alpha = float(input("alpha:\n>>"))
    l = alpha * L
    t = l * np.sqrt(2*F_max / (np.pi*E))  # [m^2]
    return t


print(myfunc(3000, 1, -50, 180))
