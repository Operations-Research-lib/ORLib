import math as m
import numpy as np


# -------------------------M/M/1-------------------------------------------------
def mm1_model_info(lam, miu):
    """Calculates the basic information of a M/M/1 queueing system.
    param lam: Arrival rate
    param miu: Service rate
    return: system_info = [rho, Lq, L, Wq, W]"""
    system_info = np.zeros(5)
    system_info[0] = lam / miu
    system_info[1] = pow(lam, 2) / (miu * (miu - lam))
    system_info[2] = lam / (miu - lam)
    system_info[3] = lam / (miu * (miu - lam))
    system_info[4] = 1 / (miu - lam)
    print("rho (Utilization Factor) = ", system_info[0])
    print("Lq (Excpected Lenght of Queue) = ", system_info[1])
    print("L (Excepected number of clientes in the system) = ", system_info[2])
    print("Wq (Expected waiting time in queue) = ", system_info[3])
    print("W (Expected waiting time in system) = ", system_info[4])
    return system_info


def mm1_model_compute_Pn(lam, miu, n):
    """Calculates the probability of n clients in a M/M/1 queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param n: Number of clients
    return: Pn = Probability of n clients"""
    rho = lam / miu
    Pn = 0
    if n == 0:
        Pn = 1 - rho
    elif n > 0:
        Pn = (1 - rho) * pow(rho, n)
    print(f"Probability of {n} clients= ", Pn)
    return Pn


# -------------------------M/M/s-------------------------------------------------
def mms_model_compute_Pzero(lam, miu, s):
    """Calculates the probability of zero clients in a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    return: P_0 = Probability of zero clients"""
    P_0 = 0
    rho = lam / (s * miu)
    for i in range(s):
        P_0 += pow(rho, i) / m.factorial(i)
        P_0 += (pow(rho, s) / (m.factorial(s))) * (1 / (1 - rho))
    P_0 = 1 / P_0
    return P_0


def mms_model_info(lam, miu, s):
    """Calculates the basic information of a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    return: system_info = [rho, Lq, L, Wq, W]"""
    # calculate rho to determine if the system is stable
    rho = lam / (s * miu)
    if rho >= 1:
        print("The system is unstable")
        return None
    # calculate the basic information
    system_info = np.zeros(6)
    system_info[0] = rho
    # calculate probability of zero clients
    system_info[1] = mms_model_compute_Pzero(lam, miu, s)
    # calculate expected lenght of queue
    system_info[2] = (system_info[1] * pow(rho, s)) / (m.factorial(s) * pow(1 - rho, 2))
    # calculate the expected waiting time in queue
    system_info[3] = system_info[2] / lam
    # calculate the expected waiting time in system
    system_info[4] = system_info[3] + (1 / miu)
    # calculate the expected number of clients in the system
    system_info[5] = system_info[2] + (lam / miu)
    print("rho (Utilization Factor) = ", system_info[0])
    print("P_0 (Probability of zero clients) = ", system_info[1])
    print("Lq (Excpected Lenght of Queue) = ", system_info[2])
    print("Wq (Expected waiting time in queue) = ", system_info[3])
    print("W (Expected waiting time in system) = ", system_info[4])
    print("L (Excepected number of clientes in the system) = ", system_info[5])
    return system_info


def mms_model_compute_Pn(lam, miu, s, n):
    """Calculates the probability of n clients in a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    param n: Number of clients
    return: Pn = Probability of n clients"""
    P_0 = mms_model_compute_Pzero(lam, miu, s)
    if n == 0:
        Pn = P_0
    elif n > 0:
        if 0 <= n and n <= s:
            Pn = P_0 * (pow((lam / miu), n) / m.factorial(n))
        elif n > s:
            Pn = P_0 * (pow((lam / miu), n) / (m.factorial(s) * pow(s, n - s)))
    print(f"Probability of {n} clients= ", Pn)
    return Pn
