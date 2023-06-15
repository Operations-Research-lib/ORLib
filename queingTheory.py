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
    return Pn


# -------------------------M/M/s-------------------------------------------------
def mms_model_compute_Pzero(lam, miu, s):
    """Calculates the probability of zero clients in a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    return: Pzero = Probability of zero clients"""
    Pzero = 0
    rho = lam / (s * miu)
    for i in range(s):
        Pzero += pow(rho, i) / m.factorial(i)
        Pzero += (pow(rho, s) / (m.factorial(s))) * (1 / (1 - rho))
    Pzero = 1 / Pzero
    return Pzero


def mms_model_info(lam, miu, s):
    """Calculates the basic information of a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    return: system_info = [rho, Pzero, Lq, L, Wq, W]"""
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
    return system_info


def mms_model_compute_Pn(lam, miu, s, n):
    """Calculates the probability of n clients in a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    param n: Number of clients
    return: Pn = Probability of n clients"""
    Pzero = mms_model_compute_Pzero(lam, miu, s)
    Pn = 0
    if n == 0:
        Pn = Pzero
    elif n > 0:
        if 0 <= n and n <= s:
            Pn = Pzero * (pow((lam / miu), n) / m.factorial(n))
        elif n > s:
            Pn = Pzero * (pow((lam / miu), n) / (m.factorial(s) * pow(s, n - s)))
    return Pn


# -------------------------M/M/1/K----------------------------------------------
def mm1k_model_compute_Pzero(lam, miu, k):
    """Calculates the probability of zero clients in a M/M/1/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    return: Pzero = Probability of zero clients"""
    rho = lam / miu
    Pzero = 0
    if rho == 1:
        Pzero = 1 / (k + 1)
    else:
        Pzero = (1 - rho) / (1 - pow(rho, k + 1))
    return Pzero


def mm1k_model_compute_Pn(lam, miu, k, n):
    # TODO(@AboTresol): should we check if k is greater than n?
    """Calculates the probability of n clients in a M/M/1/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param n: Number of clients
    return: Pn = Probability of n clients"""
    rho = lam / miu
    Pn = 0
    if rho == 1:
        Pn = 1 / (k + 1)
    else:
        Pn = ((1 - rho) / (1 - pow(rho, k + 1))) * (pow(rho, n))
    return Pn


def mm1k_model_compute_L(rho, k):
    """Calculates the expected number of clients in a M/M/1/K queueing system.
    param rho: Utilization factor
    param k: Capacity of the system
    return: L = Expected number of clients"""
    L = 0
    if rho == 1:
        L = k / 2
    else:
        L = (rho / (1 - rho)) - ((((k + 1) * pow(rho, k + 1))) / (1 - pow(rho, k + 1)))
    return L


def mm1k_model_compute_avarage_lambda(lam, miu, k):
    """Calculates the avarage arrival rate in a M/M/1/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    return: avarage_lambda = Avarage arrival rate"""
    rho = lam / miu
    Pk = mm1k_model_compute_Pn(lam, miu, k, k)
    avarage_lambda = lam * (1 - Pk)
    return avarage_lambda


def mm1k_model_info(lam, miu, k):
    """Calculates the basic information of a M/M/1/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param n: Number of clients
    return: system_info = [rho, Pzero, L, Lq, Wq, W]"""
    system_info = np.zeros(6)
    avarage_lambda = mm1k_model_compute_avarage_lambda(lam, miu, k)
    # calculate rho, remember that mm1k models are always stable
    system_info[0] = lam / miu
    # calculate probability of zero clients
    system_info[1] = mm1k_model_compute_Pzero(lam, miu, k)
    # calculate L expected number of clients in the system
    system_info[2] = mm1k_model_compute_L(system_info[0], k)
    # calculate the Lq expected lenght of queue
    system_info[3] = system_info[2] - (1 - system_info[1])
    # calculate the Wq expected waiting time in queue
    system_info[4] = system_info[3] / avarage_lambda
    # calculate the W expected waiting time in system
    system_info[5] = system_info[2] / avarage_lambda
    return system_info


# -------------------------M/M/s/K----------------------------------------------
def mmsk_model_compute_Pzero(lam, miu, k, s):
    """Computes the probability of zero clients in the mmsk model
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: Pzero = Probability of zero clients"""
    Pzero = 0
    rho = lam / miu
    first_sum = 0
    second_sum = 0
    for n in range(s):
        first_sum += pow(rho, n) / m.factorial(n)
        first_sum += pow(rho, s) / m.factorial(s)
    for n in range(s + 1, k):
        second_sum += pow(lam / (s * miu), n - s)
    Pzero = 1 / (first_sum * second_sum)
    return Pzero


def mmsk_model_compute_Pn(lam, miu, k, s, n):
    """Computes the probability of n clients in the mmsk model
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    param n: number of clients
    return: Pn = Probability of n clients"""
    Pn = 0
    Pzero = mmsk_model_compute_Pzero(lam, miu, k, s)
    if n < s:
        Pn = (pow((lam / miu), n) / m.factorial(n)) * Pzero
    elif n >= s and n <= k:
        Pn = pow((lam / miu), n) / (m.factorial(s) * pow(s, n - s))
    else:
        Pn = 0
    return Pn


def mmsk_model_compute_Lq(lam, miu, k, s):
    """Calculates the expected number of clients in queue in a M/M/s/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: Lq = Expected number of clients in queue"""
    Pzero = mmsk_model_compute_Pzero(lam, miu, k, s)
    rho = lam / (s * miu)
    Lq = (pow(rho, s + 1) * Pzero) / (m.factorial(s) * pow(1 - rho, 2))
    Lq = Lq * (1 - pow(rho, k - s) - ((k - s) * pow(rho, k - s) * (1 - rho)))
    return Lq


def mmsk_model_compute_L(lam, miu, k, s):
    """Calculates the expected number of clients in a M/M/s/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: L = Expected number of clients"""
    L = 0
    Lq = mmsk_model_compute_Lq(lam, miu, k, s)
    for n in range(s - 1):
        first_sum += n * mmsk_model_compute_Pn(lam, miu, k, s, n) + Lq
        second_sum = 0
        for n in range(s - 1):
            second_sum += mmsk_model_compute_Pn(lam, miu, k, s, n)
        L = first_sum + s * (1 - second_sum)
    return L


def mmsk_model_compute_W(lam, miu, k, s):
    """Calculates the expected waiting time in system in a M/M/s/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: W = Expected waiting time in system"""
    W = 0
    L = mmsk_model_compute_L(lam, miu, k, s)
    W = L / mm1k_model_compute_avarage_lambda(lam, miu, k)
    return W


def mmsk_model_compute_Wq(lam, miu, k, s):
    """Calculates the expected waiting time in queue in a M/M/s/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: Wq = Expected waiting time in queue"""
    Wq = 0
    Lq = mmsk_model_compute_Lq(lam, miu, k, s)
    Wq = Lq / mm1k_model_compute_avarage_lambda(lam, miu, k)
    return Wq


def mmsk_mode_info(lam, miu, k, s):
    """Calculates the system information of a M/M/s/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: system_info = [rho, Pzero, L, Lq, Wq, W]"""
    system_info = np.zeros(6)
    # calculate rho
    system_info[0] = lam / (s * miu)
    # calculate probability of zero clients
    system_info[1] = mmsk_model_compute_Pzero(lam, miu, k, s)
    # calculate L expected number of clients in the system
    system_info[2] = mmsk_model_compute_L(lam, miu, k, s)
    # calculate Lq expected number of clients in queue
    system_info[3] = mmsk_model_compute_Lq(lam, miu, k, s)
    # calculate the Wq expected waiting time in queue
    system_info[4] = mmsk_model_compute_Wq(lam, miu, k, s)
    # calculate the W expected waiting time in system
    system_info[5] = mmsk_model_compute_W(lam, miu, k, s)
    return system_info
