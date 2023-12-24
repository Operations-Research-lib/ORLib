import math as m
import numpy as np


# -------------------------M/M/1-------------------------------------------------
# all the functions in this section are based on the equations of the M/M/1 model
# as described in the book "Introduction to Operations Research" by Hillier and Lieberman
# note: all functions have been tested and work properly
UNSTABLE_MESSAGE = "The system is unstable"
NEGATIVE_INPUT_ERROR = "Lamda and miu cannot be negative or 0"


def mm1_model_info(lam, miu):
    """computes the basic information of a M/M/1 queueing system.
    param lam: Arrival rate
    param miu: Service rate
    return: system_info = [rho, Lq, L, Wq, W]"""
    if lam == miu:
        print(UNSTABLE_MESSAGE)
        return None
    system_info = np.zeros(5)
    system_info[0] = mm1_model_compute_rho(lam, miu)
    system_info[1] = mm1_model_compute_Lq(lam, miu)
    system_info[2] = mm1_model_compute_L(lam, miu)
    system_info[3] = mm1_model_compute_Wq(lam, miu)
    system_info[4] = mm1_model_compute_W(lam, miu)
    return system_info


def mm1_model_compute_rho(lam, miu):
    """computes the utilization factor of a M/M/1 queueing system.
    param lam: Arrival rate
    param miu: Service rate
    return: rho = Utilization factor"""
    rho = lam / miu
    return rho


def mm1_model_compute_Lq(lam, miu):
    """computes the average number of clients in the queue of a M/M/1 queueing system.
    param lam: Arrival rate
    param miu: Service rate
    return: Lq = Average number of clients in the queue"""
    if lam == miu:
        print("The system is unstable")
        return None
    Lq = pow(lam, 2) / (miu * (miu - lam))
    return Lq


def mm1_model_compute_L(lam, miu):
    """computes the average number of clients in the system of a M/M/1 queueing system.
    param lam: Arrival rate
    param miu: Service rate
    return: L = Average number of clients in the system"""
    if lam == miu:
        print("The system is unstable")
        return None
    L = lam / (miu - lam)
    return L


def mm1_model_compute_Wq(lam, miu):
    """computes the average waiting time in the queue of a M/M/1 queueing system.
    param lam: Arrival rate
    param miu: Service rate
    return: Wq = Average waiting time in the queue"""
    if lam == miu:
        print("The system is unstable")
        return None
    Wq = lam / (miu * (miu - lam))
    return Wq


def mm1_model_compute_W(lam, miu):
    """computes the average waiting time in the system of a M/M/1 queueing system.
    param lam: Arrival rate
    param miu: Service rate
    return: W = Average waiting time in the system"""
    if lam == miu:
        print("The system is unstable")
        return None
    W = 1 / (miu - lam)
    return W


def mm1_model_compute_Pn(lam, miu, n):
    """computes the probability of n clients in a M/M/1 queueing system.
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


# -------------------------M/M/s------------------------------------------------
# all the functions in this section are based on the equations of the M/M/s model
# as described in the book "Introduction to Operations Research" by Hillier and Lieberman
# note: all functions have been tested and work properly


def mms_model_compute_p_zero(lam, miu, s):
    """computes the probability of zero clients in a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    return: Pzero = Probability of zero clients"""
    if lam <= 0 or miu <= 0 or s <= 0:
        raise NameError(NEGATIVE_INPUT_ERROR)
    rho = lam / (s * miu)
    if rho >= 1:
        raise NameError(UNSTABLE_MESSAGE)
    denominator = 0
    for n in range(s):
        denominator += pow((lam / miu), n) / m.factorial(n)
    denominator += (pow((lam / miu), s) / (m.factorial(s))) * (1 / (1 - rho))
    p_zero = 1 / denominator
    return p_zero


def mms_model_info(lam, miu, s):
    """computes the basic information of a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    return: system_info = [rho, Pzero, Lq, L, Wq, W]"""
    # calculate rho to determine if the system is stable
    if lam <= 0 or miu <= 0 or s <= 0:
        raise NameError(NEGATIVE_INPUT_ERROR)
    rho = lam / (s * miu)
    if rho >= 1:
        raise NameError(UNSTABLE_MESSAGE)
    # calculate the basic information
    system_info = np.zeros(6)
    system_info[0] = rho
    # calculate probability of zero clients
    system_info[1] = mms_model_compute_p_zero(lam, miu, s)
    # calculate expected lenght of queue
    system_info[2] = mms_model_compute_lq(lam, miu, s)
    # calculate the expected waiting time in queue
    system_info[3] = mms_model_compute_wq(lam, miu, s)
    # calculate the expected waiting time in system
    system_info[4] = mms_model_compute_w(lam, miu, s)
    # calculate the expected number of clients in the system
    system_info[5] = mms_model_compute_l(lam, miu, s)
    return system_info


def mms_model_compute_lq(lam, miu, s):
    """computes the average number of clients in the queue of a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    return: Lq = Average number of clients in the queue"""
    if lam <= 0 or miu <= 0 or s <= 0:
        raise NameError(NEGATIVE_INPUT_ERROR)
    rho = lam / (s * miu)
    if rho >= 1:
        raise NameError(UNSTABLE_MESSAGE)
    p_zero = mms_model_compute_p_zero(lam, miu, s)
    lq = (p_zero * pow((lam / miu), s) * rho) / (m.factorial(s) * pow(1 - rho, 2))
    return lq


def mms_model_compute_wq(lam, miu, s):
    """computes the average waiting time in the queue of a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    return: Wq = Average waiting time in the queue"""
    lq = mms_model_compute_lq(lam, miu, s)
    wq = lq / lam
    return wq


def mms_model_compute_w(lam, miu, s):
    """computes the average waiting time in the system of a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    return: W = Average waiting time in the system"""
    wq = mms_model_compute_wq(lam, miu, s)
    w = wq + (1 / miu)
    return w


def mms_model_compute_l(lam, miu, s):
    """computes the average number of clients in the system of a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    return: L = Average number of clients in the system"""
    lq = mms_model_compute_lq(lam, miu, s)
    l = lq + (lam / miu)
    return l


def mms_model_compute_p_n(lam, miu, s, n):
    """computes the probability of n clients in a M/M/s queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param s: Number of servers
    param n: Number of clients
    return: Pn = Probability of n clients"""
    p_zero = mms_model_compute_p_zero(lam, miu, s)
    p_n = 0
    if n == 0:
        p_n = p_zero
    elif n > 0:
        if 0 <= n and n <= s:
            p_n = p_zero * (pow((lam / miu), n) / m.factorial(n))
        elif n > s:
            p_n = p_zero * (pow((lam / miu), n) / (m.factorial(s) * pow(s, n - s)))
    return p_n


# -------------------------M/M/1/K----------------------------------------------
# all the functions in this section are based on the equations of the M/M/s model
# as described in the book "Introduction to Operations Research" by Hillier and Lieberman
# note: all functions have been tested and work properly
def mm1k_model_compute_Pzero(lam, miu, k):
    """computes the probability of zero clients in a M/M/1/K queueing system.
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
    """computes the probability of n clients in a M/M/1/K queueing system.
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


def mm1k_model_compute_L(lam, miu, k):
    """computes the expected number of clients in a M/M/1/K queueing system.
    param rho: Utilization factor
    param k: Capacity of the system
    return: L = Expected number of clients"""
    L = 0
    rho = lam / miu
    if rho == 1:
        L = k / 2
    else:
        L = (rho / (1 - rho)) - ((k + 1) * pow(rho, k + 1) / (1 - pow(rho, k + 1)))
    return L


def mm1k_model_compute_Lq(lam, miu, k):
    """computes the expected number of clients in the queue of a M/M/1/K queueing system.
    param rho: Utilization factor
    param k: Capacity of the system
    return: Lq = Expected number of clients in the queue"""
    Lq = 0
    L = mm1k_model_compute_L(lam, miu, k)
    Pzero = mm1k_model_compute_Pzero(lam, miu, k)
    Lq = L - (1 - Pzero)
    return Lq


def mm1k_model_compute_W(lam, miu, k):
    """computes the expected waiting time in a M/M/1/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    return: W = Expected waiting time"""
    W = 0
    L = mm1k_model_compute_L(lam, miu, k)
    W = L / mm1k_model_compute_average_lambda(lam, miu, k)
    return W


def mm1k_model_compute_Wq(lam, miu, k):
    """computes the expected waiting time in the queue of a M/M/1/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    return: Wq = Expected waiting time in the queue"""
    Wq = 0
    Lq = mm1k_model_compute_Lq(lam, miu, k)
    Wq = Lq / mm1k_model_compute_average_lambda(lam, miu, k)
    return Wq


def mm1k_model_compute_average_lambda(lam, miu, k):
    """computes the average arrival rate in a M/M/1/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    return: average_lambda = average arrival rate"""
    Pk = mm1k_model_compute_Pn(lam, miu, k, k)
    average_lambda = lam * (1 - Pk)
    return average_lambda


def mm1k_model_info(lam, miu, k):
    """computes the basic information of a M/M/1/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param n: Number of clients
    return: system_info = [rho, Pzero, L, Lq, Wq, W]"""
    system_info = np.zeros(6)
    # calculate rho, remember that mm1k models are always stable
    system_info[0] = lam / miu
    # calculate probability of zero clients
    system_info[1] = mm1k_model_compute_Pzero(lam, miu, k)
    # calculate L expected number of clients in the system
    system_info[2] = mm1k_model_compute_L(lam, miu, k)
    # calculate the Lq expected lenght of queue
    system_info[3] = mm1k_model_compute_Lq(lam, miu, k)
    # calculate the Wq expected waiting time in queue
    system_info[4] = mm1k_model_compute_Wq(lam, miu, k)
    # calculate the W expected waiting time in system
    system_info[5] = mm1k_model_compute_W(lam, miu, k)
    return system_info


# -------------------------M/M/s/K----------------------------------------------
# all the functions in this section are based on the equations of the M/M/s model
# as described in the book "Introduction to Operations Research" by Hillier and Lieberman
# note: all functions have been tested and work properly
# note: do not use these functions with s = 1, use the mm1k functions instead
def mmsk_model_compute_p_zero(lam, miu, k, s):
    """Computes the probability of zero clients in the mmsk model
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: Pzero = Probability of zero clients"""
    rho = lam / (s * miu)
    first_sum = sum(pow((lam / miu), n) / m.factorial(n) for n in range(s + 1))
    constant = pow((lam / miu), s) / m.factorial(s)
    second_sum = sum(pow(rho, n - s) for n in range(s + 1, k + 1))
    p_zero = 1 / (first_sum + (constant * second_sum))
    return p_zero


def mmsk_model_compute_average_lambda(lam, miu, k, s):
    """computes the average arrival rate in a M/M/s/k queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: average_lambda = average arrival rate"""
    p_k = mmsk_model_compute_p_n(lam, miu, k, s, k)
    average_lambda = lam * (1 - p_k)
    return average_lambda


def mmsk_model_compute_p_n(lam, miu, k, s, n):
    """Computes the probability of n clients in the mmsk model
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    param n: number of clients
    return: Pn = Probability of n clients"""
    p_n = 0
    p_zero = mmsk_model_compute_p_zero(lam, miu, k, s)
    if n < s:
        p_n = (pow((lam / miu), n) / m.factorial(n)) * p_zero
    elif s <= n <= k:
        p_n = (pow((lam / miu), n) / (m.factorial(s) * pow(s, n - s))) * p_zero
    elif n > k:
        p_n = 0
    return p_n


def mmsk_model_compute_lq(lam, miu, k, s):
    """computes the expected number of clients in queue in a M/M/s/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: Lq = Expected number of clients in queue"""
    p_zero = mmsk_model_compute_p_zero(lam, miu, k, s)
    rho = lam / (s * miu)
    lq = (pow((lam / miu), s) * p_zero * rho) / (m.factorial(s) * pow(1 - rho, 2))
    lq *= 1 - pow(rho, k - s) - ((k - s) * pow(rho, k - s) * (1 - rho))
    return lq


def mmsk_model_compute_l(lam, miu, k, s):
    """computes the expected number of clients in a M/M/s/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: L = Expected number of clients"""
    l = 0
    lq = mmsk_model_compute_lq(lam, miu, k, s)
    first_sum = 0
    second_sum = 0
    for n in range(s):
        p_n = mmsk_model_compute_p_n(lam, miu, k, s, n)
        first_sum += n * p_n
        second_sum += p_n
    l = first_sum + lq + s * (1 - second_sum)
    return l


def mmsk_model_compute_w(lam, miu, k, s):
    """computes the expected waiting time in system in a M/M/s/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: W = Expected waiting time in system"""
    l = mmsk_model_compute_l(lam, miu, k, s)
    w = l / mmsk_model_compute_average_lambda(lam, miu, k, s)
    return w


def mmsk_model_compute_wq(lam, miu, k, s):
    """computes the expected waiting time in queue in a M/M/s/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: Wq = Expected waiting time in queue"""
    wq = 0
    lq = mmsk_model_compute_lq(lam, miu, k, s)
    wq = lq / mmsk_model_compute_average_lambda(lam, miu, k, s)
    return wq


def mmsk_model_info(lam, miu, k, s):
    """computes the system information of a M/M/s/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: system_info = [rho, Pzero, L, Lq, Wq, W]"""
    system_info = np.zeros(6)
    # calculate rho
    system_info[0] = lam / (s * miu)
    # calculate probability of zero clients
    system_info[1] = mmsk_model_compute_p_zero(lam, miu, k, s)
    # calculate L expected number of clients in the system
    system_info[2] = mmsk_model_compute_l(lam, miu, k, s)
    # calculate Lq expected number of clients in queue
    system_info[3] = mmsk_model_compute_lq(lam, miu, k, s)
    # calculate the Wq expected waiting time in queue
    system_info[4] = mmsk_model_compute_wq(lam, miu, k, s)
    # calculate the W expected waiting time in system
    system_info[5] = mmsk_model_compute_w(lam, miu, k, s)
    return system_info


# -------------General Birth-Death Model----------------------------------------
# all functions are based on equations from Introduction to Operations Research
# by Hillier and Lieberman ed.11 chapter 17
# Note: all functions have been tested and work properly
def birth_death_model_compute_Cn(lamdas, mius, n):
    """Computes the probability of n clients in the birth-death model
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    param n: number of clients
    return: Cn = Probability of n clients"""
    Cn = 0
    if n == 0:
        Cn = 1
    else:
        Cn = (lamdas[n - 1] / mius[n - 1]) * birth_death_model_compute_Cn(
            lamdas, mius, n - 1
        )
    return Cn


def birth_death_model_compute_Pzero(lamdas, mius):
    """Computes the probability of zero clients in the birth-death model
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: Pzero = Probability of zero clients"""
    Pzero = 0
    sum_Cn = 0
    for n in range(len(lamdas) + 1):
        Cn = birth_death_model_compute_Cn(lamdas, mius, n)
        sum_Cn += Cn
    Pzero = 1 / (sum_Cn)
    return Pzero


def birth_death_model_compute_Pn(lamdas, mius, n):
    """Computes the probability of n clients in the birth-death model
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    param n: number of clients
    return: Pn = Probability of n clients"""
    Pn = 0
    if n >= 0 and n <= len(lamdas):
        Pn = birth_death_model_compute_Cn(
            lamdas, mius, n
        ) * birth_death_model_compute_Pzero(lamdas, mius)
    else:
        Pn = 0
    return Pn


def birth_death_model_compute_L(lamdas, mius):
    """Calculates the expected number of clients in a birth-death queueing system.
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: L = Expected number of clients"""
    L = 0
    for n in range(len(lamdas) + 1):
        L += n * birth_death_model_compute_Pn(lamdas, mius, n)
    return L


def birth_death_model_compute_Lq(lamdas, mius, s):
    """Calculates the expected number of clients in queue in a birth-death queueing system.
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: Lq = Expected number of clients in queue"""
    Lq = 0
    for n in range(s + 1, len(lamdas) + 1):
        Lq += (n - s) * birth_death_model_compute_Pn(lamdas, mius, n)
    return Lq


def birth_death_model_compute_average_lambda(lamdas, mius):
    """Calculates the average arrival rate in a birth-death queueing system.
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: average_lambda = average arrival rate"""
    average_lambda = 0
    for n in range(len(lamdas)):
        average_lambda += lamdas[n] * birth_death_model_compute_Pn(lamdas, mius, n)
    return average_lambda


def birth_death_model_compute_W(lamdas, mius):
    """Calculates the expected waiting time in system in a birth-death queueing system.
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: W = Expected waiting time in system"""
    W = 0
    L = birth_death_model_compute_L(lamdas, mius)
    average_lambda = birth_death_model_compute_average_lambda(lamdas, mius)
    W = L / average_lambda
    return W


def birth_death_model_compute_Wq(lamdas, mius, s):
    """Calculates the expected waiting time in queue in a birth-death queueing system.
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: Wq = Expected waiting time in queue"""
    Wq = 0
    Lq = birth_death_model_compute_Lq(lamdas, mius, s)
    average_lambda = birth_death_model_compute_average_lambda(lamdas, mius)
    Wq = Lq / average_lambda
    return Wq


def birth_death_model_info(lamdas, mius, s):
    """Computes the system information of a birth-death queueing system
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    param s: number of servers
    param upper_bound: upper bound of the system
    return: system_info = [Pzero, L, Lq, W, Wq, average_lambda]"""
    system_info = np.zeros(6)
    # Pzero
    system_info[0] = birth_death_model_compute_Pzero(lamdas, mius)
    # L
    system_info[1] = birth_death_model_compute_L(lamdas, mius)
    # Lq
    system_info[2] = birth_death_model_compute_Lq(lamdas, mius, s)
    # W
    system_info[3] = birth_death_model_compute_W(lamdas, mius)
    # Wq
    system_info[4] = birth_death_model_compute_Wq(lamdas, mius, s)
    # average_lambda
    system_info[5] = birth_death_model_compute_average_lambda(lamdas, mius)


# -------------General Queueing theory formulas---------------------------------
# all functions are based on equations from Introduction to Operations Research
# by Hillier and Lieberman ed.11 chapter 17
# Note: all functions have been tested and work properly
def queueing_theory_compute_L(probabilities):
    """ "computes the expected number of clients in the system
    based on a vector of probabilities
    param: probabilities = vector of probabilities
    return: L = expected number of clients"""
    L = 0
    for n in range(len(probabilities)):
        L += n * probabilities[n]
    return L


def queueing_theory_compute_Lq(probabilities, s):
    """ "computes the expected number of clients in the queue
    based on a vector of probabilities
    param: probabilities = vector of probabilities
    param: s = number of servers
    return: Lq = expected number of clients in the queue"""
    Lq = 0
    for n in range(s, len(probabilities)):
        Lq += (n - s) * probabilities[n]
    return Lq


def queueing_theory_compute_W_from_Wq(Wq, miu):
    """ "computes the expected waiting time in the system
    based on a vector of probabilities
    param: Wq = expected waiting time in the queue
    param: miu = average service rate
    return: W = expected waiting time in the system"""
    W = Wq + 1 / miu
    return W


def queueing_theory_compute_W_from_L(L, lam):
    """ "computes the expected waiting time in the system
    based on a vector of probabilities
    param: L = expected number of clients in the system
    param: lam = average arrival rate
    return: W = expected waiting time in the system"""
    W = L / lam
    return W


def queueing_theory_compute_Wq(Lq, lam):
    """ "computes the expected waiting time in the queue
    based on a vector of probabilities
    param: Lq = expected number of clients in the queue
    param: lam = average arrival rate
    return: Wq = expected waiting time in the queue"""
    Wq = Lq / lam
    return Wq


def queueing_theory_compute_Ls_from_probs(probabilities, s):
    """computes expected number of customers being served in the system
    based on a vector of probabilities
    param: probabilities = vector of probabilities
    param: s = number of servers
    return: Ls = expected number of customers being served in the system"""
    L = queueing_theory_compute_L(probabilities)
    Lq = queueing_theory_compute_Lq(probabilities, s)
    Ls = L - Lq
    return Ls


def queueing_theory_compute_Ls(L, Lq):
    """computes expected number of customers being served in the system
    param: L = expected number of clients in the system
    param: Lq = expected number of clients in the queue
    return: Ls = expected number of customers being served in the system"""
    Ls = L - Lq
    return Ls


def queueing_theory_compute_service_time_from_waits(W, Wq):
    """ "computes the expected average waiting time to complete a service
    param: W = expected waiting time in the system
    param: Wq = expected waiting time in the queue
    return: expected_service_time"""
    average_service_time = W - Wq
    return average_service_time


def queueing_theory_compute_service_rate_from_waits(W, Wq):
    """ "computes the expected average waiting time to complete a service
    param: W = expected waiting time in the system
    param: Wq = expected waiting time in the queue
    return: expected_service_time"""
    average_service_rate = 1 / (W - Wq)
    return average_service_rate


def queueing_theory_compute_service_time_from_miu(miu):
    """ "computes the expected average waiting time to complete a service
    param: miu = average service rate
    return: expected_service_time"""
    expected_service_time = 1 / miu
    return expected_service_time


def queueing_theory_compute_interarrival_time(lam):
    """computes the expected average waiting time to complete a service
    param: lam = average arrival rate
    return: Ws = expected waiting time in the system"""
    expected_interarrival_time = 1 / lam
    return expected_interarrival_time


def queueing_theory_compute_rho(lam, miu, s=1):
    """Computes the utilization factor of the system
    param: lam = average arrival rate
    param: miu = average service rate
    param: s = number of servers
    return: rho = utilization factor"""
    rho = lam / (s * miu)
    return rho


def queueing_theory_get_miu_from_waits(W, Wq):
    """computes the average service rate from the waiting times
    param: W = expected waiting time in the system
    param: Wq = expected waiting time in the queue
    return: miu = average service rate"""
    miu = 1 / (W - Wq)
    return miu


# -------------Poisson distribution---------------------------------------------
# note: all functions have been tested and work properly
def poisson_distribution_events_in_interval(x, lam, t):
    """ "computes the probability of x arrivals in time t
    based on the poisson distribution
    param: x = number of events
    param: lam = average number of events per unit time
    param: t = time interval"""
    prob = 0
    prob = (lam * t) ** x * np.exp(-lam * t) / m.factorial(x)
    return prob


def poisson_distribution_mass_function(lam, k):
    """ "computes the probability of k arrivals in time t
    based on the poisson distribution
    param: lam = average number of events per unit time
    param: k = number of events
    return: prob = probability of k arrivals in time t"""
    prob = 0
    prob = (lam) ** k * np.exp(-lam) / m.factorial(k)
    return prob


def poisson_distribution_sum_P(lam, t, up_to):
    """ "computes the probability of up_to arrivals in time t
    based on the poisson distribution
    param: lam = average number of events per unit time
    param: t = time interval
    param: up_to = number of events
    return: prob = probability of up_to arrivals in time t"""
    prob = 0
    for x in range(up_to + 1):
        prob += poisson_distribution_events_in_interval(x, lam, t)
    return prob


def poison_distribution_expected_events(lam, t):
    """ "computes the expected number of events in time t
    based on the poisson distribution
    param: lam = average number of events per unit time
    param: t = time interval
    return: expected_events = expected number of events in time t"""
    expected_events = 0
    expected_events = lam * t
    return expected_events


# -------------Exponential distribution-----------------------------------------
# note: all functions have been tested and work properly
def exponential_distribution_compute_P(t, alpha):
    """ "computes the probability of a service time greater than t
    based on the exponential distribution
    param: t = time interval
    param: alpha = average number of events per unit time
    return: prob = probability of a service time greater than t"""
    prob = 0
    if t >= 0:
        prob = alpha * np.exp(-alpha * t)
    elif t < 0:
        prob = 0
    return prob


def exponential_cumulative_distribution(t, alpha, less):
    """ "computes the commulative probability of an exponential distribution
    param: t = time interval
    param: alpha = average number of events per unit time
    param: less = chooses if less equal or bigger than is desired
    return: prob = probability of a service time greater than t"""
    prob = np.exp(-alpha * t)
    if less:
        prob = 1 - prob
    return prob
