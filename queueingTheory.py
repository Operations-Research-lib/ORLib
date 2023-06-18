import math as m
import numpy as np


# -------------------------M/M/1-------------------------------------------------
def mm1_model_info(lam, miu):
    """computes the basic information of a M/M/1 queueing system.
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


# -------------------------M/M/s-------------------------------------------------
def mms_model_compute_Pzero(lam, miu, s):
    """computes the probability of zero clients in a M/M/s queueing system.
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
    """computes the basic information of a M/M/s queueing system.
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
    """computes the probability of n clients in a M/M/s queueing system.
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


def mm1k_model_compute_L(rho, k):
    """computes the expected number of clients in a M/M/1/K queueing system.
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
    """computes the avarage arrival rate in a M/M/1/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    return: avarage_lambda = Avarage arrival rate"""
    rho = lam / miu
    Pk = mm1k_model_compute_Pn(lam, miu, k, k)
    avarage_lambda = lam * (1 - Pk)
    return avarage_lambda


def mm1k_model_info(lam, miu, k):
    """computes the basic information of a M/M/1/K queueing system.
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
    """computes the expected number of clients in queue in a M/M/s/K queueing system.
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
    """computes the expected number of clients in a M/M/s/K queueing system.
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
    """computes the expected waiting time in system in a M/M/s/K queueing system.
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
    """computes the expected waiting time in queue in a M/M/s/K queueing system.
    param lam: Arrival rate
    param miu: Service rate
    param k: Capacity of the system
    param s: number of servers
    return: Wq = Expected waiting time in queue"""
    Wq = 0
    Lq = mmsk_model_compute_Lq(lam, miu, k, s)
    Wq = Lq / mm1k_model_compute_avarage_lambda(lam, miu, k)
    return Wq


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


# -------------General Birth-Death Model-----------------------------------------
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


def birth_death_model_compute_Pzero(lamdas, mius, upper_bound):
    """Computes the probability of zero clients in the birth-death model
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: Pzero = Probability of zero clients"""
    Pzero = 0
    first_sum = 0
    second_sum = 0

    Pzero = 0
    sum_Cn = 0
    for n in range(upper_bound):
        Cn = birth_death_model_compute_Cn(lamdas, mius, n)
        sum_Cn += Cn
    Pzero = 1 / (sum_Cn)
    return Pzero


def birth_death_model_compute_Pn(lamdas, mius, upper_bound, n):
    """Computes the probability of n clients in the birth-death model
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    param n: number of clients
    return: Pn = Probability of n clients"""
    Pn = 0
    if n >= 0 and n <= upper_bound:
        Pn = birth_death_model_compute_Cn(
            lamdas, mius, n
        ) * birth_death_model_compute_Pzero(lamdas, mius, upper_bound)
    else:
        Pn = 0
    return Pn


def birth_death_model_compute_L(lamdas, mius, upper_bound):
    """Calculates the expected number of clients in a birth-death queueing system.
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: L = Expected number of clients"""
    L = 0
    for n in range(upper_bound):
        L += n * birth_death_model_compute_Pn(lamdas, mius, upper_bound, n)
    return L


def birth_death_model_compute_Lq(lamdas, mius, s, upper_bound):
    """Calculates the expected number of clients in queue in a birth-death queueing system.
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: Lq = Expected number of clients in queue"""
    Lq = 0
    for n in range(1, upper_bound):
        Lq += (n - s) * birth_death_model_compute_Pn(lamdas, mius, upper_bound, n)
    return Lq


def birth_death_model_compute_average_lambda(lamdas, mius, upper_bound):
    """Calculates the average arrival rate in a birth-death queueing system.
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: average_lambda = average arrival rate"""
    average_lambda = 0
    for n in range(upper_bound):
        average_lambda += lamdas[n] * birth_death_model_compute_Pn(
            lamdas, mius, upper_bound, n
        )
    return average_lambda


def birth_death_model_compute_W(lamdas, mius, upper_bound):
    """Calculates the expected waiting time in system in a birth-death queueing system.
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: W = Expected waiting time in system"""
    W = 0
    L = birth_death_model_compute_L(lamdas, mius, upper_bound)
    average_lambda = birth_death_model_compute_average_lambda(lamdas, mius, upper_bound)
    W = L / average_lambda
    return W


def birth_death_model_compute_Wq(lamdas, mius, s, upper_bound):
    """Calculates the expected waiting time in queue in a birth-death queueing system.
    param lamdas: an array of Arrival rates. Starting from 0 to n
    param mius: an array of Service rates. Starting from 1 to n
    return: Wq = Expected waiting time in queue"""
    Wq = 0
    Lq = birth_death_model_compute_Lq(lamdas, mius, s, upper_bound)
    average_lambda = birth_death_model_compute_average_lambda(lamdas, mius, upper_bound)
    Wq = Lq / average_lambda
    return Wq


# -------------General Queueing theory formulas---------------------------------
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


def queueing_theory_compute_W(Wq, miu):
    """ "computes the expected waiting time in the system
    based on a vector of probabilities
    param: Wq = expected waiting time in the queue
    param: miu = average service rate
    return: W = expected waiting time in the system"""
    W = Wq + 1 / miu
    return W


def queueing_theory_compute_W(L, lam):
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


def queueing_theory_compute_Ls(probabilities, s):
    """ "computes the expected number of clients in the system
    based on a vector of probabilities
    param: probabilities = vector of probabilities
    param: s = number of servers
    return: Ls = expected number of clients in the system"""
    L = queueing_theory_compute_L(probabilities)
    Lq = queueing_theory_compute_Lq(probabilities, s)
    Ls = L - Lq
    return Ls


def queueing_theory_compute_Ls(L, Lq):
    """ "computes the expected number of clients in the system
    based on a vector of probabilities
    param: L = expected number of clients in the system
    param: Lq = expected number of clients in the queue
    return: Ls = expected number of clients in the system"""
    Ls = L - Lq
    return Ls


# -------------Poisson distribution---------------------------------
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
    for x in range(up_to):
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


# -------------Exponential distribution---------------------------------
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
