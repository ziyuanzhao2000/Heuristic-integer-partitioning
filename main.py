import numpy as np
import time
from math import floor, e


#  Karmarkar-Karp algorithm
def KK(arr):
    arr[::-1].sort()  # sort in descending order, in place
    counter = len(arr)
    while (counter > 0):
        absdiff = abs(arr[0] - arr[1])
        arr[0] = abs(arr[0] - arr[1])
        arr[1] = 0
        arr[::-1].sort()  # this is not the most efficient! but i'm lazy
        if absdiff:
            counter -= 1
        else:
            counter -= 2
    return arr[0]


def residue(S, A, arr_len, repn):
    if repn == 'sign':
        return abs(np.dot(S, A))
    elif repn == 'prepartition':
        # first convert back to standard form
        Ap = np.zeros(arr_len)
        for j in range(arr_len):
            Ap[S[j]] += A[j]
        # do KK on the standard form
        return KK(Ap)


def random_solution(arr_len, repn):
    if repn == 'sign':
        return np.random.randint(2, size=arr_len) * 2 - 1
    elif repn == 'prepartition':
        return np.random.randint(0, arr_len, size=arr_len)


def repeated_random(A, max_iter, repn):
    # generate a random solution
    arr_len = len(A)
    S = random_solution(arr_len, repn)
    res_S = residue(S, A, arr_len, repn)

    # repeat and see if get better solution
    for iter in range(max_iter):
        Sp = random_solution(arr_len, repn)
        if residue(Sp, A, arr_len, repn) < res_S:
            S = Sp
            res_S = residue(S, A, arr_len, repn)

    return S


def random_neighbor(S, index_list, arr_len, repn):
    Snew = np.copy(S)
    if repn == 'sign':
        indices = np.random.choice(index_list, 2, replace=False)
        Snew[indices[0]] *= -1
        if np.random.random() > 0.5:
            Snew[indices[1]] *= -1
    elif repn == 'prepartition':
        i, j = np.random.randint(0, arr_len), np.random.randint(0, arr_len)
        while S[i] == j:  # because we don't want p_i = j, given by problem spec
            i, j = np.random.randint(0, arr_len), np.random.randint(0, arr_len)  # regenerate
        S[i] = j
    return Snew


def hill_climbing(A, max_iter, repn, S_init=None):
    # generate a random solution
    arr_len = len(A)
    index_list = list(range(arr_len))
    if S_init is not None:
        S = S_init
    else:
        S = random_solution(arr_len, repn)

    res_S = residue(S, A, arr_len, repn)

    for iter in range(max_iter):
        Sp = random_neighbor(S, index_list, arr_len, repn)
        if residue(Sp, A, arr_len, repn) < res_S:
            S = np.copy(Sp)
            res_S = residue(S, A, arr_len, repn)
    return S


def annealing(A, max_iter, T, repn, S_init=None):
    # generate a random solution
    arr_len = len(A)
    index_list = list(range(arr_len))

    # This is for debugging, when we feed the same initial solutions to
    # the hill climbing and annealing algorithms
    if S_init is not None:
        S = S_init
    else:
        S = random_solution(arr_len, repn)

    Sbest = np.copy(S)
    res_Sbest = residue(Sbest, A, arr_len, repn)
    for iter in range(max_iter):
        Sp = random_neighbor(S, index_list, arr_len, repn)
        res_Sp = residue(Sp, A, arr_len, repn)
        res_S = residue(S, A, arr_len, repn)
        if res_Sp < res_S:
            S = np.copy(Sp)
        elif np.random.random() < e ** (-(res_Sp - res_S) / T(iter)):
            S = np.copy(Sp)
        if res_S < res_Sbest:
            Sbest = np.copy(S)
            res_Sbest = residue(Sbest, A, arr_len, repn)

    return Sbest


def T(iter):
    return 10 ** 10 * 0.8 ** floor(iter / 300)


def random_sequence(length, upper_bound):
    return np.random.randint(1, upper_bound + 1, size=length)


def all_residues(repn, max_iters, n_ints, upper_bound):
    test = random_sequence(n_ints, upper_bound)
    init_solution = random_solution(n_ints, repn)
    random_residue = residue(init_solution, test, n_ints, repn)
    repeated_random_residue = residue(repeated_random(test, max_iters, repn), test, n_ints, repn)
    hill_climbing_residue = residue(hill_climbing(test, max_iters, repn), test, n_ints, repn)
    annealing_residue = residue(annealing(test, max_iters, T, repn), test, n_ints, repn)
    return random_residue, repeated_random_residue, hill_climbing_residue, annealing_residue


def estimated_average_residues(repn, max_iters, n_instances, n_ints, upper_bound):
    total_random_residue = total_repeated_random_residue = total_hill_climbing_residue = total_annealing_residue = 0
    for repetition in range(n_instances):
        r1, r2, r3, r4 = all_residues(repn, max_iters, n_ints, upper_bound)
        total_random_residue += r1
        total_repeated_random_residue += r2
        total_hill_climbing_residue += r3
        total_annealing_residue += r4
    return total_random_residue // n_instances, total_repeated_random_residue // n_instances, \
           total_hill_climbing_residue // n_instances, total_annealing_residue // n_instances


max_iters = 25000
np.random.seed(125)

start = time.time()
print(estimated_average_residues('prepartition', 25000, 50, 100, 10 ** 12))
end = time.time()
print("Time elapsed:", end - start)
