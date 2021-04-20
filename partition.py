import os

import numpy as np
import time
from math import floor, e
from kk import KK
import sys



def residue(S, A, arr_len, repn):
    if repn == 'sign':
        return int(abs(np.dot(S, A)))
    elif repn == 'prepartition':
        # first convert back to standard form
        Ap = np.zeros(arr_len)
        for j in range(arr_len):
            Ap[S[j]] += A[j]
        # do KK on the standard form
        return int(KK(Ap))


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


def random_neighbor(S, arr_len, repn):
    Snew = np.copy(S)
    if repn == 'sign':
        i,j = 0,0
        while i == j:
            i, j = np.random.randint(0, arr_len), np.random.randint(0, arr_len)
        Snew[i] = -Snew[i]
        if np.random.random() < 0.5:
            Snew[j] = -Snew[j]
    elif repn == 'prepartition':
        i, j = np.random.randint(0, arr_len), np.random.randint(0, arr_len)
        while Snew[i] == j:  # because we don't want p_i = j, given by problem spec
            i, j = np.random.randint(0, arr_len), np.random.randint(0, arr_len)  # regenerate
        Snew[i] = j
    return Snew


def hill_climbing(A, max_iter, repn, S_init=None):
    # generate a random solution
    arr_len = len(A)
    if S_init is not None:
        S = S_init
    else:
        S = random_solution(arr_len, repn)

    res_S = residue(S, A, arr_len, repn)

    for iter in range(max_iter):
        Sp = random_neighbor(S, arr_len, repn)
        res_Sp = residue(Sp, A, arr_len, repn)
        if res_Sp < res_S:
            S = np.copy(Sp)
            res_S = res_Sp
    return S


def annealing(A, max_iter, T, repn, S_init=None):
    arr_len = len(A)

    # This is for debugging, when we feed the same initial solutions to
    # the hill climbing and annealing algorithms
    if S_init is not None:
        S = S_init
    else: # generate a random solution
        S = random_solution(arr_len, repn)

    Sbest = np.copy(S)
    res_Sbest = residue(Sbest, A, arr_len, repn)
    for iter in range(max_iter):
        Sp = random_neighbor(S, arr_len, repn)
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
    seq = random_sequence(n_ints, upper_bound)
    init_solution = random_solution(n_ints, repn)
    kk_residue = KK(seq)
    repeated_random_residue = residue(repeated_random(seq, max_iters, repn), seq, n_ints, repn)
    hill_climbing_residue = residue(hill_climbing(seq, max_iters, repn), seq, n_ints, repn)
    annealing_residue = residue(annealing(seq, max_iters, T, repn), seq, n_ints, repn)
    return kk_residue, repeated_random_residue, hill_climbing_residue, annealing_residue


def estimated_average_residues(repn, max_iters, n_instances, n_ints, upper_bound):
    total_kk_residue = total_repeated_random_residue = total_hill_climbing_residue = total_annealing_residue = 0
    for repetition in range(n_instances):
        r1, r2, r3, r4 = all_residues(repn, max_iters, n_ints, upper_bound)
        total_kk_residue += r1
        total_repeated_random_residue += r2
        total_hill_climbing_residue += r3
        total_annealing_residue += r4
    return total_kk_residue // n_instances, total_repeated_random_residue // n_instances, \
           total_hill_climbing_residue // n_instances, total_annealing_residue // n_instances

def read_from_file(filename):
    out = []
    with open(filename) as file:
        for line in file.readlines():
            out.append(int(line.rstrip()))
    file.close()
    return out

def generate_and_write(filename):
    pass

# max_iters = 25000
# np.random.seed(123)
#
# start = time.time()
# print(estimated_average_residues('sign', 25000, 100, 100, 10 ** 12))
# end = time.time()
# print("Time elapsed:", end - start)

if __name__ == "__main__":
    flag, code, input_path = sys.argv[1:]
    # print(flag, code, input_path)
    # print(os.getcwd())
    # exit()
    seq = read_from_file(input_path)
    max_iters = 25000
    n_ints = 100
    code = int(code)
    if code == 0:
        print(KK(seq))
    elif code == 1:
        print(residue(repeated_random(seq, max_iters, "sign"), seq, n_ints, "sign"))
    elif code == 2:
        print(residue(hill_climbing(seq, max_iters, "sign"), seq, n_ints, "sign"))
    elif code == 3:
        print(residue(annealing(seq, max_iters, T, "sign"), seq, n_ints, "sign"))
    elif code == 11:
        print(residue(repeated_random(seq, max_iters, "prepartition"), seq, n_ints, "prepartition"))
    elif code == 12:
        print(residue(hill_climbing(seq, max_iters, "prepartition"), seq, n_ints, "prepartition"))
    elif code == 13:
        print(residue(annealing(seq, max_iters, T, "prepartition"), seq, n_ints, "prepartition"))
    else:
        pass

