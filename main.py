import numpy as np

#  Karmarkar-Karp algorithm
def KK(arr):
    arr[::-1].sort() # sort in descending order, in place
    counter = len(arr)
    while(counter > 0):
        absdiff = abs(arr[0]-arr[1])
        arr[0] = abs(arr[0]-arr[1])
        arr[1] = 0
        arr[::-1].sort()   # this is not the most efficient! but i'm lazy
        if absdiff:
            counter -= 1
        else:
            counter -= 2
    return arr[0]

def residue(S,A):
    return abs(np.dot(S, A))

def random_solution(arr_len):
    return np.random.randint(2, size=arr_len) * 2 - 1



def repeated_random(A, max_iter, repn):
    if repn == 'sign':
        # generate a random solution
        arr_len = len(A)
        S = random_solution(arr_len)

        # repeat and see if get better solution
        for iter in range(max_iter):
            Sp = random_solution(arr_len)
            if residue(Sp, A) < residue(S, A):
                S = Sp
            print(S)

        return S
    elif repn == 'prepartition':
        pass

def random_neighbor(S, repn):
    if repn == 'sign':
        pass
    elif repn == 'prepartition':
        pass


def hill_climbing(A, max_iter, repn):
    if repn == 'sign':
        # generate a random solution
        arr_len = len(A)
        S = random_solution(arr_len)
        for iter in range(max_iter):
            Sp = random_neighbor(S, repn)
            if residue(Sp, A) < residue(S, A):
                S = Sp
        return S
    elif repn == 'prepartition':
        pass

def annealing(A, max_iter, T, repn):
    if repn == 'sign':
        pass
    elif repn == 'prepartition':
        pass


test = np.array([10,8,7,6,5])
repeated_random(test,100)