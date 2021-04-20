#  Karmarkar-Karp algorithm
from bisect import insort_left
import sys
def KK(arr):
    # arr[::-1].sort()  # sort in descending order, in place
    # counter = len(arr)
    # while (counter > 0):
    #     absdiff = abs(arr[0] - arr[1])
    #     arr[0] = abs(arr[0] - arr[1])
    #     arr[1] = 0
    #     arr[::-1].sort()  # this is not the most efficient! but i'm lazy
    #     if absdiff:
    #         counter -= 1
    #     else:
    #         counter -= 2
    # return arr[0]

    # first convert to python list then use bisect.insort_left
    l = list(arr)
    l.sort()
    while len(l) > 1:
        insort_left(l, abs(l.pop() - l.pop())) # pop two (largest) from the right, find abs diff, then insert back
    return l[0]

def read_from_file(filename):
    out = []
    with open(filename) as file:
        for line in file.readlines():
            out.append(int(line.rstrip()))
    file.close()
    return out

if __name__ == "__main__":
    print(KK(read_from_file(sys.argv[1]+".txt")))