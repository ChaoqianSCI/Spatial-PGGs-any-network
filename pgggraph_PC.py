
import numpy as np  # 生成备用矩阵nei
from numba import njit, prange
import random
import math
import scipy.io as sio
import datetime
import scipy.stats
from scipy.io import loadmat




TT = 400000
multitime = 1000000000

w = 0.01

ne1 = 1; ne2 = 5


@njit(parallel=True)
def parallel_loop_pgg():

    data = np.zeros((ne1 + 1, ne2 + 1))
    for needle1 in range(1, ne1 + 1):

        for needle2 in range(1, ne2 + 1):
            r = limi + (needle2 - 3) * 1

            ta = np.zeros((multitime + 1))
            for time in prange(1, multitime + 1):

                S = np.zeros((N + 1))
                P = np.zeros((N + 1))

                i = random.randint(1, N)
                S[i] = 1

                for t in range(2, TT + 1):
                    for MCS in range(1, N + 1):
                        i1 = random.randint(1, N)
                        i2 = nei[i1, random.randint(2, numnei[i1] + 1)]

                        # pairwise comparison
                        if S[i1] != S[i2]:

                            for k in [i1, i2]:
                                P[k] = 0
                                for j in range(1, numnei[k] + 1 + 1):
                                    jj = nei[k, j]
                                    nC = np.sum(S[nei[jj, 1:numnei[jj] + 1 + 1]])
                                    P[k] = P[k] + r * nC / (numnei[jj] + 1) - S[k]
                                P[k] = P[k] / (numnei[k] + 1)

                            if random.random() < 1 / (1 + np.exp(w * (P[i1] - P[i2]))):
                                S[i1] = S[i2]

                    aveC = 0
                    for i in range(1, N + 1):
                        aveC = aveC + S[i]
                    if aveC == 0 or aveC == N or t == TT:
                        ta[time] = aveC / N
                        break

            data[needle1, needle2] = np.sum(ta[1:multitime + 1]) / multitime
            print(r, data[needle1, needle2])

    return data



if __name__ == "__main__":

    # star graph, or input other network structure by yourself
    n = 9

    N = n + 1
    nei = np.zeros((N + 1, n + 1 + 1))
    numnei = np.zeros((N + 1))
    for j in range(1, n + 1 + 1):
        nei[1, j] = j
    numnei[1] = n
    for i in range(2, N + 1):
        nei[i, 1] = i
        nei[i, 2] = 1
        numnei[i] = 1
    nei = nei.astype(np.int64)
    numnei = numnei.astype(np.int64)

    limi = -(4 * (3 * n - 1) * (n + 1)) / (- 3 * np.power(n, 2) + 2 * n + 1)


    data = np.zeros((ne1 + 1, ne2 + 1))
    data = parallel_loop_pgg()


    save_fn = 'data.mat'
    sio.savemat(save_fn, {'data': data})