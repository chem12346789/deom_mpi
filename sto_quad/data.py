import time
import pandas as pd
import numpy as np
from functools import partial
import itertools

fastread = partial(pd.read_csv, header=None, sep='\s+', dtype=np.float64)


def fastread_np(str):
    return fastread(str).to_numpy()


def date_stamp():
    return time.strftime("%b_%d_%H_%M_%S_%Y", time.localtime()) + '.pdf'


dir_1 = './{}/prop-rho.dat_0.01_{}_False_'
# num_line1_l = []
# alp2_l = [0.1, -0.1, 0.05, -0.05]

# for alp2, trun in itertools.product(alp2_l, [1e20]):
#     # for alp2, trun in itertools.product([0.1], [1e2, 1e3]):
#     data = fastread_np(dir_1.format(0, alp2) + str(0))

#     data_sum1 = np.zeros_like(data[:, 0], dtype=complex)
#     data_sum2 = np.zeros_like(data[:, 0], dtype=complex)
#     num_line = 1
#     index = 1
#     num_line1 = 0

#     for j in range(10):
#         for i in range(0, 10000):
#             if (i % 1000 == 0):
#                 print(alp2, j, i // 1000)
#             data = fastread_np(dir_1.format(j, alp2) + str(i))
#             if (data[-1, :] < trun).all():
#                 data_sum1 += (data[:, 1]) + 1.j * data[:, 2]
#                 data_sum2 += (data[:, 7]) + 1.j * data[:, 8]
#                 num_line1 += 1
#     num_line1_l.append(num_line1)
#     np.savetxt('data_sum1_{}_{}'.format(alp2, trun), data_sum1 / num_line1)
#     np.savetxt('data_sum2_{}_{}'.format(alp2, trun), data_sum2 / num_line1)
# np.savetxt("num_line1", np.array([num_line1_l, alp2_l]).T)

dir_2 = './{}/prop-rho.dat_0.01_{}_True_'

alp2_l = [0.05, -0.05]
num_line1_g_l = []

for alp2, trun in itertools.product(alp2_l, [1e20]):
    data = fastread_np(dir_2.format(0, alp2) + str(0))
    data_sum1_g = np.zeros_like(data[:, 0], dtype=complex)
    data_sum2_g = np.zeros_like(data[:, 0], dtype=complex)
    num_line_g = 1
    index = 1
    num_line1_g = 0

    for j in range(10):
        print(j)
        for i in range(0, 10000):
            num_line1_g += 1
            data = fastread_np(dir_2.format(j, alp2) + str(i))
            data_tr_g = (data[:, 1] +
                         data[:, 7]) + 1.j * (data[:, 2] + data[:, 8])
            data_sum1_g += (data[:, 1] + 1.j * data[:, 2]) / data_tr_g
            data_sum2_g += (data[:, 7] + 1.j * data[:, 8]) / data_tr_g
    num_line1_g_l.append(num_line1_g)
    np.savetxt('data_sum1_g_{}_{}'.format(alp2, trun), data_sum1_g / num_line1_g)
    np.savetxt('data_sum2_g_{}_{}'.format(alp2, trun), data_sum2_g / num_line1_g)
np.savetxt("num_line1", np.array([num_line1_g_l, alp2_l]).T)
