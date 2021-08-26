import math
import re
import os
import time
import numpy as np
import sympy as sp
import itertools
import json
import matplotlib.pyplot as plt


def convert(o):
    if isinstance(o, np.int64): return int(o)
    raise TypeError


def fBose(x, pole, resi):
    return 1 / x + 0.5 + sum(2.0 * resi[i] * x / (x**2 + pole[i]**2)
                             for i in range(len(pole)))


def tseig(D, E):
    mat = np.diag(E, -1) + np.diag(D, 0) + np.diag(E, 1)
    return -np.sort(-np.linalg.eigvalsh(mat))


def MSD(N, BoseFermi=1):
    if BoseFermi == 1:
        pole = np.array([2 * (i + 1) * np.pi for i in range(N)])
        resi = np.ones(N, dtype=float)
        return pole, resi
    elif BoseFermi == 2:
        pole = np.array([(2 * i + 1) * np.pi for i in range(N)])
        resi = np.ones(N, dtype=float)
        return pole, resi


def PSD(N, BoseFermi=1, pade=1):
    if N < 0 or BoseFermi < 1 or BoseFermi > 2 or pade < 0 or pade > 3:
        raise ValueError("N or BoseFermi or pade has wrong value!")

    if pade == 0:
        return MSD(N, BoseFermi)
    elif pade == 1 or pade == 2:
        pole, resi = [], []
        if N > 0:
            M = 2 * N + pade // 2
            temp = 3.0 if BoseFermi == 1 else 1.0
            diag = np.zeros(M, dtype=float)
            doff = np.array([
                1.0 / math.sqrt((temp + 2.0 * i) * (temp + 2.0 * (i + 1)))
                for i in range(M - 1)
            ])
            pole = 2.0 / tseig(diag, doff)[:N]
            pol2 = np.array([x * x for x in pole])
            M -= 1
            temp = 5.0 if BoseFermi == 1 else 3.0
            diag = np.zeros(M, dtype=float)
            doff = np.array([
                1.0 / math.sqrt((temp + 2.0 * i) * (temp + 2.0 * (i + 1)))
                for i in range(M - 1)
            ])
            M //= 2
            eig2 = np.power(2.0 / tseig(diag, doff)[:M], 2)
            scaling = 0.0
            if BoseFermi == 1:
                scaling = N*(2.0*N+3.0) if pade == 1 else 1.0 / \
                    (4.0*(N+1.0)*(2.0*N+3.0))
            elif BoseFermi == 2:
                scaling = N*(2.0*N+1.0) if pade == 1 else 1.0 / \
                    (4.0*(N+1.0)*(2.0*N+1.0))
            resi = np.zeros(N, dtype=float)
            for j in range(N):
                if pade == 2:
                    temp = 0.5 * scaling * (eig2[j] - pol2[j])
                elif pade == 1:
                    if j == N - 1:
                        temp = 0.5 * scaling
                    else:
                        temp = 0.5*scaling * \
                            (eig2[j]-pol2[j])/(pol2[N-1]-pol2[j])
                for k in range(M):
                    temp *= (eig2[k]-pol2[j]) / \
                        (pol2[k]-pol2[j]) if k != j else 1.0
                resi[j] = temp
        rn, tn = 0.0, 0.0
        if BoseFermi == 1 and pade == 2:
            rn = 1.0 / (4.0 * (N + 1.0) * (2.0 * N + 3.0))
        return pole, resi
    elif pade == 3:
        Np1 = N + 1
        temp = 3.0 if BoseFermi == 1 else 1.0
        d = np.empty(2 * Np1, dtype=float)
        d[0] = 0.25 / temp
        d[-1] = -4.0 * (N + 1.0) * (N + 1.0) * (temp + 2 * N) * (
            temp + 2 * N) * (temp + 4 * N + 2.0)
        for i in range(1, Np1):
            d[2*i-1] = -4.0*i*i*(temp+2.0*i-2.0) * \
                (temp+2.0*i-2.0)*(temp+4.0*i-2.0)
            d[2 * i] = -0.25 * (temp + 4.0 * i) / i / (i + 1) / (
                temp + 2.0 * i - 2.0) / (temp + 2.0 * i)
        sumd2 = np.empty(Np1, dtype=float)
        sumd2[0] = d[1]
        for i in range(1, Np1):
            sumd2[i] = sumd2[i - 1] + d[2 * i + 1]
        tn = 0.25 / sumd2[-1]
        rn = sum(d[2 * i] * (4.0 * tn *
                             (sumd2[-1] - sumd2[i - 1]))**2 if i > 0 else d[2 *
                                                                            i]
                 for i in range(Np1))
        M = 2 * N + 1
        diag = np.zeros(M, dtype=float)
        doff = np.array(
            [1.0 / math.sqrt(d[i + 1] * d[i + 2]) for i in range(M - 1)])
        pole = 2.0 / tseig(diag, doff)[:N]
        resi = np.zeros(N, dtype=float)
        for j in range(N):
            scaling = pole[j] * pole[j]
            r0, t1 = 0.0, 0.25 / d[1]
            eta0, eta1, eta2 = 0.0, 0.5, 0.0
            for i in range(Np1):
                r1 = t1 if (i == j
                            or i == N) else t1 / (pole[i] * pole[i] - scaling)
                r2 = 2.0*math.sqrt(abs(r1)) if r1 > 0 else - \
                    2.0*math.sqrt(abs(r1))
                r1 = 2.0 * math.sqrt(abs(r1))
                eta2 = d[2 * i] * r1 * eta1 - 0.25 * r1 * r0 * scaling * eta0
                eta0 = eta1
                eta1 = eta2
                eta2 = d[2 * i +
                         1] * r2 * eta1 - 0.25 * r2 * r1 * scaling * eta0
                eta0 = eta1
                eta1 = eta2
                r0 = r2
                if i != N:
                    t1 = sumd2[i] / sumd2[i + 1]
            resi[j] = eta2
        return pole, resi


def arma_print(ndarray):

    shape = ndarray.shape
    dimen = len(shape)

    if dimen == 1:

        if issubclass(type(ndarray[0]), np.int_):
            print('ARMA_MAT_TXT_IS004\n%d %d' % (shape[0], 1))
            for row in ndarray:
                print('%d' % row)
        elif issubclass(type(ndarray[0]), float):
            print('ARMA_MAT_TXT_FN008\n%d %d' % (shape[0], 1))
            for row in ndarray:
                print('%.8e' % row)
        elif issubclass(type(ndarray[0]), complex):
            print('ARMA_MAT_TXT_FC016\n%d %d' % (shape[0], 1))
            for row in ndarray:
                print('(%.8e,%-.8e)' % (row.real, row.imag))

    elif dimen == 2:

        if issubclass(type(ndarray[0, 0]), np.int_):
            print('ARMA_MAT_TXT_IS004\n%d %d' % (shape[0], shape[1]))
            for row in ndarray:
                print(' '.join('%d' % x for x in row))
        elif issubclass(type(ndarray[0, 0]), float):
            print('ARMA_MAT_TXT_FN008\n%d %d' % (shape[0], shape[1]))
            for row in ndarray:
                print(' '.join('%.8e' % x for x in row))
        elif issubclass(type(ndarray[0, 0]), complex):
            print('ARMA_MAT_TXT_FC016\n%d %d' % (shape[0], shape[1]))
            for row in ndarray:
                print(' '.join('(%.8e,%-.8e)' % (x.real, x.imag) for x in row))

    elif dimen == 3:

        if issubclass(type(ndarray[0, 0, 0]), np.int_):
            print('ARMA_CUB_TXT_IS004\n%d %d %d' %
                  (shape[1], shape[2], shape[0]))
            for slc in ndarray:
                for row in slc:
                    print(' '.join('%d' % x for x in row))
        elif issubclass(type(ndarray[0, 0, 0]), float):
            print('ARMA_CUB_TXT_FN008\n%d %d %d' %
                  (shape[1], shape[2], shape[0]))
            for slc in ndarray:
                for row in slc:
                    print(' '.join('%-.8e' % x for x in row))
        elif issubclass(type(ndarray[0, 0, 0]), complex):
            print('ARMA_CUB_TXT_FC016\n%d %d %d' %
                  (shape[1], shape[2], shape[0]))
            for slc in ndarray:
                for row in slc:
                    print(' '.join('(%.8e,%-.8e)' % (x.real, x.imag)
                                   for x in row))


def arma_write(ndarray, filename):

    shape = ndarray.shape
    dimen = len(shape)

    with open(filename, 'w') as f:
        if dimen == 1:
            if issubclass(type(ndarray[0]), np.int_):
                print('ARMA_MAT_TXT_IS004\n%d %d' % (shape[0], 1), file=f)
                for row in ndarray:
                    print('%d' % row, file=f)
            elif issubclass(type(ndarray[0]), float):
                print('ARMA_MAT_TXT_FN008\n%d %d' % (shape[0], 1), file=f)
                for row in ndarray:
                    print('%.8e' % row, file=f)
            elif issubclass(type(ndarray[0]), complex):
                print('ARMA_MAT_TXT_FC016\n%d %d' % (shape[0], 1), file=f)
                for row in ndarray:
                    print('(%.8e,%-.8e)' % (row.real, row.imag), file=f)

        elif dimen == 2:

            if issubclass(type(ndarray[0, 0]), np.int_):
                print('ARMA_MAT_TXT_IS004\n%d %d' % (shape[0], shape[1]),
                      file=f)
                for row in ndarray:
                    print(' '.join('%d' % x for x in row), file=f)
            elif issubclass(type(ndarray[0, 0]), float):
                print('ARMA_MAT_TXT_FN008\n%d %d' % (shape[0], shape[1]),
                      file=f)
                for row in ndarray:
                    print(' '.join('%.8e' % x for x in row), file=f)
            elif issubclass(type(ndarray[0, 0]), complex):
                print('ARMA_MAT_TXT_FC016\n%d %d' % (shape[0], shape[1]),
                      file=f)
                for row in ndarray:
                    print(' '.join('(%.8e,%-.8e)' % (x.real, x.imag)
                                   for x in row),
                          file=f)

        elif dimen == 3:

            if issubclass(type(ndarray[0, 0, 0]), np.int_):
                print('ARMA_CUB_TXT_IS004\n%d %d %d' %
                      (shape[1], shape[2], shape[0]),
                      file=f)
                for slc in ndarray:
                    for row in slc:
                        print(' '.join('%d' % x for x in row))
            elif issubclass(type(ndarray[0, 0, 0]), float):
                print('ARMA_CUB_TXT_FN008\n%d %d %d' %
                      (shape[1], shape[2], shape[0]),
                      file=f)
                for slc in ndarray:
                    for row in slc:
                        print(' '.join('%-.8e' % x for x in row), file=f)
            elif issubclass(type(ndarray[0, 0, 0]), complex):
                print('ARMA_CUB_TXT_FC016\n%d %d %d' %
                      (shape[1], shape[2], shape[0]),
                      file=f)
                for slc in ndarray:
                    for row in slc:
                        print(' '.join('(%.8e,%-.8e)' % (x.real, x.imag)
                                       for x in row),
                              file=f)


# in this script, we can decompose any given spectrum, but the sympy format is must been given
# do u like haskell?
# sympy[spe(def by sympy)], dict[sp_para_dict], dict[para_dict], dict[npsd],
# dict[pade] >> np.array[etal], np.array[etar],np.array[etaa], np.array[expn]


def decompose_spe(spe, sp_para_dict, para_dict, condition_dict, npsd, pade=1):
    numer, denom = sp.cancel(
        sp.factor(sp.cancel(
            spe.subs(condition_dict)).as_real_imag()[1])).as_numer_denom()
    numer_get_para = (sp.factor(numer)).subs(sp_para_dict)
    denom_get_para = (sp.factor(denom)).subs(sp_para_dict)

    poles = sp.nroots(denom_get_para)
    float(sp.re(poles[0]))

    expn = []
    poles_allplane = np.array([])
    for i in poles:
        i = complex(i)
        if i.imag < 0:
            expn.append(i * 1.J)
        poles_allplane = np.append(poles_allplane, i)

    etal = []
    etar = []
    etaa = []

    expn = np.array(expn)

    expn_imag_sort = np.argsort(np.abs(np.imag(expn)))[::-1]
    expn_imag = np.sort(np.abs(np.imag(expn)))[::-1]

    expn_val_cc = expn[expn_imag_sort[expn_imag != 0]]
    # expn_arg_cc = expn_imag_sort[expn_imag != 0]
    expn_val_n_cc = expn[expn_imag_sort[expn_imag == 0]]
    # expn_arg_n_cc = expn_imag_sort[expn_imag == 0]

    expn = list(expn[expn_imag_sort])
    pole, resi = PSD(npsd, 1, pade)
    beta = para_dict['beta']
    temp = 1 / beta

    for ii in range(0, len(expn_val_cc), 2):
        etal.append(
            complex(
                sp.N((-2.j * numer_get_para /
                      np.multiply.reduce(w_sp - poles_allplane[np.abs(
                          poles_allplane + 1.J * expn_val_cc[ii]) > 1e-14])
                      ).subs({w_sp: -1.j * expn_val_cc[ii]}) *
                     fBose(-1.J * expn_val_cc[ii] / temp, pole, resi))))

        etal.append(
            complex(
                sp.N((-2.j * numer_get_para /
                      np.multiply.reduce(w_sp - poles_allplane[np.abs(
                          poles_allplane + 1.J * expn_val_cc[ii + 1]) > 1e-14])
                      ).subs({w_sp: -1.j * expn_val_cc[ii + 1]}) *
                     fBose(-1.J * expn_val_cc[ii + 1] / temp, pole, resi))))

        etar.append(np.conj(etal[-1]))
        etar.append(np.conj(etal[-2]))
        etaa.append(np.sqrt(np.abs(etal[-2]) * np.abs(etar[-2])))
        etaa.append(np.sqrt(np.abs(etal[-1]) * np.abs(etar[-1])))

    for ii in range(len(expn_val_n_cc)):
        etal.append(
            complex(
                sp.N((-2.j * numer_get_para /
                      np.multiply.reduce(w_sp - poles_allplane[np.abs(
                          poles_allplane + 1.J * expn_val_n_cc[ii]) > 1e-14])
                      ).subs({w_sp: -1.j * expn_val_n_cc[ii]}) *
                     fBose(-1.J * expn_val_n_cc[ii] / temp, pole, resi))))
        etar.append(np.conj(etal[-1]))
        etaa.append(np.sqrt(np.abs(etal[-1]) * np.abs(etar[-1])))

    f = numer_get_para / np.multiply.reduce(w_sp - poles_allplane)
    f = sp.lambdify(w_sp, f)

    for inma in range(len(pole)):
        zomg = -1.J * pole[inma] * temp
        jsum = np.sum(f(zomg))
        expn.append(pole[inma] * temp)
        etal.append(-2.J * resi[inma] * temp * jsum)
        etar.append(np.conj(etal[-1]))
        etaa.append(np.abs(etal[-1]))

    etal = np.array(etal)
    etar = np.array(etar)
    etaa = np.array(etaa)
    expn = np.array(expn)

    #     print('%%%%%%#######%%%%%%')
    #     print(etal)
    #     print(etar)
    #     print(etaa)
    #     print(expn)
    arma_write(expn, 'inp_expn.mat')
    arma_write(etal, 'inp_etal.mat')
    arma_write(etar, 'inp_etar.mat')
    arma_write(etaa, 'inp_etaa.mat')

    #     arma_print(expn)
    #     arma_print(etal)
    #     arma_print(etar)
    #     arma_print(etaa)

    return etal, etar, etaa, expn


def fit_J(w, expn, etal):
    res = np.zeros_like(w, dtype=complex)
    for i in range(len(etal)):
        res += etal[i] / (expn[i] - 1.j * w)
    return res


def INDEX3(i, j, k, mum):
    return mum * mum * i + mum * j + k


# for vib_key, S_0, npsd, lmax in itertools.product([0, 0.5, 1.0], [0.25, 0.5], [1, 2], [40, 50, 60]):
for vib_a, S_0, npsd, lmax, dt in itertools.product([0], [0.1], [2], [1],
                                                    [0.01]):
    # for vib_a, S_0, npsd, lmax in itertools.product([0, 1], [0.1, 0.2, 0.3, 0.4], [1, 2], [40, 50]):
    nmod = 1
    omgs = 1.0  # in 200cm-1
    # Deltesip = 100  # in 20000cm-1
    Deltesip = 50  # in 10000cm-1
    temp = 300 * 0.69503457 / 200  # in 298K
    gams = temp
    beta = 1 / temp
    lamdxx = 0.1  # in 20cm-1
    D2 = 2.0 * S_0
    # etas = omgs * 2.0
    # lams = D2 / 2.0
    # lamd0 = omgs * D2 / 2.0
    g1 = 0.5
    g2 = 0.5
    D_val = np.sqrt(2 * S_0)
    lamdxy = D_val * lamdxx / g1
    lamdyy = D_val * lamdxy / g2
    # eta = 1.0
    w_sp, lamd1_sp, gams1_sp, lamd2_sp, gams2_sp, lamd3_sp, gams3_sp, D_sp, omgs_sp, beta_sp = sp.symbols(
        r"\omega , \lambda_1, \gamma, \lambda_2, \gamma, \lambda_3, \gamma, D, \Omega_{s}, \beta",
        real=True)

    phixx_sp = 2 * lamd1_sp * gams1_sp / (gams1_sp - sp.I * w_sp)
    phixy_sp = 2 * lamd2_sp * gams2_sp / (gams2_sp - sp.I * w_sp)
    phiyy_sp = 2 * lamd3_sp * gams3_sp / (gams3_sp - sp.I * w_sp)

    barD_sp = phixy_sp.subs({w_sp: 0
                             }) - (omgs_sp + phixx_sp.subs({w_sp: 0})) * D_sp
    xqq_sp = omgs_sp / (omgs_sp * omgs_sp - w_sp * w_sp - omgs_sp *
                        (phixx_sp - phixx_sp.subs({w_sp: 0})))

    # spe_vib_sp = phixx_sp
    # spe_vib_b_sp = phixy_sp
    # spe_b_sp = phiyy_sp
    # spe_sp = spe_vib_sp + spe_vib_b_sp + spe_vib_b_sp + spe_b_sp

    spe_vib_sp = barD_sp * barD_sp * xqq_sp
    spe_vib_b_sp = 1 / barD_sp * (D_sp * phixx_sp - phixy_sp) * spe_vib_sp
    spe_b_sp = D_sp * D_sp * phixx_sp - 2 * D_sp * phixy_sp + phiyy_sp + \
        1 / barD_sp * (D_sp * phixx_sp - phixy_sp) * spe_vib_b_sp
    spe_sp = D_sp*D_sp*phixx_sp-2*D_sp*phixy_sp+phiyy_sp + \
        (barD_sp+D_sp*phixx_sp-phixy_sp) * \
        (barD_sp+D_sp*phixx_sp-phixy_sp)*xqq_sp

    sp_para_dict = {
        lamd1_sp: lamdxx,
        lamd2_sp: lamdxx,
        lamd3_sp: lamdxx,
        gams1_sp: gams,
        gams2_sp: gams,
        gams3_sp: gams,
        D_sp: D_val,
        omgs_sp: omgs
    }

    condition_dict = {
        lamd2_sp: lamd1_sp,
        lamd3_sp: lamd1_sp,
    }
    para_dict = {'beta': beta}

    etal1, etar1, etaa1, expn1 = decompose_spe(spe_vib_sp, sp_para_dict,
                                               para_dict, condition_dict, npsd)
    etal2, etar2, etaa2, expn2 = decompose_spe(spe_vib_b_sp, sp_para_dict,
                                               para_dict, condition_dict, npsd)
    etal3, etar3, etaa3, expn3 = decompose_spe(spe_b_sp, sp_para_dict,
                                               para_dict, condition_dict, npsd)

    expn = expn1

    for i in np.append(expn2, expn3):
        if ((len(np.argwhere(np.abs(i - expn) < 1e-8))) == 0):
            expn = np.append(expn, i)
            # print(i, expn, np.argwhere(np.abs(i - expn) < 1e-14),
            #   expn[np.argwhere(np.abs(i - expn) < 1e-14)[0][0]]-i)
            # print(expn)
    # print(expn)

    n_spe = 4

    etar = np.zeros((n_spe, len(expn)), dtype=complex)
    etal = np.zeros((n_spe, len(expn)), dtype=complex)
    etaa = np.zeros((n_spe, len(expn)), dtype=float)

    expn_all = [expn1, expn2, expn2, expn3]
    etar_all = [etar1, etar2, etar2, etar3]
    etal_all = [etal1, etal2, etal2, etal3]
    etaa_all = [etaa1, etaa2, etaa2, etaa3]

    for i in range(n_spe):
        for j in range(len(expn)):
            arg = np.argwhere(np.abs(expn[j] - expn_all[i]) < 1e-8)
            if len(arg):
                etar[i, j] = etar_all[i][arg[0]]
                etal[i, j] = etal_all[i][arg[0]]
                etaa[i, j] = etaa_all[i][arg[0]]
            else:
                etar[i, j] = 0
                etal[i, j] = 0
                etaa[i, j] = 0

    n_not_zero = 0
    etar_n = np.zeros((n_spe, len(expn)), dtype=complex)
    etal_n = np.zeros((n_spe, len(expn)), dtype=complex)
    etaa_n = np.zeros((n_spe, len(expn)), dtype=float)
    expn_n = np.zeros(len(expn), dtype=complex)
    for i in range(len(expn)):
        is_zero = False
        for j in range(n_spe):
            is_zero = is_zero or (np.abs(etaa[j, i]) > 1e-8)
        if (is_zero):
            for j in range(n_spe):
                etar_n[j, n_not_zero] = etar[j, i]
                etal_n[j, n_not_zero] = etal[j, i]
                etaa_n[j, n_not_zero] = etaa[j, i]
                expn_n[n_not_zero] = expn[i]
            n_not_zero = n_not_zero + 1

    etar = np.ndarray.flatten(np.transpose(etar_n[:, :n_not_zero]))
    etal = np.ndarray.flatten(np.transpose(etal_n[:, :n_not_zero]))
    etaa = np.sum(etaa_n[:, :n_not_zero], axis=0)
    expn = expn_n[:n_not_zero]
    mode = np.append(np.zeros_like(expn, dtype=int),
                     1 + np.zeros_like(expn, dtype=int))

    etar = np.append(etar, etar)
    etal = np.append(etal, etal)
    etaa = np.append(etaa, etaa)
    expn = np.append(expn, expn)

    arma_write(etar, 'inp_etar.mat')
    arma_write(etal, 'inp_etal.mat')
    arma_write(etaa, 'inp_etaa.mat')
    arma_write(mode, 'inp_mode.mat')
    arma_write(expn, 'inp_expn.mat')

    arma_write(
        etar, 'inp_etar_{}_{}_{}_{}.mat'.format(str(S_0), str(vib_key),
                                                str(lmax), str(npsd)))
    arma_write(
        etal, 'inp_etal_{}_{}_{}_{}.mat'.format(str(S_0), str(vib_key),
                                                str(lmax), str(npsd)))
    arma_write(
        etaa, 'inp_etaa_{}_{}_{}_{}.mat'.format(str(S_0), str(vib_key),
                                                str(lmax), str(npsd)))
    arma_write(
        mode, 'inp_mode_{}_{}_{}_{}.mat'.format(str(S_0), str(vib_key),
                                                str(lmax), str(npsd)))
    arma_write(
        expn, 'inp_expn_{}_{}_{}_{}.mat'.format(str(S_0), str(vib_key),
                                                str(lmax), str(npsd)))

    lamd1, gams1 = sp_para_dict[lamd1_sp], sp_para_dict[gams1_sp]
    lamd2, gams2 = sp_para_dict[lamd2_sp], sp_para_dict[gams2_sp]
    lamd3, gams3 = sp_para_dict[lamd3_sp], sp_para_dict[gams3_sp]
    beta = para_dict['beta']
    D = sp_para_dict[D_sp]
    omgs = sp_para_dict[omgs_sp]
    len_ = 10000
    spe_wid = 100
    w = np.append(np.linspace(0, spe_wid, len_),
                  np.linspace(-spe_wid, 0, len_))

    phixx = 2 * lamd1 * gams1 / (gams1 - 1.j * w)
    phixy = 2 * lamd2 * gams2 / (gams2 - 1.j * w)
    phiyy = 2 * lamd3 * gams3 / (gams3 - 1.j * w)

    barD = phixy[0] - (omgs + phixx[0]) * D
    xqq = omgs / (omgs * omgs - w * w - omgs * (phixx - phixx[0]))

    spe_vib = barD * barD * xqq
    spe_vib_b = 1 / barD * (D * phixx - phixy) * spe_vib
    spe_b = D * D * phixx - 2 * D * phixy + phiyy + 1 / barD * (
        D * phixx - phixy) * spe_vib_b
    # print(sp_para_dict)
    # print(npsd, barD, phixx[0])
    # arma_print(etar)
    # arma_print(etal)
    # arma_print(etaa)
    # arma_print(mode)
    # arma_print(expn)

    etal_check_0 = np.zeros_like(expn, dtype=complex)
    etal_check_1 = np.zeros_like(expn, dtype=complex)
    etal_check_3 = np.zeros_like(expn, dtype=complex)

    for i in range(len(expn)):
        etal_check_0[i] = etal[INDEX3(i, 0, 0, 2)]
        etal_check_1[i] = etal[INDEX3(i, 1, 0, 2)]
        etal_check_3[i] = etal[INDEX3(i, 1, 1, 2)]

    plt.plot(w[:len_], (np.imag(spe_vib) / (1 - np.exp(-beta * w)))[:len_],
             'r',
             label='spe_vib')
    plt.plot(w[len_:], (np.imag(spe_vib) / (1 - np.exp(-beta * w)))[len_:],
             'r')
    plt.plot(w[:len_],
             fit_J(w, expn, etal_check_0).real[:len_] / 2,
             '-.b',
             label='spe_vib')
    plt.plot(w[len_:], fit_J(w, expn, etal_check_0).real[len_:] / 2, '-.b')
    plt.plot(w[:len_],
             fit_J(w, expn1, etal1).real[:len_],
             ':k',
             label='spe_vib')
    plt.plot(w[len_:], fit_J(w, expn1, etal1).real[len_:], ':k')
    plt.xlim(-2, 2)

    plt.plot(w[:len_], (np.imag(spe_vib_b) / (1 - np.exp(-beta * w)))[:len_],
             'g',
             label='spe_vib_b')
    plt.plot(w[len_:], (np.imag(spe_vib_b) / (1 - np.exp(-beta * w)))[len_:],
             'g')
    plt.plot(w[:len_],
             fit_J(w, expn, etal_check_1).real[:len_] / 2,
             '-.b',
             label='spe_vib_b')
    plt.plot(w[len_:], fit_J(w, expn, etal_check_1).real[len_:] / 2, '-.b')
    plt.plot(w[:len_],
             fit_J(w, expn2, etal2).real[:len_],
             ':k',
             label='spe_vib_b')
    plt.plot(w[len_:], fit_J(w, expn2, etal2).real[len_:], ':k')
    plt.xlim(-2, 2)

    plt.plot(w[:len_], (np.imag(spe_b) / (1 - np.exp(-beta * w)))[:len_],
             'y',
             label='spe_b')
    plt.plot(w[len_:], (np.imag(spe_b) / (1 - np.exp(-beta * w)))[len_:], 'y')
    plt.plot(w[:len_],
             fit_J(w, expn, etal_check_3).real[:len_] / 2,
             '-.b',
             label='spe_b')
    plt.plot(w[len_:], fit_J(w, expn, etal_check_3).real[len_:] / 2, '-.b')
    plt.plot(w[:len_], fit_J(w, expn3, etal3).real[:len_], ':k', label='spe_b')
    plt.plot(w[len_:], fit_J(w, expn3, etal3).real[len_:], ':k')
    plt.xlim(-2, 2)

    plt.legend(loc='best')
    plt.savefig('spe_{}_{}_{}_{}.pdf'.format(str(S_0), str(vib_a), str(lmax),
                                             str(npsd)))
    plt.clf()
    # spe_vib = phixx
    # spe_vib_b = phixy
    # spe_b = phiyy

    spe = spe_vib + spe_vib_b + spe_vib_b + spe_b

    # syst
    print(spe[0], D2 / 2)
    lamd = 1 / 2 * spe[0]  # reorganization energy
    hams = np.zeros((2, 2), dtype=complex)
    hams[0, 0] = (2) / 2
    hams[1, 1] = -(2) / 2

    qmds = np.zeros((nmod, 2, 2), dtype=complex)
    qmds[0, 1, 0] = 1.0
    qmds[0, 0, 1] = 1.0

    # hidx
    trun = 0
    nmax = 100000
    ferr = 1e-18

    # arma_write(hame, jsonInit['td-rhot']['hameFile'])

    sdip = np.zeros((2, 2), dtype=float)
    sdip[0, 1] = sdip[1, 0] = 1.0
    arma_write(sdip, 'inp_sdip.mat')
    bdip = 10 * np.ones(12 + 2 * npsd, dtype=float)
    arma_write(bdip, 'inp_bdip.mat')

    pdip = np.zeros((nmod, 2, 2), dtype=float)
    pdip[0, 0, 0] = pdip[0, 1, 1] = 1.0
    arma_write(pdip, 'inp_pdip.mat')

    # proprho
    jsonInit = {
        "nk": 1,
        "nmax": nmax,
        "lmax": lmax,
        "n_threads": 8,
        "nsys": 2,
        "nham": 2,
        "ferr": ferr,
        "ti": -np.pi * 80,
        "tf": np.pi * 80,
        "nind": len(expn1),
        "nmod": len(expn1),
        "nmodmax": nmod,
        "dt": dt,
        "nt": 20,
        "staticErr": 0,
        "inistate": 0,
        "expn": {
            "real": list(np.real(expn1.flatten())),
            "imag": list(np.imag(expn1.flatten()))
        },
        "ham1": {
            "real": list(np.real(hams.flatten())),
            "imag": list(np.imag(hams.flatten()))
        },
        "qmd1": {
            "real": list(np.real(qmds.flatten())),
            "imag": list(np.imag(qmds.flatten()))
        },
        "modLabel": {
            "real": list(np.real(mode.flatten())),
            "imag": list(np.imag(mode.flatten()))
        },
        "delt_res": {
            "real": list(np.zeros(nmod)),
            "imag": list(np.zeros(nmod))
        },
        "coef_abs": {
            "real": list(np.real(etaa1.flatten())),
            "imag": list(np.imag(etaa1.flatten()))
        },
        "coef_lft": {
            "real": list(np.real(etal1.flatten())),
            "imag": list(np.imag(etal1.flatten()))
        },
        "coef_rht": {
            "real": list(np.real(etar1.flatten())),
            "imag": list(np.imag(etar1.flatten()))
        },
        "sdip_cub": {
            "real": list(np.real(sdip.flatten())),
            "imag": list(np.imag(sdip.flatten()))
        },
        "bdip_cub": {
            "real": list(np.real(bdip.flatten())),
            "imag": list(np.imag(bdip.flatten()))
        },
        "pdip_cub": {
            "real": list(np.real(pdip.flatten())),
            "imag": list(np.imag(pdip.flatten()))
        },
        "ampl": 1,
        "freq": 2,
        "sigm": 1.5 * np.pi,
        "varp": 0.0
    }

    with open('input.json', 'w') as f:
        json.dump(jsonInit, f, indent=4, default=convert)