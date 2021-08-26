import math
import re
import os
from sys import flags
import time
import numpy as np
import sympy as sp
import itertools
import json
import matplotlib.pyplot as plt
from scipy.linalg import sqrtm


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
    numer, denom = sp.cancel(sp.factor(sp.cancel(
        spe.subs(condition_dict)))).as_numer_denom()
    numer_get_para = (sp.factor(numer)).subs(sp_para_dict)
    denom_get_para = (sp.factor(denom)).subs(sp_para_dict)
    print(numer_get_para, "$$$$$$", denom_get_para)
    poles = sp.nroots(denom_get_para)
    float(sp.re(poles[0]))
    print(poles)

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
for girsanov, npsd, lmax, (alp1, alp2), dt in itertools.product([False, True], [1],
                                                                [1000],
                                                                [(1, 0.1),
                                                                 (1, -0.1),
                                                                 (1, 0.05),
                                                                 (1, -0.05)],
                                                                [0.01]):
    # for vib_a, S_0, npsd, lmax in itertools.product([0, 1], [0.1, 0.2, 0.3, 0.4], [1, 2], [40, 50]):
    nmod = 1
    temp = 1  # in 298K
    beta = 1 / temp
    omgs = 1
    zeta = 4  # in 20cm-1
    w_sp, zeta_sp, omgs_sp, beta_sp = sp.symbols(
        r"\omega , \zeta_{b}, \omega_{b}, \beta", real=True)

    phixx_sp = w_sp * omgs_sp * zeta_sp / (
        (omgs_sp * omgs_sp - w_sp * w_sp) *
        (omgs_sp * omgs_sp - w_sp * w_sp) + zeta_sp * zeta_sp * w_sp * w_sp)

    spe_vib_sp = phixx_sp

    sp_para_dict = {zeta_sp: zeta, omgs_sp: omgs}

    condition_dict = {}
    para_dict = {'beta': beta}

    etal, etar, etaa, expn = decompose_spe(spe_vib_sp, sp_para_dict, para_dict,
                                           condition_dict, npsd)
    mode = np.zeros_like(expn, dtype=int)

    beta = para_dict['beta']

    len_ = 10000
    spe_wid = 100
    w = np.append(np.linspace(0, spe_wid, len_),
                  np.linspace(-spe_wid, 0, len_))

    phixx = w * omgs * zeta / ((omgs * omgs - w * w) *
                               (omgs * omgs - w * w) + zeta * zeta * w * w)

    spe_vib = phixx

    plt.plot(w[:len_], ((spe_vib) / (1 - np.exp(-beta * w)))[:len_],
             'r',
             label='spe_vib')
    plt.plot(w[len_:], ((spe_vib) / (1 - np.exp(-beta * w)))[len_:],
             'r')
    plt.plot(w[:len_], fit_J(w, expn, etal).real[:len_], ':k', label='spe_vib')
    plt.plot(w[len_:], fit_J(w, expn, etal).real[len_:], ':k')
    plt.xlim(-20, 50)
    plt.legend(loc='best')
    plt.savefig('spe_{}_{}.pdf'.format(str(lmax), str(npsd)))
    plt.clf()

    # syst
    hams = np.zeros((2, 2), dtype=complex)
    hams[0, 0] = -1
    hams[1, 1] = 1
    hams[0, 1] = 1
    hams[1, 0] = 1

    qmds = np.zeros((nmod, 2, 2), dtype=complex)
    qmds[0, 0, 0] = 0
    qmds[0, 0, 1] = 0
    qmds[0, 1, 0] = 0
    qmds[0, 1, 1] = 1

    # hidx
    trun = 0
    nmax = 5000000
    ferr = 1e-14

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
        "ti": 0,
        "tf": 20,
        "nind": len(expn),
        "nmod": len(expn),
        "nmodmax": nmod,
        "dt": dt,
        "nt": 1,
        "staticErr": 0,
        "inistate": 0,
        "alp1": alp1,
        "alp2": alp2,
        "expn": {
            "real": list(np.real(expn.flatten())),
            "imag": list(np.imag(expn.flatten()))
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
            "real": list(np.real(etaa.flatten())),
            "imag": list(np.imag(etaa.flatten()))
        },
        "coef_lft": {
            "real": list(np.real(etal.flatten())),
            "imag": list(np.imag(etal.flatten()))
        },
        "coef_rht": {
            "real": list(np.real(etar.flatten())),
            "imag": list(np.imag(etar.flatten()))
        },
        "frho": "",
        "girsanov": girsanov
    }

    for i in range(10000):
        jsonInit["frho"] = 'prop-rho.dat_' + str(dt) + '_' + str(alp2) + '_' + str(girsanov) + '_' + str(i)
        with open('input.json', 'w') as f:
            json.dump(jsonInit, f, indent=4, default=convert)
        cmd = r'. ./r.sh'
        os.system(cmd)
