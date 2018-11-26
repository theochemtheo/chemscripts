#!/usr/bin/python3

'''
This script takes the output from a pysoc calculation and
calculates the mixing parameter between states from
2nd order perturbation theory
'''

import numpy as np
import argparse
import parseFUNCS as pF


def lambda_print(S, T):
    mixing = lambdas[S, T]
    if mixing > 0.1:
        print('S{}-T{}:\t{:.6e}\tWarning, lambda is large, perturbation theory is probably invalid here!'.format(S + 1, T + 1, mixing))
    else:
        print('S{}-T{}:\t{:.6e}'.format(S + 1, T + 1, mixing))


parser = argparse.ArgumentParser(description="This script calculates mixing parameters")
specify = parser.add_mutually_exclusive_group()
specify.add_argument("-s", dest="state", metavar="N", help="Only print the Nth state", required=False, type=int)
specify.add_argument("-r", dest="range", metavar="N", help="Print up to the Nth state", required=False, type=int)
args = vars(parser.parse_args())


ene_out = np.genfromtxt('./ene_out.dat')[:, 1:]
soc_out = pF.cmtoeV(np.genfromtxt('./soc_out.dat', usecols=3))
n_states = ene_out.shape[-1]
lambdas = np.zeros([n_states, n_states])
for S in range(n_states):
    for T in range(n_states):
        deltaE_ST = np.abs(ene_out[0, S] - ene_out[1, T])
        soc_ST = soc_out[n_states * (S + 1) + T]
        lambdas[S, T] = (soc_ST / deltaE_ST)**2  # 2nd order perturbation theory


if args["state"]:
    for T in range(n_states):
        lambda_print(args["state"] - 1, T)
elif args["range"]:
    for S in range(args["range"]):
        for T in range(n_states):
            lambda_print(S, T)
else:
    for S in range(n_states):
        for T in range(n_states):
            lambda_print(S, T)
