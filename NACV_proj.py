#!/usr/bin/python3

import numpy as np
from numpy import linalg as LA


NACVs = {}
NACVs['T1toT2'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T1-to-T2.NACV.mat')
NACVs['T1toT3'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T1-to-T3.NACV.mat')
NACVs['T1toT4'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T1-to-T4.NACV.mat')
NACVs['T1toT5'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T1-to-T5.NACV.mat')
NACVs['T1toT6'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T1-to-T6.NACV.mat')
NACVs['T2toT3'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T2-to-T3.NACV.mat')
NACVs['T2toT4'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T2-to-T4.NACV.mat')
NACVs['T2toT5'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T2-to-T5.NACV.mat')
NACVs['T2toT6'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T2-to-T6.NACV.mat')
NACVs['T3toT4'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T3-to-T4.NACV.mat')
NACVs['T3toT5'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T3-to-T5.NACV.mat')
NACVs['T3toT6'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T3-to-T6.NACV.mat')
NACVs['T4toT5'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T4-to-T5.NACV.mat')
NACVs['T4toT6'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T4-to-T6.NACV.mat')
NACVs['T5toT6'] = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.T5-to-T6.NACV.mat')

nNACVs = {}
for key in NACVs.keys():
    nNACVs[key] = np.divide(NACVs[key], LA.norm(NACVs[key], 'fro'))

natoms = len(NACVs['T1toT2'])

freqlist = np.genfromtxt('B3LYP_CRENBL-Pt-6311Gdp-HCNOPS_CPCM_S0_NAPme_CC_Pt2PBu3_CC_Ph_CH2_PTZ_freq-hpmodes.freqlist')

freqvecs = {}
print("Reading modes")
for i in range(len(freqlist)):
    thismode = np.int(freqlist[i][0])
    freqvecs[thismode] = np.genfromtxt('B3LYP_CRENBL-Pt-6-311Gdp-HCNOPS_CPCM_S0_NAPme_CC_Pt2PBu3_CC_Ph_CH2_PTZ_freq-hpmodes.m{}.dat'.format(thismode))
print("Modes read")

for pairs in NACVs.keys():
    similarities = {}
    for mode in freqvecs.keys():
        sim = 0
        for atom in range(natoms):
            sim += np.abs(np.dot(nNACVs[pairs][atom, :], freqvecs[mode][atom, :]))
        similarities[mode] = sim
    combined = np.column_stack((freqlist, np.array(list(similarities.values()))))
    np.savetxt('B3LYP_CRENBL-Pt-6311Gdp_CPCM_TDA6_T_S0_MeNAP-CC-PtPBu3-CC-Ph-CH2-PTZ_NACME-6_c10+t14_175-974.{}.NACV.scaled.dist'.format(pairs), combined, delimiter='\t')
