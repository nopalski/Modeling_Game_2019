# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 20:22:45 2019

@author: Nic
"""

import tellurium as _te
import pandas as _pd
import lmfit as _lmfit
from math import sqrt as _sqrt, floor as _floor, ceil as _ceil

_mRNA_exp = None
_mRNA_exp_np = None
_prot_exp = None
_prot_exp_np = None

_unknown_model_part = """
    # mRNA translation
    => mRNA1; L1 + Vm1*(%v1%);
    => mRNA2; L2 + Vm2*(%v2%);
    => mRNA3; L3 + Vm3*(%v3%);
    => mRNA4; L4 + Vm4*(%v4%);
    => mRNA5; L5 + Vm5*(%v5%);
    => mRNA6; L6 + Vm6*(%v6%);
    => mRNA7; L7 + Vm7*(%v7%);
    => mRNA8; L8 + Vm8*(%v8%);
"""

_known_model_part = """  
    # mRNA degradation
    mRNA1 => ; d_mRNA1*mRNA1;
    mRNA2 => ; d_mRNA2*mRNA2;
    mRNA3 => ; d_mRNA3*mRNA3;
    mRNA4 => ; d_mRNA4*mRNA4;
    mRNA5 => ; d_mRNA5*mRNA5;
    mRNA6 => ; d_mRNA6*mRNA6;
    mRNA7 => ; d_mRNA7*mRNA7;
    mRNA8 => ; d_mRNA8*mRNA8;
    
    # protein translation
    => P1; a_protein1*mRNA1;
    => P2; a_protein2*mRNA2;
    => P3; a_protein3*mRNA3;
    => P4; a_protein4*mRNA4;
    => P5; a_protein5*mRNA5;
    => P6; a_protein6*mRNA6;
    => P7; a_protein7*mRNA7;
    => P8; a_protein8*mRNA8;
    
    # protein degradation
    P1 => ; d_protein1*P1;
    P2 => ; d_protein2*P2;
    P3 => ; d_protein3*P3;
    P4 => ; d_protein4*P4;
    P5 => ; d_protein5*P5;
    P6 => ; d_protein6*P6;
    P7 => ; d_protein7*P7;
    P8 => ; d_protein8*P8;
    
    # initial concentrations
    INPUT = 1;
    mRNA1 = 0;
    mRNA2 = 0;
    mRNA3 = 0;
    mRNA4 = 0;
    mRNA5 = 0;
    mRNA6 = 0;
    mRNA7 = 0;
    mRNA8 = 0;
    P1 = 0;
    P2 = 0;
    P3 = 0;
    P4 = 0;
    P5 = 0;
    P6 = 0;
    P7 = 0;
    P8 = 0;
    
    # parameters
    L1 = 0.0284569399501349;
    L2 = 0.01126375;
    L3 = 0.0108824170117358;
    L4 = 0.0142166568749038;
    L5 = 0.0124824037745626;
    L6 = 0.0166516087433626;
    L7 = 0.0178636133026099;
    L8 = 0.0165535374304033;
    
    Vm1 = 1.16380673480284;
    Vm2 = 0.855433454948057;
    Vm3 = 1.61383118795785;
    Vm4 = 0.916251896011744;
    Vm5 = 0;                                 # unknown
    Vm6 = 0.889302076744445;
    Vm7 = 0;                                 # unknown
    Vm8 = 0.87881445337468;
    
    d_mRNA1 = 0.600013;
    d_mRNA2 = 0.607263145327485;
    d_mRNA3 = 1.4096553751623;
    d_mRNA4 = 1.19069657063437;
    d_mRNA5 = 0.911653907722886;
    d_mRNA6 = 0.67287496171115;
    d_mRNA7 = 0.618711430584466;
    d_mRNA8 = 1.17483328491068;
    
    a_protein1 = 0.089105566248939;
    a_protein2 = 0.0825104648147814;
    a_protein3 = 0.118672807163739;
    a_protein4 = 0.0862964088164644;
    a_protein5 = 0.106293056264931;
    a_protein6 = 0.0890528700251159;
    a_protein7 = 0.0764169841455256;
    a_protein8 = 0.103749989801903;
    
    d_protein1 = 0.01576525;
    d_protein2 = 0.0100753359178861;
    d_protein3 = 0.0165270958726424;
    d_protein4 = 0.0205716618573404;
    d_protein5 = 0.0180685727313577;
    d_protein6 = 0.0178004316181647;
    d_protein7 = 0.0206180615545929;
    d_protein8 = 0.0131749080364666;
    
    # mRNA1 (fully known)
    H1 = 4.52340391321994;
    K1_1 = 0.0269204907071558;
    K2_1 = 0.0169635567504703;
    K3_1 = 0.0114278645720656;
    
    # mRNA2 (partially known)
    H2 = 3.21939257313515;
    K1_2 = 0.0170170903653747;
    K2_2 = 0;                                # unknown
    K3_2 = 0;                                # unknown
    
    # mRNA3 (partially known)
    H3 = 4.57189341195625;
    K1_3 = 0.0133069236136431;
    K2_3 = 0;                                # unknown
    K3_3 = 0;                                # unknown
    
    # mRNA4 (partially known)
    H4 = 5.00512303222327;
    K1_4 = 0.0179894288457716;
    K2_4 = 0;                                # unknown
    K3_4 = 0;                                # unknown
    
    # mRNA 5 (fully unknown)
    H5 = 0;                                  # unknown
    K1_5 = 0;                                # unknown
    K2_5 = 0;                                # unknown
    K3_5 = 0;
    
    # mRNA6 (fully known)
    H6 = 5.58112408673455;
    K1_6 = 0.0139445776013774;
    K2_6 = 0.0121764364668572;
    K3_6 = 0;                                # not used
    
    # mRNA7 (fully unknown)
    H7 = 0;                                  # unknown
    K1_7 = 0;                                # unknown
    K2_7 = 0;                                # unknown
    K3_7 = 0;                                # unknown
    
    # mRNA8 (partially known)
    H8 = 2.17775388441324;
    K1_8 = 0.0168599518440462;
    K2_8 = 0;                                # unknown
    K3_8 = 0;                                # unknown
"""

_sub_model = """
    # translation
    => mRNA; L + Vm*(%r%) - d_mRNA*mRNA%mRNA_num%;

    # parameters
    L = %L%;
    Vm = %Vm%;
    H = %H%;
    K1 = %K1_%;
    K2 = %K2_%;
    K3 = %K3_%;
    d_mRNA = %d_mRNA%;
                
    # proteins
    P%p1_num% = 0;
    P%p2_num% = 0;
    
    # mRNA
    mRNA%mRNA_num% = 0;
        
    # events
    %p1_events%
    %p2_events%
    %mRNA_events%
    
    # fake
    Q = 0;
"""

_known_parameters = None

_known_regulators = {
    'mRNA1': 'iA4A',
    'mRNA2': '4A',
    'mRNA3': '6A',
    'mRNA4': '2R',
    'mRNA5': None,
    'mRNA6': '1R7A',
    'mRNA7': None,
    'mRNA8': '1R'
}

_all_regs = None

_rate_laws = {
    'A': 'K1_%mRNA_num%*P%P_num%^H%mRNA_num% / (1 + K1_%mRNA_num%*P%P_num%^H%mRNA_num%)',
    'R': '1 / (1 + K1_%mRNA_num%*P%P_num%^H%mRNA_num%)',
    'AA': '(K1_%mRNA_num%*P%P1_num%^H%mRNA_num% + K2_%mRNA_num%*P%P2_num%^H%mRNA_num% + K1_%mRNA_num%*K3_%mRNA_num%*P%P1_num%^H%mRNA_num%*P%P2_num%^H%mRNA_num%) / (1 + K1_%mRNA_num%*P%P1_num%^H%mRNA_num% + K2_%mRNA_num%*P%P2_num%^H%mRNA_num% + K1_%mRNA_num%*K3_%mRNA_num%*P%P1_num%^H%mRNA_num%*P%P2_num%^H%mRNA_num%)',
    'RR': '1 / (1 + K1_%mRNA_num%*P%P1_num%^H%mRNA_num% + K2_%mRNA_num%*P%P2_num%^H%mRNA_num% + K1_%mRNA_num%*K3_%mRNA_num%*P%P1_num%^H%mRNA_num%*P%P2_num%^H%mRNA_num%)',
    # P1 is activator and P2 is repressor
    'AR': 'K1_%mRNA_num%*P%P1_num%^H%mRNA_num% / (1 + K1_%mRNA_num%*P%P1_num%^H%mRNA_num% + K2_%mRNA_num%*P%P2_num%^H%mRNA_num% + K1_%mRNA_num%*K2_%mRNA_num%*P%P1_num%^H%mRNA_num%*P%P2_num%^H%mRNA_num%)',
    # P2 is activator and P1 is repressor
    'RA': 'K1_%mRNA_num%*P%P2_num%^H%mRNA_num% / (1 + K1_%mRNA_num%*P%P2_num%^H%mRNA_num% + K2_%mRNA_num%*P%P1_num%^H%mRNA_num% + K1_%mRNA_num%*K2_%mRNA_num%*P%P1_num%^H%mRNA_num%*P%P2_num%^H%mRNA_num%)'
}

_param_ranges = {
    'V': [0.5, 2.0],
    'K': [0.01, 0.03],
    'H': [2, 8]
}

_sim_time = 1200

_num_points = 121

def get_sim_time():
    return _sim_time

def get_num_points():
    return _num_points

def set_prot_data(filename):
    global _prot_exp
    global _prot_exp_np
    _prot_exp = _pd.read_csv(filename)
    _prot_exp_np = _prot_exp.to_numpy()

def set_mRNA_data(filename):
    global _mRNA_exp
    global _mRNA_exp_np
    _mRNA_exp = _pd.read_csv(filename)
    _mRNA_exp_np = _mRNA_exp.to_numpy()

def get_sub_model():
    return _sub_model

def get_known_params():
    return _known_parameters.copy()

def get_known_regs():
    return _known_regulators.copy()

def _generate_all_regs():
    prot_nums = ['1', '2', '3', '4', '5', '6', '7', '8']
    reg_types = ['A', 'R']
    pairs = []
    for prot_num in prot_nums:
        for reg_type in reg_types:
            pairs.append(prot_num + reg_type)
    all_regs = pairs.copy()
    for i in range(len(pairs)):
        for j in range(i + 1, len(pairs)):
            if pairs[i][0] != pairs[j][0]:
                all_regs.append(pairs[i] + pairs[j])
    return all_regs

def generate_regs(known_reg):
    if known_reg is None:
        return _all_regs.copy()
    elif len(known_reg) == 4:
        return [known_reg]
    else:
        return [reg for reg in _all_regs if known_reg in reg]

def generate_mRNA_regs(mRNA_num):
    mRNA_name = 'mRNA' + str(mRNA_num)
    return generate_regs(_known_regulators[mRNA_name])
    
def _generate_events(data):
    all_events = ""
    name = data.columns[1]
    event_template = 'at time > %t%: %name% = %val%; '.replace('%name%', name)
    for i in range(1, data.shape[0] - 1):
        all_events = all_events + event_template.replace('%t%', str(data['time'][i])).replace('%val%', str(data[name][i]))
    return all_events

def _generate_prot_events(prot_num):
    if prot_num == 'i':
        return ''
    else:
        prot_name = "P" + str(prot_num)
        return _generate_events(_prot_exp[['time', prot_name]])

def _generate_mRNA_events(mRNA_num):
    mRNA_name = 'mRNA' + str(mRNA_num)
    return _generate_events(_mRNA_exp[['time', mRNA_name]])


def _make_rate(mRNA_num, regs):
    p1_num = regs[0]
    p2_num = None
    reg_type = regs[1]
    if len(regs) == 4:
        p2_num = regs[2]
        reg_type = reg_type + regs[3]
    rate_law = _rate_laws[reg_type]
    if p2_num is None:
        rate_law = rate_law.replace('%P_num%', p1_num)
    else:
        rate_law = rate_law.replace('%P1_num%', p1_num).replace('%P2_num%', p2_num)
    if mRNA_num is not None:
        rate_law = rate_law.replace('%mRNA_num%', str(mRNA_num))
    else:
        rate_law = rate_law.replace('%mRNA_num%', '').replace('_', '')
    rate_law = rate_law.replace('Pi', 'INPUT')
    return rate_law

def _generate_sub_models(mRNA_num, override_reg=None, use_mRNA_events=False):
    mRNA_num = str(mRNA_num)
    sub_models = []
    if override_reg is not None:
        reg_combos = generate_regs(override_reg)
    else:
        mRNA_name = 'mRNA' + mRNA_num
        reg_combos = generate_regs(_known_regulators[mRNA_name])
    for reg in reg_combos:
        sub_models.append(generate_sub_model(mRNA_num, reg))
    return sub_models

def generate_sub_model(mRNA_num, reg, use_mRNA_events=False):
    mRNA_num = str(mRNA_num)
    sub_model = _sub_model
    # put in rate law
    rate_law = _make_rate(None, reg)
    sub_model = sub_model.replace('%r%', rate_law)
    # put in parameters
    for param in _known_parameters:
        if param[-1] == mRNA_num:
            rep = '%' + param[0:-1] + '%'
            sub_model = sub_model.replace(rep, str(_known_parameters[param]))
    # put in proteins
    sub_model = sub_model.replace('%p1_num%', reg[0])
    if len(reg) == 4:
        sub_model = sub_model.replace('%p2_num%', reg[2])
    else:
        sub_model = sub_model.replace('%p2_num%', '')
    # put in mRNA events
    if use_mRNA_events:
        mRNA_events = _generate_mRNA_events(mRNA_num)
        sub_model = sub_model.replace('%mRNA_events%', mRNA_events)
        sub_model = sub_model.replace('%mRNA_num%', mRNA_num)
    else:
        sub_model = sub_model.replace('%mRNA_events%', '').replace('%mRNA_num%', '')
    # put in protein events
    p1_events = _generate_prot_events(reg[0])
    sub_model = sub_model.replace('%p1_events%', p1_events)
    if len(reg) == 4:
        p2_events = _generate_prot_events(reg[2])
        sub_model = sub_model.replace('%p2_events%', p2_events)
    else:
        sub_model = sub_model.replace('%p2_events%', '')
    sub_model = sub_model.replace('P = 0;', '')
    return sub_model

def _min(parameters, rr, mRNA_num=None, selections=None):
    fitter = _lmfit.Minimizer(_make_objective_fnc(rr, mRNA_num=mRNA_num, selections=selections), parameters)
    fit = fitter.minimize(method='differential_evolution')
    print('took', fit.nfev, ' function evals')
    return fit.params.valuesdict(), fit.chisqr

def min_full(ant_model, selections=None):
    rr = _te.loada(ant_model)
    parameters = _make_parameters_full(ant_model)
    return _min(parameters, rr, selections=selections)

def min_sub(ant_model, mRNA_num):
    rr = _te.loada(ant_model)
    parameters = _make_parameters_sub(ant_model, rr)
    return _min(parameters, rr, mRNA_num=mRNA_num)
    
def _make_parameters(ant_model, param_dict):
    parameters = _lmfit.Parameters()
    for param in param_dict:
        # all unknown parameters are initialized to zero
        # all parameters being used will be present at least twice (once in an initialization and at least once in the rate law)
        start = param[0]
        # print(param + ': ' + str(param_dict[param] == 0) + ', ' + str(ant_model.count(param) > 1) + ', ' + str(start in _param_ranges.keys()))
        if param_dict[param] == 0 and ant_model.count(param) > 1 and start in _param_ranges.keys():
            parameters.add(param, min=_param_ranges[start][0], max=_param_ranges[start][1])
    if (len(parameters.valuesdict())) == 0:
        parameters.add('Q', min=0, max=0.0000001)
    #print(parameters.valuesdict().keys())
    return parameters

def _make_parameters_full(ant_model):
    return _make_parameters(ant_model, _known_parameters)
    
def _make_parameters_sub(ant_model, rr_model):
    param_dict = dict(zip(rr_model.getGlobalParameterIds(), rr_model.getGlobalParameterValues()))
    return _make_parameters(ant_model, param_dict)

def _make_objective_fnc(rr, mRNA_num=None, selections=None):
    if selections is not None:
        rr.timeCourseSelections = selections
    def objective_fnc(p):
        rr.reset()
        params = p.valuesdict()
        for param in params:
            rr[param] = params[param]
        sim_data = rr.simulate(0, _sim_time, _num_points)
        # col_names = [name.replace('[', '').replace(']','') for name in sim_data.colnames]
        # col_names = [name[1:-1] for name in sim_data.colnames]
        # print(sim_data.colnames)
        # sim_data = _pd.DataFrame(sim_data, columns=col_names)
        sum_squares = 0
        if mRNA_num is None:    
            for col in sim_data.colnames:                
                if col[1] == 'P':
                    idx = int(col[-2])
                    sum_squares = sum_squares + ((_prot_exp_np[:,idx] - sim_data[col]) ** 2).sum()
                elif col[1] == 'm':
                    idx = int(col[-2])
                    sum_squares = sum_squares + ((_mRNA_exp_np[:,idx] - sim_data[col]) ** 2).sum()
        else:
            # mRNA_name = 'mRNA' + str(mRNA_num)
            sum_squares = ((_mRNA_exp_np[:,mRNA_num] - sim_data['[mRNA]']) ** 2).sum()
        return sum_squares
    return objective_fnc

def get_dims(num):
    sq = _sqrt(num)
    cols = _ceil(sq)
    rows = _floor(sq)
    if cols*rows < num:
        rows = rows + 1
    return (rows, cols)

def make_model_combos(guesses, cur_combos=[[]], cur_mRNA_num=1):
    if cur_mRNA_num <= 8:
        mRNA_name = 'mRNA' + str(cur_mRNA_num)
        new_combos = []
        for combo in cur_combos:
            for guess in guesses[mRNA_name]:
                combo_copy = combo.copy()
                combo_copy.append(guess)
                new_combos.append(combo_copy)
        return make_model_combos(guesses, new_combos, cur_mRNA_num + 1)
    else:
        return cur_combos

def make_model(regs):
    unknown_part = _unknown_model_part
    for i in range(len(regs)):
        rate = _make_rate(i + 1, regs[i])
        unknown_part = unknown_part.replace('%v' + str(i + 1) + '%', rate)
    return unknown_part + _known_model_part

def _init():
    global _known_parameters
    rr = _te.loada(_known_model_part)
    _known_parameters = dict(zip(rr.model.getGlobalParameterIds(), rr.model.getGlobalParameterValues()))
    
    global _all_regs
    _all_regs = _generate_all_regs()

_init()
    
    
