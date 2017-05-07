#!/usr/bin/env python3

### Requirements ###
# bqpjson v0.5 - pip install bqpjson
# gurobi v7.0 - http://www.gurobi.com/

import argparse, json, sys

from gurobipy import *

import bqpjson


# tolerance for determining integrability of gurobi variables
integrality_tol = 1e-6

def main(args):
    if args.input_file == None:
        data = json.load(sys.stdin)
    else:
        with open(args.input_file) as file:
            data = json.load(file)

    bqpjson.validate(data)

    if data['variable_domain'] != 'boolean':
        print('only boolean domains are supported. Given %s' % data['variable_domain'])
        quit()

    variable_ids = set(data['variable_ids'])
    variable_product_ids = set([(qt['id_tail'], qt['id_head']) for qt in data['quadratic_terms']])


    m = Model()

    if args.runtime_limit != None:
        m.setParam('TimeLimit', args.runtime_limit)

    m.setParam('Threads', args.thread_limit)

    if args.cuts != None:
        m.setParam('Cuts', args.cuts)

    #m.setParam('Cuts', 3)
    #m.setParam('MIPFocus', 1)
    #m.setParam('MIPFocus', 2)

    variable_lookup = {}
    for vid in variable_ids:
        variable_lookup[(vid,vid)] = m.addVar(lb=0, ub=1, vtype = GRB.BINARY, name='site_%04d' % vid)
    for pair in variable_product_ids:
        variable_lookup[pair] = m.addVar(lb=0, ub=1, vtype = GRB.BINARY, name='product_%04d_%04d' % (pair[0], pair[1]))
    m.update()

    for i,j in variable_product_ids:
        #m.addConstr(variable_lookup[(i,i)]*variable_lookup[(j,j)] >= variable_lookup[(i,j)]*variable_lookup[(i,j)])
        m.addConstr(variable_lookup[(i,j)] >= variable_lookup[(i,i)] + variable_lookup[(j,j)] - 1)
        m.addConstr(variable_lookup[(i,j)] <= variable_lookup[(i,i)])
        m.addConstr(variable_lookup[(i,j)] <= variable_lookup[(j,j)])
        #m.addGenConstrAnd(variable_lookup[(i,j)], [variable_lookup[(i,i)], variable_lookup[(j,j)]])
        
    obj = 0.0
    for lt in data['linear_terms']:
        i = lt['id']
        obj += lt['coeff']*variable_lookup[(i,i)]

    for qt in data['quadratic_terms']:
        i = qt['id_tail']
        j = qt['id_head']
        obj += qt['coeff']*variable_lookup[(i,j)]

    m.setObjective(obj, GRB.MINIMIZE)

    m.update()

    m._cut_count = 0
    m.optimize(cut_counter)

    if False:
        print('')
        for k,v in variable_lookup.items():
            print('%s - %f' % (v.VarName, v.X))

    if False:
        print('')

        if 'solutions' in data:
            sol = data['solutions'][0]
            print(sol['evaluation'])
        # for v in variable_ids:
        #     g_var = variable_lookup[(v,v)]
        #     print('%s - %f' % (g_var.VarName, g_var.X))


    lower_bound = m.MIPGap*m.ObjVal + m.ObjVal
    scaled_objective = data['scale']*(m.ObjVal+data['offset'])
    scaled_lower_bound = data['scale']*(lower_bound+data['offset'])
    print('BQP_DATA, %d, %d, %f, %f, %f, %f, %f, %d, %d' % (len(variable_ids), len(variable_product_ids), scaled_objective, scaled_lower_bound, m.ObjVal, lower_bound, m.Runtime, m._cut_count, m.NodeCount))


def cut_counter(model, where):
    cut_names = {
        'Clique:', 'Cover:', 'Flow cover:', 'Flow path:', 'Gomory:', 
        'GUB cover:', 'Inf proof:', 'Implied bound:', 'Lazy constraints:', 
        'Learned:', 'MIR:', 'Mod-K:', 'Network:', 'Projected Implied bound:', 
        'StrongCG:', 'User:', 'Zero half:'}
    if where == GRB.Callback.MESSAGE:
        # Message callback
        msg = model.cbGet(GRB.Callback.MSG_STRING)
        if any(name in msg for name in cut_names):
            model._cut_count += int(msg.split(':')[1])


def build_cli_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input-file', help='the data file to operate on (.json)')

    parser.add_argument('-rtl', '--runtime-limit', help='gurobi runtime limit (sec.)', type=int)
    parser.add_argument('-tl', '--thread-limit', help='gurobi thread limit', type=int, default=1)
    parser.add_argument('-cuts', help='gurobi cuts parameter', type=int)

    return parser


if __name__ == '__main__':
    parser = build_cli_parser()
    main(parser.parse_args())



