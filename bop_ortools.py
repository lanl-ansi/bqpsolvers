#!/usr/bin/env python3

### Requirements ###
# bqpjson v0.5 - pip install bqpjson
# ortools v1.5 - https://developers.google.com/optimization/ - installed and available by calling qubo

import sys, json, argparse

from ortools.linear_solver import pywraplp

import bqpjson

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

    solver = pywraplp.Solver('BOP', pywraplp.Solver.BOP_INTEGER_PROGRAMMING)
    solver.EnableOutput()
    #solver.SetSolverSpecificParametersAsString('prune_search_tree:true')
    #solver.SetSolverSpecificParametersAsString('use_random_lns:false')
    #solver.SetSolverSpecificParametersAsString('num_relaxed_vars:40')
    #solver.SetSolverSpecificParametersAsString('use_potential_one_flip_repairs_in_ls:true')
    #solver.SetSolverSpecificParametersAsString('use_lp_strong_branching:true')
    #solver.SetSolverSpecificParametersAsString('lp_max_deterministic_time:10.0')

    if args.runtime_limit != None:
        solver.SetTimeLimit(args.runtime_limit*1000)

    variable_ids = set(data['variable_ids'])
    variable_product_ids = set([(qt['id_tail'], qt['id_head']) for qt in data['quadratic_terms']])

    variable_lookup = {}
    for vid in variable_ids:
        variable_lookup[(vid,vid)] = solver.BoolVar(name='site_{:04d}'.format(vid))
    for pair in variable_product_ids:
        variable_lookup[pair] = solver.BoolVar(name='product_{:04d}_{:04d}'.format(*pair))

    # models conjunction of two binary variablies
    for i,j in variable_product_ids:
        solver.Add(variable_lookup[(i,j)] >= variable_lookup[(i,i)] + variable_lookup[(j,j)] - 1)
        solver.Add(variable_lookup[(i,j)] <= variable_lookup[(i,i)])
        solver.Add(variable_lookup[(i,j)] <= variable_lookup[(j,j)])
        # TODO is there a way to give "/\" to the solver directly?

    linear_terms = [int(lt['coeff'])*variable_lookup[(lt['id'], lt['id'])] for lt in data['linear_terms']]
    quadratic_terms = [int(qt['coeff'])*variable_lookup[(qt['id_tail'], qt['id_head'])] for qt in data['quadratic_terms']]
    obj_expr = solver.Sum(linear_terms + quadratic_terms)

    obj = solver.Minimize(obj_expr)

    solver.Solve()

    if args.show_solution:
        print('')
        for k,v in variable_lookup.items():
            print('{} - {}'.format(k, v.SolutionValue()))

    print('')
    print('obj_ub = {}'.format(solver.Objective().Value()))
    print('obj_lb = {}'.format(solver.Objective().BestBound()))
    print('walltime: {}ms'.format(solver.WallTime()))

    nodes = len(data['variable_ids'])
    edges = len(data['quadratic_terms'])

    obj_ub = solver.Objective().Value()
    obj_lb = solver.Objective().BestBound()

    node_count = 0 # iterations and nodes are not available
    cut_count = 0
    runtime = solver.WallTime()/1000.0

    scaled_objective = data['scale']*(obj_ub+data['offset'])
    scaled_lower_bound = data['scale']*(obj_lb+data['offset'])
    print('BQP_DATA, %d, %d, %f, %f, %f, %f, %f, %d, %d' % (nodes, edges, scaled_objective, scaled_lower_bound, obj_ub, obj_lb, runtime, cut_count, node_count))


def build_cli_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input-file', help='the data file to operate on (.json)')

    parser.add_argument('-rtl', '--runtime-limit', help='runtime limit (sec.)', type=int)
    parser.add_argument('-tl', '--thread-limit', help='thread limit', type=int, default=10)
    parser.add_argument('-ss', '--show-solution', help='prints the best solutoin found at termination', action='store_true', default=False)

    return parser


if __name__ == '__main__':
    parser = build_cli_parser()
    main(parser.parse_args())
