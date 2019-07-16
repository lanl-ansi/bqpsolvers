#!/usr/bin/env python3

import argparse, json, random, time, math
from collections import namedtuple

import bqpjson


Model = namedtuple('Model', ['variables', 'linear', 'quadratic', 'adjacent'])


def load_model(data):
    variable = data['variable_ids']
    linear = {lt['id']:lt['coeff'] for lt in data['linear_terms']}
    quadratic = {(qt['id_tail'],qt['id_head']):qt['coeff'] for qt in data['quadratic_terms']}
    adjacent = {}
    for qt in data['quadratic_terms']:
        i, j, coeff = qt['id_tail'], qt['id_head'], qt['coeff']
        adjacent.setdefault(i, []).append((j, coeff))
        adjacent.setdefault(j, []).append((i, coeff))
    return Model(variable, linear, quadratic, adjacent)


def make_random_assignemnt(model):
    return {var: random.choice([0.0,1.0]) for var in model.variables}


def evaluate(model, assignment):
    objective = 0.0
    for var, coeff in model.linear.items():
        objective += coeff * assignment[var]
    for (var1, var2), coeff in model.quadratic.items():
        objective += coeff * assignment[var1] * assignment[var2]
    return objective


def flip(assignment, variable):
    assignment[variable] = 1.0 - assignment[variable]


def flip_delta(model, assignment, variable):
    delta = 0.0
    difference = 1.0 - 2.0 * assignment[variable]
    if variable in model.linear:
        delta += model.linear[variable] * difference
    for i, coeff in model.adjacent[variable]:
        delta += coeff * assignment[i] * difference
    return delta


def step(model, assignment, objective):
    best_var = None
    best_delta = math.inf
    for var in model.variables:
        delta = flip_delta(model, assignment, var)
        if delta < best_delta:
            best_var = var
            best_delta = delta
    if best_delta >= 0.0:
        return None
    else:
        flip(assignment, best_var)
        return best_delta


def main(args):
    with open(args.input_file) as input_file:
        data = json.load(input_file)

    bqpjson.validate(data)

    if data['variable_domain'] != 'boolean':
        raise Exception('only boolean domains are supported. Given {}'.format(data['variable_domain']))

    model = load_model(data)
    scale, offset = data['scale'], data['offset']

    assignment = make_random_assignemnt(model)
    objective = evaluate(model, assignment)
    iterations = 1
    restarts = 0
    best_objective = math.inf
    start_time = time.process_time()
    end_time = start_time + args.runtime_limit

    while time.process_time() < end_time:
        result = step(model, assignment, objective)
        if result is None: # restart
            assignment = make_random_assignemnt(model)
            objective = evaluate(model, assignment)
            restarts += 1
        else: # move downward
            objective += result
            iterations += 1
        if objective < best_objective:
            best_objective = objective
        if args.show_objectives:
            print('objective:',  objective)
        if args.show_scaled_objectives:
            print('scaled objective:', scale * (objective + offset))

    runtime = time.process_time() - start_time
    nodes = len(model.variables)
    edges = len(model.quadratic)
    objective = best_objective
    lower_bound = - sum(abs(lt['coeff']) for lt in data['linear_terms']) - sum(abs(qt['coeff']) for qt in data['quadratic_terms']) 
    scaled_objective = scale * (objective + offset)
    scaled_lower_bound = scale * (lower_bound + offset)
    cut_count = 0
    node_count = iterations

    print()
    print('iterations:', iterations)
    print('restarts:', restarts)
    print('best objective:', objective)
    print('best scaled objective:', scaled_objective)

    print()
    print('BQP_DATA, %d, %d, %f, %f, %f, %f, %f, %d, %d' % (nodes, edges, scaled_objective, scaled_lower_bound, objective, lower_bound, runtime, cut_count, node_count))


def build_cli_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input-file', help='the data file to operate on (.json)')
    parser.add_argument('-so', '--show-objectives', help='print the objectives seen by the program', action='store_true', default=False)
    parser.add_argument('-sso', '--show-scaled-objectives', help='print the scaled objectives seen by the program', action='store_true', default=False)
    parser.add_argument('-rtl', '--runtime-limit', help='runtime limit (sec.)', type=float, default=10)
    return parser


if __name__ == '__main__':
    parser = build_cli_parser()
    main(parser.parse_args())
