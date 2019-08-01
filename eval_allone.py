#!/usr/bin/env python3

import argparse, json, random, time, math
from collections import namedtuple

import bqpjson


Model = namedtuple('Model', ['variables', 'linear', 'quadratic'])


def load_model(data):
    variables = data['variable_ids']
    linear = {lt['id']:lt['coeff'] for lt in data['linear_terms']}
    quadratic = {(qt['id_tail'],qt['id_head']):qt['coeff'] for qt in data['quadratic_terms']}
    return Model(variables, linear, quadratic)


def make_all_one_assignemnt(model):
    return {var: 1.0 for var in model.variables}


def evaluate(model, assignment):
    objective = 0.0
    for var, coeff in model.linear.items():
        objective += coeff * assignment[var]
    for (var1, var2), coeff in model.quadratic.items():
        objective += coeff * assignment[var1] * assignment[var2]
    return objective


def main(args):
    with open(args.input_file) as input_file:
        data = json.load(input_file)

    bqpjson.validate(data)

    if data['variable_domain'] != 'boolean':
        raise Exception('only boolean domains are supported. Given {}'.format(data['variable_domain']))

    model = load_model(data)
    scale, offset = data['scale'], data['offset']

    start_time = time.process_time()

    assignment = make_all_one_assignemnt(model)
    best_objective = evaluate(model, assignment)

    runtime = time.process_time() - start_time
    nodes = len(model.variables)
    edges = len(model.quadratic)
    objective = best_objective
    lower_bound = - sum(abs(lt['coeff']) for lt in data['linear_terms']) - sum(abs(qt['coeff']) for qt in data['quadratic_terms']) 
    scaled_objective = scale * (objective + offset)
    scaled_lower_bound = scale * (lower_bound + offset)
    cut_count = 0
    node_count = 1

    print('BQP_DATA, %d, %d, %f, %f, %f, %f, %f, %d, %d' % (nodes, edges, scaled_objective, scaled_lower_bound, objective, lower_bound, runtime, cut_count, node_count))


def build_cli_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input-file', help='the data file to operate on (.json)')
    return parser


if __name__ == '__main__':
    parser = build_cli_parser()
    main(parser.parse_args())
