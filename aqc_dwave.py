#!/usr/bin/env python2

### Requirements ###
# bqpjson v0.5 - pip install bqpjson
# dwave_sapi2 v2.6 - available in the qubist interface

import argparse, json, time, os, sys

from dwave_sapi2.remote import RemoteConnection
from dwave_sapi2.core import solve_qubo

import bqpjson

DEFAULT_CONFIG_FILE = '_config'

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

    # A core assumption of this solver is that the given B-QP will magically be compatable with the given D-Wave QPU
    dw_url = args.dw_url
    dw_token = args.dw_token
    dw_solver_name = args.dw_solver_name
    dw_chip_id = None

    print(args.config_file)
    print(dw_token)

    if 'dw_url' in data['metadata']:
        dw_url = data['metadata']['dw_url'].encode('ascii','ignore')
        print('using d-wave url provided in data file: %s' % dw_url)

    if 'dw_solver_name' in data['metadata']:
        dw_solver_name = data['metadata']['dw_solver_name'].encode('ascii','ignore')
        print('using d-wave solver name provided in data file: %s' % dw_solver_name)

    if 'dw_chip_id' in data['metadata']:
        dw_chip_id = data['metadata']['dw_chip_id'].encode('ascii','ignore')
        print('found d-wave chip id in data file: %s' % dw_chip_id)

    if dw_url is None or dw_token is None or dw_solver_name is None:
        print('dwave solver parameters not found')
        quit()

    if args.dw_proxy is None: 
        remote_connection = RemoteConnection(dw_url, dw_token)
    else:
        remote_connection = RemoteConnection(dw_url, dw_token, args.dw_proxy)

    solver = remote_connection.get_solver(dw_solver_name)
    if not dw_chip_id is None:
        if solver.properties['chip_id'] != dw_chip_id:
            print('WARNING: chip ids do not match.  data: %s  hardware: %s' % (dw_chip_id, solver.properties['chip_id']))

    couplers = solver.properties['couplers']
    sites = solver.properties['qubits']

    site_range = tuple(solver.properties['h_range'])
    coupler_range = tuple(solver.properties['j_range'])

    Q = {}
    #obj = data['offset']
    for lt in data['linear_terms']:
        i = lt['id']
        assert(not (i,i) in Q)
        Q[(i,i)] = lt['coeff']

    for qt in data['quadratic_terms']:
        i = qt['id_tail']
        j = qt['id_head']
        assert(not (i,j) in Q)
        Q[(i,j)] = qt['coeff']

    params = {
        'auto_scale': False,
        'num_reads': args.num_reads,
        'num_spin_reversal_transforms': args.num_reads/args.spin_reversal_transform_rate,
        'annealing_time': args.annealing_time
    }

    print('dwave parameters:')
    for k,v in params.items():
        print('  {} - {}'.format(k,v))

    t0 = time.time()
    answers = solve_qubo(solver, Q, **params)
    solve_time = time.time() - t0

    for i in range(len(answers['energies'])):
        print('%f - %d' % (answers['energies'][i], answers['num_occurrences'][i]))

    nodes = len(data['variable_ids'])
    edges = len(data['quadratic_terms'])
    
    lt_lb = -sum(abs(lt['coeff']) for lt in data['linear_terms'])
    qt_lb = -sum(abs(qt['coeff']) for qt in data['quadratic_terms']) 
    lower_bound = lt_lb+qt_lb

    best_objective = answers['energies'][0]
    best_nodes = args.num_reads
    best_runtime = answers['timing']['total_real_time']/1000000.0
    scaled_objective = data['scale']*(best_objective+data['offset'])
    scaled_lower_bound = data['scale']*(lower_bound+data['offset'])

    print('BQP_DATA, %d, %d, %f, %f, %f, %f, %f, %d, %d' % (nodes, edges, scaled_objective, scaled_lower_bound, best_objective, lower_bound, best_runtime, 0, best_nodes))


# loads a configuration file and sets up undefined CLI arguments
def load_config(args):
    config_file_path = args.config_file

    if os.path.isfile(config_file_path):
        with open(config_file_path, 'r') as config_file:
            try:
                config_data = json.load(config_file)
                for key, value in config_data.items():
                    if isinstance(value, dict) or isinstance(value, list):
                        print('invalid value for configuration key "%s", only single values are allowed' % config_file_path)
                        quit()
                    if not hasattr(args, key) or getattr(args, key) == None:
                        if isinstance(value, unicode):
                            value = value.encode('ascii','ignore')
                        setattr(args, key, value)
                    else:
                        print('skipping the configuration key "%s", it already has a value of %s' % (key, str(getattr(args, key))))
                    #print(key, value)
            except ValueError:
                print('the config file does not appear to be a valid json document: %s' % config_file_path)
                quit()
    else:
        if config_file_path != DEFAULT_CONFIG_FILE:
            print('unable to open conifguration file: %s' % config_file_path)
            quit()

    return args


def build_cli_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input-file', help='the data file to operate on (.json)')

    parser.add_argument('-cf', '--config-file', help='a configuration file for specifing common parameters', default=DEFAULT_CONFIG_FILE)

    parser.add_argument('-url', '--dw-url', help='url of the d-wave machine')
    parser.add_argument('-token', '--dw-token', help='token for accessing the d-wave machine')
    parser.add_argument('-proxy', '--dw-proxy', help='proxy for accessing the d-wave machine')
    parser.add_argument('-solver', '--dw-solver-name', help='d-wave solver to use', type=int)

    parser.add_argument('-nr', '--num-reads', help='the number of reads to take from the d-wave hardware', type=int, default=10000)
    parser.add_argument('-at', '--annealing-time', help='the annealing time of each d-wave sample', type=int, default=5)
    parser.add_argument('-srtr', '--spin-reversal-transform-rate', help='the number of reads to take before each spin reversal transform', type=int, default=100)

    return parser


if __name__ == '__main__':
    parser = build_cli_parser()
    main(load_config(parser.parse_args()))





