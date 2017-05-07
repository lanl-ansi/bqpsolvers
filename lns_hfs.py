#!/usr/bin/env python3

### Requirements ###
# bqpjson v0.5 - pip install bqpjson
# qubo (cli) - https://github.com/alex1770/QUBO-Chimera/ - commit cf627afbf501c7028659553272611a7e313da531
# docker and the hfs_alg container are required for docker-based execution

import sys, os, argparse, json, random

from subprocess import Popen
from subprocess import PIPE

from collections import namedtuple

import bqpjson

HFS_DIR = 'hfs'
Result = namedtuple('Result', ['nodes', 'objective', 'runtime'])

# NOTE: this code assumes the HFS solver (i.e. "qubo") is in available in the local path
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

    if not os.path.exists(HFS_DIR):
        os.makedirs(HFS_DIR)

    tmp_hfs_file = '{}/tmp_{}.hfs'.format(HFS_DIR, random.randint(100000,999999))
    #TODO add file exists check and resample if needed

    print('INFO: running bqp2hfs on {}'.format(args.input_file), file=sys.stderr)
    proc = Popen(['bqp2hfs'], stdout=PIPE, stderr=PIPE, stdin=open(args.input_file, 'r'))
    stdout, stderr = proc.communicate()

    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')

    print('INFO: bqp2hfs stderr', file=sys.stderr)
    print(stderr, file=sys.stderr)

    print('INFO: bqp2hfs stdout', file=sys.stderr)
    print(stdout, file=sys.stderr)

    first_line = stdout.split('\n', 1)[0]
    chimera_degree_effective = int(first_line.split()[0])
    print('INFO: found effective chimera degree {}'.format(chimera_degree_effective), file=sys.stderr)

    scale = 1.0
    offset = 0.0
    for line in stderr.split('\n'):
        if 'scaling factor' in line:
            scale = float(line.split('scaling factor')[1].split()[0])
            print('INFO: found scaling factor {}'.format(scale), file=sys.stderr)
        if 'offset' in line:
            offset = float(line.split('offset')[1].split()[0])
            print('INFO: found offset {}'.format(offset), file=sys.stderr)


    print('INFO: writing data to {}'.format(tmp_hfs_file), file=sys.stderr)
    with open(tmp_hfs_file, 'w') as hfs_file:
        hfs_file.write(stdout)

    # print(err.getvalue())

    if args.docker_run:
        # assume that the hfs_alg container is available
        volume_map = '{}:/{}'.format(os.path.abspath(HFS_DIR), HFS_DIR)
        cmd = ['docker', 'run', '-v', volume_map, 'hfs_alg']
    else:
        # assume that the qubo executable is natively accessible
        cmd = ['qubo']

    # s - seed
    # m0 - mode of operation, try to find minimum value by heuristic search
    # N - size of Chimera graph 
    cmd.extend(['-s', '0', '-m0', '-N', str(chimera_degree_effective)])

    if args.runtime_limit != None:
        # t - min run time for some modes
        # T - max run time for some modes
        cmd.extend(['-t', str(args.runtime_limit), '-T', str(args.runtime_limit+10)])
    cmd.append(tmp_hfs_file)

    print('INFO: running command {}'.format(cmd), file=sys.stderr)
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = proc.communicate()

    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')

    print('INFO: qubo stderr', file=sys.stderr)
    print(stderr, file=sys.stderr)

    print('INFO: qubo stdout', file=sys.stderr)
    print(stdout, file=sys.stderr)

    results = []
    reading_results = False
    for line in stdout.split('\n'):
        if not reading_results:
            if 'Nodes' in line and 'bv' in line and 'nsol' in line:
                reading_results = True
        else:
            parts = line.split()
            if len(parts) == 3:
                parts = (int(parts[0]), int(parts[1]), float(parts[2]))
                results.append(Result(*parts))
            else:
                reading_results = False

    print('INFO: found {} result lines'.format(len(results)), file=sys.stderr)
    assert(len(results) > 0)

    print('INFO: removing file {}'.format(tmp_hfs_file), file=sys.stderr)
    try:
        os.remove(tmp_hfs_file)
    except:
        print('WARNING: removing file {} failed'.format(tmp_hfs_file), file=sys.stderr)

    nodes = len(data['variable_ids'])
    edges = len(data['quadratic_terms'])
    
    lt_lb = -sum(abs(lt['coeff']) for lt in data['linear_terms'])/scale
    qt_lb = -sum(abs(qt['coeff']) for qt in data['quadratic_terms'])/scale 
    lower_bound = lt_lb+qt_lb

    best_objective = results[-1].objective
    best_nodes = results[-1].nodes
    best_runtime = results[-1].runtime
    scaled_objective = scale*(best_objective+offset)
    scaled_lower_bound = scale*(lower_bound+offset)

    print('BQP_DATA, %d, %d, %f, %f, %f, %f, %f, %d, %d' % (nodes, edges, scaled_objective, scaled_lower_bound, best_objective, lower_bound, best_runtime, 0, best_nodes))


def build_cli_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input-file', help='the data file to operate on (.json)')

    parser.add_argument('-dr', '--docker-run', help='run in hfs_alg docker container', action='store_true', default=False)

    parser.add_argument('-rtl', '--runtime-limit', help='runtime limit (sec.)', type=float)

    return parser


if __name__ == '__main__':
    parser = build_cli_parser()
    main(parser.parse_args())


