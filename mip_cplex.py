#!/usr/bin/env python3

### Requirements ###
# bqpjson v0.5 - pip install bqpjson
# cplex v12.7.0.0 - https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
#

### Note ###
# these are good articles to reference when using this solver
#
# @article{1306.1202,
#   Author = {Sanjeeb Dash},
#   Title = {A note on QUBO instances defined on Chimera graphs},
#   Year = {2013},
#   Eprint = {arXiv:1306.1202},
#   url = {https://arxiv.org/abs/1612.05024}
# }
#
# @Article{Billionnet2007,
#   author="Billionnet, Alain and Elloumi, Sourour",
#   title="Using a Mixed Integer Quadratic Programming Solver for the Unconstrained Quadratic 0-1 Problem",
#   journal="Mathematical Programming",
#   year="2007",
#   volume="109",
#   number="1",
#   pages="55--68",
#   issn="1436-4646",
#   doi="10.1007/s10107-005-0637-9",
#   url="http://dx.doi.org/10.1007/s10107-005-0637-9"
# }
#

import argparse, json, sys

import cplex
from cplex.exceptions import CplexSolverError
from cplex.callbacks import MIPInfoCallback

import bqpjson


class StatsCallback(MIPInfoCallback):
    def __call__(self):
        if not hasattr(self, 'cut_types'):
            self.cut_types = [
                MIPInfoCallback.cut_type.GUB_cover,
                MIPInfoCallback.cut_type.MIR,
                MIPInfoCallback.cut_type.clique,
                MIPInfoCallback.cut_type.cover,
                MIPInfoCallback.cut_type.disjunctive,
                MIPInfoCallback.cut_type.flow_cover,
                MIPInfoCallback.cut_type.flow_path,
                MIPInfoCallback.cut_type.fractional,
                MIPInfoCallback.cut_type.implied_bound,
                MIPInfoCallback.cut_type.lift_and_project,
                MIPInfoCallback.cut_type.multi_commodity_flow,
                MIPInfoCallback.cut_type.solution_pool,
                MIPInfoCallback.cut_type.table,
                MIPInfoCallback.cut_type.user,
                MIPInfoCallback.cut_type.zero_half
            ]

        self.cuts = sum(self.get_num_cuts(ct) for ct in self.cut_types)
        self.nodes = self.get_num_nodes()


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

    coefficient_lookup = {}
    for lt in data['linear_terms']:
        i = lt['id']
        coefficient_lookup[(i,i)] = lt['coeff']
    for qt in data['quadratic_terms']:
        i = qt['id_tail']
        j = qt['id_head']
        coefficient_lookup[(i,j)] = qt['coeff']

    m = cplex.Cplex()
    stats_cb = m.register_callback(StatsCallback)

    if args.runtime_limit != None:
        m.parameters.timelimit.set(args.runtime_limit)

    m.parameters.threads.set(args.thread_limit)

    variable_lookup = {}
    for vid in variable_ids:
        coeff = 0
        if (vid,vid) in coefficient_lookup:
            coeff = coefficient_lookup[(vid,vid)]

        indexes = m.variables.add(obj=[coeff], lb=[0], ub=[1], types=['B'])
        assert(len(indexes) == 1)
        variable_lookup[(vid,vid)] = next((x for x in indexes))

    for pair in variable_product_ids:
        indexes = variable_lookup[pair] = m.variables.add(
            obj=[coefficient_lookup[pair]],
            lb=[0], ub=[1], types=["B"])
        assert(len(indexes) == 1)
        variable_lookup[pair] = next((x for x in indexes))

    for i,j in variable_product_ids:
        terms = cplex.SparsePair(ind=[variable_lookup[(i,j)], variable_lookup[(i,i)], variable_lookup[(j,j)]], val=[1, -1, -1])
        m.linear_constraints.add(lin_expr=[terms], senses=['G'], rhs=[-1])

        terms = cplex.SparsePair(ind=[variable_lookup[(i,j)], variable_lookup[(i,i)]], val=[1, -1])
        m.linear_constraints.add(lin_expr=[terms], senses=['L'], rhs=[0])

        terms = cplex.SparsePair(ind=[variable_lookup[(i,j)], variable_lookup[(j,j)]], val=[1, -1])
        m.linear_constraints.add(lin_expr=[terms], senses=['L'], rhs=[0])

    spin_data = bqpjson.core.swap_variable_domain(data)
    if len(spin_data['linear_terms']) <= 0 or all(lt['coeff'] == 0.0 for lt in spin_data['linear_terms']):
        print('detected spin symmetry, adding symmetry breaking constraint')
        v1 = data['variable_ids'][0]
        terms = cplex.SparsePair(ind=[variable_lookup[(v1,v1)]], val=[1])
        m.linear_constraints.add(lin_expr=[terms], senses=['E'], rhs=[0])

    m.objective.set_sense(m.objective.sense.minimize)

    try:
        timestart = m.get_time()
        m.solve()
        runtime = m.get_time() - timestart
    except CplexSolverError as e:
        print("Exception raised during solve: " + e)
    else:
        solution = m.solution
        #print(solution.get_quality_metrics())
        print('status: {} - {}'.format(solution.get_status(), solution.get_status_string()))
        print('   gap: {}'.format(solution.MIP.get_mip_relative_gap()))
        print('  cuts: {}'.format(stats_cb.cuts))
        print(' nodes: {}'.format(stats_cb.nodes))

        if args.show_solution:
            print('')
            for k in sorted(variable_lookup):
                v = variable_lookup[k]
                print('{}: {}'.format(k, solution.get_values(v)))

        upper_bound = solution.get_objective_value()
        lower_bound = upper_bound*(1+solution.MIP.get_mip_relative_gap())
        scaled_upper_bound = data['scale']*(upper_bound+data['offset'])
        scaled_lower_bound = data['scale']*(lower_bound+data['offset'])

        cut_count = stats_cb.cuts
        node_count = stats_cb.nodes

        print('')
        print('BQP_DATA, %d, %d, %f, %f, %f, %f, %f, %d, %d' % (len(variable_ids), len(variable_product_ids), scaled_upper_bound, scaled_lower_bound, upper_bound, lower_bound, runtime, cut_count, node_count))


def build_cli_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input-file', help='the data file to operate on (.json)')

    parser.add_argument('-rtl', '--runtime-limit', help='cplex runtime limit (sec.)', type=int)
    parser.add_argument('-tl', '--thread-limit', help='cplex thread limit', type=int, default=1)
    parser.add_argument('-ss', '--show-solution', help='print the solution', action='store_true', default=False)

    return parser


if __name__ == '__main__':
    parser = build_cli_parser()
    main(parser.parse_args())

