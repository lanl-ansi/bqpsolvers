#!/usr/bin/env python3

### Requirements ###
# bqpjson v0.5 - pip install bqpjson
# gurobi v7.0 - http://www.gurobi.com/
#

import argparse, json, sys

from gurobipy import *

import bqpjson


class UnionFind:
    def __init__(self, n):
        self._parent = list(range(n))

    def union(self, a, b):
        a = self.find(a)
        b = self.find(b)
        self._parent[a] = b

    def find(self, a):
        while a != self._parent[a]:
            self._parent[a] = self._parent[self._parent[a]]
            a = self._parent[a]
        return a


def ordered(i, j):
    return (i, j) if i < j else (j, i)


class Ising:
    def __init__(self, data):
        assert data['variable_domain'] == 'spin'
        self.variables = set(data['variable_ids'])
        self.linear = {lt['id']: lt['coeff'] for lt in data['linear_terms']}
        self.quadratic = {(qt['id_tail'],qt['id_head']):qt['coeff'] for qt in data['quadratic_terms']}
        self.adjacent = {}
        for i, j in self.quadratic:
            self.adjacent.setdefault(i, set()).add(j)
            self.adjacent.setdefault(j, set()).add(i)
        self.scale = data['scale']
        self.offset = data['offset']

    def merge(self, var1, var2):
        self.linear[var1] += self.linear[var2]
        del self.linear[var2]
        if ordered(var1, var2) in self.quadratic:
            self.offset += self.quadratic[ordered(var1, var2)]
            del self.quadratic[ordered(var1, var2)]
        self.adjacent[var2].discard(var1)
        for v in self.adjacent[var2]:
            self.adjacent[v].remove(var2)
            self.adjacent[v].add(var1)
            self.quadratic[ordered(var1, v)] = self.quadratic.get(ordered(var1, v), 0) + self.quadratic[ordered(var2, v)]
            del self.quadratic[ordered(v, var2)]
        self.adjacent[var1] |= self.adjacent[var2]
        del self.adjacent[var2]
        self.variables.remove(var2)

    def to_boolean_domain(self):
        offset = self.offset
        linear = {}
        quadratic = {}
        for var, coeff in self.linear.items():
            linear[var] = 2.0 * coeff
            offset -= coeff
        for (i, j), coeff in self.quadratic.items():
            if (i, j) not in quadratic:
                quadratic[i, j] = 0.0
            if i not in linear:
                linear[i] = 0.0
                linear[j] = 0.0
            quadratic[i, j] = quadratic.get((i, j), 0.0) + 4.0 * coeff
            linear[i] = linear.get(i, 0.0) - 2.0 * coeff
            linear[j] = linear.get(j, 0.0) - 2.0 * coeff
            offset += coeff
        return self.variables, linear, quadratic, self.scale, offset


def main(args):
    if args.input_file == None:
        data = json.load(sys.stdin)
    else:
        with open(args.input_file) as input_file:
            data = json.load(input_file)

    bqpjson.validate(data)

    if data['variable_domain'] != 'boolean':
        print('only boolean domains are supported. Given %s' % data['variable_domain'])
        quit()

    ising = Ising(bqpjson.swap_variable_domain(data))
    uf = UnionFind(max(ising.variables) + 1)
    for (i, j), coeff in ising.quadratic.items():
        if coeff * ising.scale <= -args.chain_strength:
            uf.union(i, j)

    groups = {}
    group_of_var = {}
    for var in ising.variables:
        group = uf.find(var)
        if group not in groups:
            groups[group] = []
        groups[group].append(var)
        group_of_var[var] = group

    for group, members in groups.items():
        for var in members:
            if group != var:
                ising.merge(group, var)

    print('Reduced number of variables: {}'.format(len(ising.variables)))
    print('Reduced number of edges: {}'.format(len(ising.quadratic)))
    print()

    variables, linear, quadratic, scale, offset = ising.to_boolean_domain()

    m = Model()

    if args.runtime_limit != None:
        m.setParam('TimeLimit', args.runtime_limit)

    m.setParam('Threads', args.thread_limit)

    if args.cuts != None:
        m.setParam('Cuts', args.cuts)

    gurobi_variables= {}
    for var in variables:
        gurobi_variables[var, var] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name='group_%04d' % var)
    for pair in quadratic:
        gurobi_variables[pair] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name='product_%04d_%04d' % pair)
    m.update()

    for i,j in quadratic:
        m.addConstr(gurobi_variables[(i,j)] >= gurobi_variables[(i,i)] + gurobi_variables[(j,j)] - 1)
        m.addConstr(gurobi_variables[(i,j)] <= gurobi_variables[(i,i)])
        m.addConstr(gurobi_variables[(i,j)] <= gurobi_variables[(j,j)])

    obj = 0.0
    for i, coeff in linear.items():
        obj += coeff * gurobi_variables[i, i]

    for (i, j), coeff in quadratic.items():
        obj += coeff * gurobi_variables[i, j]

    m.setObjective(obj, GRB.MINIMIZE)

    m.update()

    m._cut_count = 0
    m.optimize(cut_counter)

    if args.show_solution:
        print('')
        for group, members in groups.items():
            for var in members:
                v = gurobi_variables[group]
                print('site_{:04d}: {}'.format(var, v.X))

    lower_bound = m.MIPGap*m.ObjVal + m.ObjVal
    scaled_objective = scale * (m.ObjVal + offset)
    scaled_lower_bound = scale * (lower_bound + offset)

    print('')
    print('BQP_DATA, %d, %d, %f, %f, %f, %f, %f, %d, %d' % (len(variables), len(quadratic), scaled_objective, scaled_lower_bound, m.ObjVal, lower_bound, m.Runtime, m._cut_count, m.NodeCount))


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
    parser.add_argument('-cs', '--chain_strength', help='chain strength threshold', type=float, default=1.0)
    parser.add_argument('-ss', '--show-solution', help='print the solution', action='store_true', default=False)
    parser.add_argument('-rtl', '--runtime-limit', help='gurobi runtime limit (sec.)', type=float)
    parser.add_argument('-tl', '--thread-limit', help='gurobi thread limit', type=int, default=1)
    parser.add_argument('-cuts', help='gurobi cuts parameter', type=int)
    return parser


if __name__ == '__main__':
    parser = build_cli_parser()
    main(parser.parse_args())
