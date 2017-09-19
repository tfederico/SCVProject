#/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import json
import maxflow
from ortools.linear_solver import pywraplp
import os

class CCGSolver():
    '''
    Class that encapsulates solution method for a Minimum Weighted Vertex
    Cover Problem on the CCG. The class allows to:

    - Inspect the CCG
    - Apply the Nemhauser Totter reduction (a form of kernelization)
    - Heuristically fix some variables
    - Solve the problem to optimality (conditioned to the fixed variables)
    '''

    def __init__(self, ccgdata):
        '''
        Build a CCGSolver starting from a JSON description of the CCG. The
        internal fields of the class are:
        - offset: basic, fixed value that should be added to the cost returned
          by this solver in order to obtain a solution value for the original
          Markov Field MAP problem
        - nodes: set containing all nodes in the CCG. This includes both
          "main" nodes (corresponding to the original Markov Field variables)
          and auxiliary nodes, added during the construction of the CCG
        - wgts: dictionary with the weight of each node
        - arcs: all the arcs in the CCG
        - aux: list of the auxiliary (non-thorn) nodes in the CCG. This is a 
          subset of the "nodes" field
        - thorns: list of the "thorn" auxiliary nodes in the CCG. This is a
          subset of the "nodes" field
        - included: set of nodes that must be included in the cover. This list
          is populated automatically by the kernelization methods and by the
          optimal solver. It can be modified manually via the "include"  method
          to make ad-hoc assignments
        - excluded: set of node that must be excluded from the cover. See the
          comments for the "included" field
        - model_time: total time spent for building the optimization models
        - sol_time: total time spent for solving the optimization problems
        - solvalue: value of the best solution found. This is available only
          once all node variables have been assigned
        - solbound: best known bound on the cost. NOTE: if ad-hoc assignments
          are made, the bound take them into account. By doing so, it ceases
          to be a valid bound for the original problem
        '''
        self.offset = ccgdata['offset']
        self.wgts = {int(i): ccgdata['wgts'][i] for i in ccgdata['wgts']}
        self.arcs = ccgdata['arcs']
        self.aux = set(ccgdata['aux'])
        self.thorns = set(ccgdata['thorns'])
        # Set of all nodes (provided for convenience)
        self.nodes = set(int(v) for v in self.wgts.keys())

        # Set of necessarily included nodes
        # NOTE: this is also used to store the effect of the kernelizatio and
        # the problem solutions
        self.included = set([])
        # Set of necessarily excluded nodes
        # NOTE: this is also used to store the effect of the kernelizatio and
        # the problem solutions
        self.excluded = set([])

        # Prepare a field to store the elapsed time
        self.model_time = 0
        self.sol_time = 0

        # Data about solution value and bound
        self.solvalue = None
        self.solbound = None

    def include(self, node):
        '''
        Force the inclusion of a node in the cover
        '''
        self.included.add(node)

    def exclude(self, node):
        '''
        Force the exclusion of a node in the cover
        '''
        self.excluded.add(node)

    def solve(self, solver='GLOP', kernelize=True, time_limit=None):
        '''
        Solve the MWVC Problem or apply the Nemhauser-Totter reduction.
        Parameters:
        - solver: the linear programming (or MILP) solver to be used. By
          default this is GLOP (linear programming only), but CLP (linear
          programming) or CBC (Mixed Integer Linear Programming) can also be
          used. In particular CLP is usually faster than GLOP and CBC can also
          attempt to solve the problem to optimality
        - kernelize: if True, the NT reduction is applied to fix as many
          variables as possible. If False, then the solver attempt to find the
          optimal solution. The solver must be CBC in this case
        - time_limit: a time limit for the process
        '''
        # Build an LP solver
        if solver == 'GLOP' and kernelize:
            slv = pywraplp.Solver('CCG Solver',
                                  pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)
        elif solver == 'CLP' and kernelize:
            slv = pywraplp.Solver('CCG Solver',
                                  pywraplp.Solver.CLP_LINEAR_PROGRAMMING)
        elif solver == 'CBC':
            slv = pywraplp.Solver('CCG Solver',
                                  pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
        else:
            s = 'Invalid parameter combination: solver="%s", kernelize="%s"' %\
                (solver, 'True' if kernelize else 'False')
            print(s)
            return False
        # Apply a time limit
        if time_limit is not None:
            slv.set_time_limit(time_limit)
        # Define the variables
        if kernelize:
            x = {i:slv.NumVar(0, 1, 'x%d' % i) for i in self.nodes
                 if i not in self.included and i not in self.excluded}
        else:
            x = {i:slv.IntVar(0, 1, 'x%d' % i) for i in self.nodes
                 if i not in self.included and i not in self.excluded}
        # Define the constraints (taking into account the forced choices)
        for i, j in self.arcs:
            # If any variable is already included, skip the constraint
            if i in self.included or j in self.included:
                pass
            # If i is excluded, than j must be included
            elif i in self.excluded:
                slv.Add(x[j] >= 1)
            # If j is excluded, than i must be included
            elif j in self.excluded:
                slv.Add(x[i] >= 1)
            # Otherwise, post a normal constraint
            else:
                slv.Add(x[i] + x[j] >= 1)
        # Define the problem objective
        obj = slv.Objective()
        for i in self.nodes:
            if i not in self.included and i not in self.excluded:
                obj.SetCoefficient(x[i], self.wgts[i])
        obj.SetMinimization()
        # Keep track of passing time
        self.model_time += slv.WallTime()
        # print('Kernelization model built in %.3f sec' % slv.WallTime())
        # Solve the LP problem
        start = slv.WallTime()
        status = slv.Solve()
        if status == slv.OPTIMAL:
            time = slv.WallTime() - start
            self.sol_time += time
            # print('Kernelization solved in %.3f sec' % time)
        elif status == slv.FEASIBLE:
            print('WARNING: Could not solve the problem to optimality')
        else:
            print('ERROR: Could not solve the problem')
            return False
        # Perform the actual kernelization
        for i in self.nodes:
            if i not in self.included and i not in self.excluded:
                val = x[i].SolutionValue()
                if abs(val - 0) < 1e-5:
                    self.excluded.add(i)
                elif abs(val - 1) < 1e-5:
                    self.included.add(i)
        # Compute the solution bound
        if kernelize:
            self.solbound = obj.Value()
        else:
            self.solbound = obj.BestBound()
        # Obtain the solution value, if available
        if not kernelize or \
           len(self.included) + len(self.excluded) == len(self.nodes):
            self.solvalue = obj.Value()
        return True


if __name__ == '__main__':
    # Load all file names in the benchmars
    target = 'benchmarks'
    files = [f for f in os.listdir(target) if f.endswith('.ccg')]

    # Focus on a specific instance, to provide an example
    # files = [files[0]]

    # Process all benchmark files
    for fname in files:
        fullname = os.path.join(target, fname)
        # Load data
        print('=' * 78)
        print('Loading %s' % fullname)
        with open(fullname) as fp:
            ccgdata = json.load(fp)

        # Define a time limit
        time_limit=10000

        # ====================================================================
        # Example 1: build a solver an apply a kernelization
        # ====================================================================

        print('-' * 78)
        ccgsolver = CCGSolver(ccgdata)

        # Perform one round of kernelization
        print('FIRST ROUND OF KERNELIZATION')
        oldfixed = ccgsolver.included.union(ccgsolver.excluded)
        ccgsolver.solve(solver='GLOP', kernelize=True, time_limit=time_limit)
        # Obtain some stats
        model_time = ccgsolver.model_time
        sol_time = ccgsolver.sol_time
        print('Kernelization model built in %.3f sec' % model_time)
        print('Kernelization model solved in %.3f sec' % sol_time)
        print('Current bound: %f' % ccgsolver.solbound)
        if ccgsolver.solvalue is not None:
            print('Current solution: %f' % ccgsolver.solvalue)
        newfixed = ccgsolver.included.union(ccgsolver.excluded)
        difffixed= newfixed - oldfixed
        remaining = ccgsolver.nodes - newfixed
        print('Fixed %d variables, %d remaining' %
               (len(difffixed), len(remaining)))

        print()
        print('FIXING THE FIRST VARIABLE')
        remaining = list(remaining) # So that we can use indexing
        ccgsolver.include(remaining[0])

        # Perform one round of kernelization
        print()
        print('SECOND ROUND OF KERNELIZATION')
        oldfixed = ccgsolver.included.union(ccgsolver.excluded)
        ccgsolver.solve(solver='GLOP', kernelize=True, time_limit=time_limit)
        # Obtain some stats
        model_time = ccgsolver.model_time
        sol_time = ccgsolver.sol_time
        print('Kernelization model built in %.3f sec' % model_time)
        print('Kernelization model solved in %.3f sec' % sol_time)
        print('Bound, given the heuristic assignemt: %f' % ccgsolver.solbound)
        if ccgsolver.solvalue is not None:
            print('Current solution: %f' % ccgsolver.solvalue)
        newfixed = ccgsolver.included.union(ccgsolver.excluded)
        difffixed= newfixed - oldfixed
        remaining = ccgsolver.nodes - newfixed
        print('Fixed %d variables, %d remaining' %
               (len(difffixed), len(remaining)))


        # ====================================================================
        # Example 2: solve to optimality (not always viable)
        # NOTE: this requires to have CBC installed
        # ====================================================================

        print('-' * 78)
        ccgsolver = CCGSolver(ccgdata)

        print('FIRST ROUND OF KERNELIZATION')
        oldfixed = ccgsolver.included.union(ccgsolver.excluded)
        ccgsolver.solve(solver='CLP', kernelize=True, time_limit=time_limit)
        # Obtain some stats
        model_time = ccgsolver.model_time
        sol_time = ccgsolver.sol_time
        print('Kernelization model built in %.3f sec' % model_time)
        print('Kernelization model solved in %.3f sec' % sol_time)
        print('Current bound: %f' % ccgsolver.solbound)
        if ccgsolver.solvalue is not None:
            print('Current solution: %f' % ccgsolver.solvalue)
        newfixed = ccgsolver.included.union(ccgsolver.excluded)
        difffixed= newfixed - oldfixed
        remaining = ccgsolver.nodes - newfixed
        print('Fixed %d variables, %d remaining' %
               (len(difffixed), len(remaining)))

        print('SOLVING TO OPTIMALITY')
        ccgsolver.solve(solver='CBC', kernelize=False, time_limit=time_limit)
        # Obtain some stats
        model_time = ccgsolver.model_time
        sol_time = ccgsolver.sol_time
        print('MWVCP model built in %.3f sec' % model_time)
        print('MWVCP model solved in %.3f sec' % sol_time)
        print('Current bound: %f' % ccgsolver.solbound)
        if ccgsolver.solvalue is not None:
            print('Current solution value: %f' % ccgsolver.solvalue)

