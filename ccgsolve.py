# -*- coding: utf-8 -*-

from __future__ import print_function

import json
#import maxflow
from ortools.linear_solver import pywraplp
import os
import sys
import argparse
import operator
import math

def mean(list):
    return float(sum(list))/ len(list)

def std(list):
    avg = mean(list)
    total = 0.0
    for element in list:
        total = total + (element - avg)**2
    total = float(total)
    return math.sqrt(total/len(list))



def nodeToInclude(remaining, euristica):

    remaining = list(remaining) # So that we can use indexing
    toInclude = None

    # peso dei nodi di tutti i nodi che rimangono da fissare
    remaininWeights = {}
    for node in remaining:
        weight = ccgsolver.wgts[node]
        remaininWeights[node] = weight

    if euristica == 'maxPeso':
        # nodo che rimane da fissare con peso massimo
        toInclude = max(remaininWeights.iteritems(), key=operator.itemgetter(1))[0]

    elif(euristica == 'nodeMinPeso'):
        toInclude = min(remaininWeights.iteritems(), key=operator.itemgetter(1))[0]

    elif euristica == 'thorn' or euristica == 'thornMaxPeso' or euristica == 'thornMinPeso':
        remainingThorns = [node for node in remaining if node in ccgsolver.thorns]
        # peso dei nodi dei thorn che rimangono da fissare
        remaininWeightsThorns = {}
        for node in remainingThorns:
            weight = ccgsolver.wgts[node]
            remaininWeightsThorns[node] = weight
        if euristica == 'thorn':
            toInclude = remainingThorns[0]
        elif euristica == 'thornMaxPeso':
            # thorn che rimane da fissare con peso max
            toInclude = max(remaininWeightsThorns.iteritems(), key=operator.itemgetter(1))[0]
        else:
            toInclude = min(remaininWeightsThorns.iteritems(), key=operator.itemgetter(1))[0]
    elif euristica == 'aux' or euristica == 'auxMaxPeso' or euristica == 'auxMinPeso':
        # aux che rimangono da fissare
        remainingAux = [node for node in remaining if node in ccgsolver.aux]
        # peso dei nodi degli aux che rimangono da fissare
        remaininWeightsAux = {}

        for node in remainingAux:
            weight = ccgsolver.wgts[node]
            remaininWeightsAux[node] = weight
        if euristica == 'aux' :
            toInclude = remainingAux[0]
        elif euristica == 'auxMaxPeso':
            # aux che rimane da fissare con peso max
            toInclude = max(remaininWeightsAux.iteritems(), key=operator.itemgetter(1))[0]
        else:
            toInclude = min(remaininWeightsAux.iteritems(), key=operator.itemgetter(1))[0]
    elif euristica == 'noAuxNoThorn' or euristica == 'noAuxNoThornMaxPeso' or euristica == 'noAuxNoThornMinPeso':
        #thorn che rimangono da fissare
        remainingThorns = [node for node in remaining if node in ccgsolver.thorns]
        # aux che rimangono da fissare
        remainingAux = [node for node in remaining if node in ccgsolver.aux]

        # nè thorn nè aux che rimangono da fissare
        remainingNone = [node for node in remaining
                        if node not in remainingAux and node not in remainingThorns]
	remaininWeightsNone = {}
	for node in remainingNone:
            weight = ccgsolver.wgts[node]
            remaininWeightsNone[node] = weight

        if euristica == 'noAuxNoThorn':
            toInclude = remainingNone[0]
        elif euristica == 'noAuxNoThornMaxPeso':
            # nè thorn nè aux che rimangono da fissare che rimane da fissare con peso max e min
            toInclude = max(remaininWeightsNone.iteritems(), key=operator.itemgetter(1))[0]
        else:
            toInclude = min(remaininWeightsNone.iteritems(), key=operator.itemgetter(1))[0]

    elif nomeEuristica == 'maxArchiDaFissare' or nomeEuristica == 'minArchiDaFissare':
        # Conta il numero di archi per ogni nodo
        nodeCounts = {}
        for i, j in ccgsolver.arcs:
            if i in remaining: # se i è da fissare
                # se i è connesso ad un nodo da fissare aumentiamo il count
                if i in nodeCounts and j in remaining:
                    nodeCounts[i] = nodeCounts[i] + 1
                elif j in remaining:
                    nodeCounts[i] = 1
            # stessa cosa per j
            if j in remaining:
                if j in nodeCounts and i in remaining:
                    nodeCounts[j] += 1
                elif i in remaining:
                    nodeCounts[j] = 1

        if(euristica == 'maxArchiDaFissare'):
            # nodi da fissare con max e min numero di archi che lo connettono ad un nodo da fissare
            toInclude = max(nodeCounts.iteritems(), key=operator.itemgetter(1))[0]
        else:
            toInclude = min(remaininWeights.iteritems(), key=operator.itemgetter(1))[0]

    if toInclude == None:
        return remaining[0]
    return toInclude

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
        !- nodes: set containing all nodes in the CCG. This includes both
          "main" nodes (corresponding to the original Markov Field variables)
          and auxiliary nodes, added during the construction of the CCG
        !- wgts: dictionary with the weight of each node
        - arcs: all the arcs in the CCG
        !- aux: list of the auxiliary (non-thorn) nodes in the CCG. This is a
          subset of the "nodes" field
        !- thorns: list of the "thorn" auxiliary nodes in the CCG. This is a
          subset of the "nodes" field
        !- included: set of nodes that must be included in the cover. This list
          is populated automatically by the kernelization methods and by the
          optimal solver. It can be modified manually via the "include"  method
          to make ad-hoc assignments
        !- excluded: set of node that must be excluded from the cover. See the
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


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-eur', '--eur', help="Euristic chosen. Options: thorn, " \
    + "aux, noAuxNoThorn, maxPeso, minPeso, thornMaxPeso, thornMinPeso, " \
    + "auxMaxPeso, auxMinPeso, noAuxNoThornMaxPeso, noAuxNoThornMinPeso, maxArchiDaFissare," \
    + "minArchiDaFissare", required=True)


    args = parser.parse_args()

    nomeEuristica = args.eur

    passiTOT=[]
    mediaFissatiTOT=[]
    devStdNodiFissatiTOT=[]
    maxFissatiTOT=[]
    minFissatiTOT=[]
    noFissatiTOT=[]
    soluzioneTOT=[]



    # Load all file names in the benchmars
    target = 'benchmarks'
    files = [f for f in os.listdir(target) if f.endswith('.ccg')]

    # Focus on a specific instance, to provide an example
    #files = [files[0]]

    # Process all benchmark files
    noImprovements = 0
    minImprovements = sys.maxint
    maxImprovements = 0

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
        #ccgsolver.solve(solver='GLOP', kernelize=True, time_limit=time_limit)
        ccgsolver.solve(solver='CBC', kernelize=True)
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

        ################################################################ POCI

        varfixed = []
        while len(list(remaining)) != 0:

            node = nodeToInclude(remaining, nomeEuristica)

            ccgsolver.include(node)

            # Perform one round of kernelization
            print()
            print('SECOND ROUND OF KERNELIZATION')
            oldfixed = ccgsolver.included.union(ccgsolver.excluded)
            #ccgsolver.solve(solver='GLOP', kernelize=True, time_limit=time_limit)
            ccgsolver.solve(solver='CBC', kernelize=True)

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
            if len(difffixed) == 0:
                noImprovements = noImprovements + 1

            if len(difffixed) > maxImprovements:
                maxImprovements = len(difffixed)

            if len(difffixed) < minImprovements:
                minImprovements = len(difffixed)

            varfixed.append(len(difffixed))


        '''
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
            print('Current solution value: %f' % ccgsolver.solvalue)'''

        # nome euristica
        # numero passi
        # valore soluzione
        # numero medio var fissate
        # dev std numero medio var fissate
        # num max var fissate
        # num min var fissate
        # numero di iterazioni che non fissano




        avg_varfixed = mean(varfixed)
        #print('Total nodes', len(ccgsolver.nodes))

        #print('thorns', len(ccgsolver.thorns))
        #print('aux', len(ccgsolver.aux))

        std_varfixed = std(varfixed)

        file = open("results/" + nomeEuristica + '_' + fname + '.txt', 'w')
        file.write('Euristica: ' + nomeEuristica + '\n')
        file.write('Problema: ' + fname + '\n')
        file.write('Passi: ' + str(len(varfixed)) + '\n')
        file.write('Media nodi fissati: ' + str(avg_varfixed) + '\n')
        file.write('Dev std dei nodi fissati: ' + str(std_varfixed) + '\n')
        file.write('Max nodi fissati: ' + str(maxImprovements) + '\n')
        file.write('Min nodi fissati: ' + str(minImprovements) + '\n')
        file.write('Iterazioni che non fissano nodi: ' + str(noImprovements) + '\n')
        file.write('Soluzione: ' + str(ccgsolver.solvalue) + '\n')
        file.close()

        passiTOT.append(len(varfixed))
        mediaFissatiTOT.append(avg_varfixed)
        devStdNodiFissatiTOT.append(std_varfixed)
        maxFissatiTOT.append(maxImprovements)
        minFissatiTOT.append(minImprovements)
        noFissatiTOT.append(noImprovements)
        soluzioneTOT.append(ccgsolver.solvalue)


file = open("results/" + nomeEuristica + '.txt', 'w')
file.write('Media passi: ' + str(mean(passiTOT)) + '\n')
file.write('Media fissati: ' + str(mean(mediaFissatiTOT)) + '\n')
file.write('Media max nodi fissati: ' + str(mean(maxFissatiTOT)) + '\n')
file.write('Media min nodi fissati: ' + str(mean(minFissatiTOT)) + '\n')
file.write('Media nodi no fissati: ' + str(mean(noFissatiTOT)) + '\n')
file.write('Media soluzioni tot: ' + str(mean(soluzioneTOT)) + '\n')

file.write('Dev std passi: ' + str(std(passiTOT)) + '\n')
file.write('Dev std fissati: ' + str(std(mediaFissatiTOT)) + '\n')
file.write('Dev std max nodi fissati: ' + str(std(maxFissatiTOT)) + '\n')
file.write('Dev std min nodi fissati: ' + str(std(minFissatiTOT)) + '\n')
file.write('Dev std nodi no fissati: ' + str(std(noFissatiTOT)) + '\n')
file.write('Dev std soluzioni tot: ' + str(std(soluzioneTOT)) + '\n')

file.write('Max passi: ' + str(max(passiTOT)) + '\n')
file.write('Max fissati: ' + str(max(mediaFissatiTOT)) + '\n')
file.write('Max max nodi fissati: ' + str(max(maxFissatiTOT)) + '\n')
file.write('Max min nodi fissati: ' + str(max(minFissatiTOT)) + '\n')
file.write('Max nodi no fissati: ' + str(max(noFissatiTOT)) + '\n')
file.write('Max soluzioni tot: ' + str(max(soluzioneTOT)) + '\n')

file.write('Min passi: ' + str(min(passiTOT)) + '\n')
file.write('Min fissati: ' + str(min(mediaFissatiTOT)) + '\n')
file.write('Min max nodi fissati: ' + str(min(maxFissatiTOT)) + '\n')
file.write('Min min nodi fissati: ' + str(min(minFissatiTOT)) + '\n')
file.write('Min nodi no fissati: ' + str(min(noFissatiTOT)) + '\n')
file.write('Min soluzioni tot: ' + str(min(soluzioneTOT)) + '\n')

file.close()
