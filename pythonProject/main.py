from __future__ import print_function
from ortools.linear_solver import pywraplp
import networkx as nx
import matplotlib.pyplot as plt


# from operator import itemgetter


def main():
    WP = 0  # max number of waypoints allowed per demand
    pos = {}
    M = 100  # a constant large enough
    nodes = []
    links = []
    demands = []
    segments = []

    # Create the mip solver with the CBC backend.
    solver = pywraplp.Solver.CreateSolver('simple_mip_program', 'CBC')
    infinity = solver.infinity()

    def init():
        nonlocal M, segments
        M = max(sum([d[2] for d in demands]), 2 * len(links), 100)  # a constant large enough
        # pairs of nodes p,q where at least one of p or q is a source/terminal
        segments = [(p, q) for p in nodes for q in nodes if p != q]
        for i, (s, t, d) in enumerate(demands):
            demands[i] = (s, t, d, i)

    def example1():
        # Example 1
        nonlocal nodes, links, demands
        nodes = [1, 2, 3, 4]
        links = [(1, 2, 1), (1, 3, 2), (2, 4, 0.5), (3, 4, 2.5), (2, 3, 0.5)]  # (u,v,capacity)
        demands = [(1, 4, 1, 0), (1, 4, 1, 1)]  # (s,t,d,i)   Todo: check s != t

    def example2():
        nonlocal nodes, links, demands
        nodes = [1, 3, 2, 4, 5, 6, 141, 142, 231, 232, 361, 362, 461, 462]
        # pos = {1: (20, 20), 2: (30, 20), 3: (25, 10), 4: (25, 30), 5: (23, 20), 6: (27, 20)}
        # links = [(1, 3, 1), (1, 4, 2), (2, 3, 2), (2, 4, 1), (3, 5, 1), (3, 6, 2), (4, 5, 1), (4, 6, 2)]  # (u,v,capacity)
        links = [(1, 3, 1),
                 (1, 141, 1), (1, 142, 1), (141, 4, 1), (142, 4, 1),
                 (2, 231, 1), (2, 232, 1), (231, 3, 1), (232, 3, 1),
                 (2, 4, 1), (3, 5, 1),
                 (3, 361, 1), (3, 362, 1), (361, 6, 1), (362, 6, 1),
                 (4, 5, 1),
                 (4, 461, 1), (4, 462, 1), (461, 6, 1), (462, 6, 1)]  # (u,v,capacity)
        demands = [(1, 5, 1, 0), (1, 6, 2, 1), (2, 5, 1, 2), (2, 6, 2, 3)]  # (s,t,d,i)   Todo: check s != t

    def ex_uniform_capacities():
        nonlocal nodes, links, demands
        nodes = [1, 2, 3, 4]
        links = [(1, 2, 1), (1, 3, 1), (1, 4, 1), (2, 4, 1), (3, 4, 1)]
        demands = [(1, 4, 2), (1, 4, 1.5), (2, 4, 0.5)]

    def example_grid(m, n):
        for i in range(m):
            nodes.append(i)
            for j in range(m):
                if j == i:
                    continue
                links.append((i, j, 1))

    def draw():
        G = nx.DiGraph(directed=True)
        # G.add_nodes_from(nodes)
        for l in links:
            G.add_edge(l[0], l[1], capacity=l[2])
        # print("Nodes of graph: ")
        # print(G.nodes())
        # print("Edges of graph: ")
        # print(G.edges())

        options = {
            # 'node_color': 'blue',
            'node_size': 300,
            'width': 1,
            'arrowstyle': '->',
            'arrowsize': 20,
        }
        edge_labels = dict([((u, v,), 'c=' + str(data['capacity']))
                            for u, v, data in G.edges(data=True)])
        nonlocal pos
        if not pos:
            pos = nx.circular_layout(G)
        nx.draw_networkx(G, pos, arrows=True, with_labels=True, **options)
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels);
        plt.show()

    def x(v, l):
        l_str = str(l[0]) + ',' + str(l[1])
        return 'x^' + str(v) + '_{' + l_str + '}'

    def f(p, q, l):
        l_str = str(l[0]) + ',' + str(l[1])
        pq = str(p) + ',' + str(q)
        return 'f^{' + pq + '}_{' + l_str + '}'

    def f_v(t, v):
        return 'f^' + str(t) + '_' + str(v)

    def seg(s, t, i, p, q):
        return 'S^{' + str((s, t)) + '_' + str(i) + '}_{' + str((p, q)) + '}'

    def dist(v, t):
        return 'd_{' + str(v) + ',' + str(t) + '}'

    def w(l):
        return 'w_{' + str(l[0]) + ',' + str(l[1]) + '}'

    def Z(l):
        return 'z_{' + str(l[0]) + ',' + str(l[1]) + '}'

    def I(l):
        return 'I_{' + str(l[0]) + ',' + str(l[1]) + '}'

    def variables():
        init()
        # objective variables
        solver.NumVar(0.0, M, 'L')

        for p, q in segments:
            for l in links:
                solver.NumVar(0.0, M, f(p, q, l))

        for s, t, d, i in demands:
            for p, q in segments:
                solver.BoolVar(seg(s, t, i, p, q))

        for v in nodes:
            for l in links:
                solver.BoolVar(x(v, l))

        for p, q in segments:
            solver.NumVar(0.0, M, f_v(p, q))

        for t in nodes:
            for v in nodes:
                solver.NumVar(0.0, infinity, dist(v, t))

        for l in links:
            solver.NumVar(0.0, 1.0, w(l))

        print('Number of variables =', solver.NumVariables())

    # [END variables]

    #  flow constraints, per segment
    def flows():
        for p, q in segments:
            for v in nodes:
                f_pq_out = []
                f_pq_in = []
                for l in links:
                    if l[0] == v:
                        f_pq_out.append(solver.LookupVariable(f(p, q, l)))
                    elif l[1] == v:
                        f_pq_in.append(solver.LookupVariable(f(p, q, l)))
                s_sti_pq = []
                for s, t, d, i in demands:
                    s_sti_pq.append(solver.LookupVariable(seg(s, t, i, p, q)) * d)
                d_pq = solver.Sum(s_sti_pq);  # total flow on (p,q)

                if v == p:
                    solver.Add(solver.Sum(f_pq_in) - solver.Sum(f_pq_out) == -d_pq)  # flow source
                elif v == q:
                    solver.Add(solver.Sum(f_pq_in) - solver.Sum(f_pq_out) == d_pq)  # flow sink
                else:
                    solver.Add(solver.Sum(f_pq_in) == solver.Sum(f_pq_out))  # flow conservation

    # segment constraints, per demand
    def segments_paths():
        for s, t, d, i in demands:
            for v in nodes:
                s_stiv_out = []
                s_stiv_in = []
                for p, q in segments:
                    if v == p:
                        s_stiv_out.append(solver.LookupVariable(seg(s, t, i, p, q)))
                    elif v == q:
                        s_stiv_in.append(solver.LookupVariable(seg(s, t, i, p, q)))

                if v == s:
                    solver.Add(solver.Sum(s_stiv_in) - solver.Sum(s_stiv_out) == -1)
                elif v == t:
                    solver.Add(solver.Sum(s_stiv_in) - solver.Sum(s_stiv_out) == 1)
                else:
                    solver.Add(solver.Sum(s_stiv_out) == solver.Sum(s_stiv_in))

        # number of segments, per demand
        for s, t, d, i in demands:
            s_sti_pq = []
            for p, q in segments:
                s_sti_pq.append(solver.LookupVariable(seg(s, t, i, p, q)))
            solver.Add(solver.Sum(s_sti_pq) <= WP + 1)

    # capacity constraints
    def capacity():
        L = solver.LookupVariable('L')
        for l in links:
            f_pql = []
            for p, q in segments:
                f_pql.append(solver.LookupVariable(f(p, q, l)))
            solver.Add(solver.Sum(f_pql) <= L * l[2])

    # shortest path tree constraints
    def SP_tree():
        for p, q in segments:
            for l in links:
                f_pql = solver.LookupVariable(f(p, q, l))
                x_ql = solver.LookupVariable(x(q, l))
                solver.Add(f_pql <= M * x_ql)

    # even split constraints
    def even_split():
        for v, t in segments:  # all (v,t) pairs, v!=t
            v_out = [l for l in links if l[0] == v]  # out-links from v
            for l in v_out:
                f_ptl = []
                for p in nodes:
                    if p == t: continue
                    f_ptl.append(solver.LookupVariable(f(p, t, l)))

                f_tv = solver.LookupVariable(f_v(t, v))
                x_tl = solver.LookupVariable(x(t, l))
                solver.Add(solver.Sum(f_ptl) <= f_tv)
                solver.Add(f_tv - solver.Sum(f_ptl) <= M * (1 - x_tl))
                # print(x_tl, f_ptl, sep='<=')

    # weight-setting constraints
    def weights():
        # for v in nodes:
        #     solver.Add(solver.LookupVariable(dist(v, v)) == 0)

        for l in links:
            u = l[0]
            v = l[1]
            for t in nodes:
                if t == u: continue
                d_vt = solver.LookupVariable(dist(v, t))
                d_ut = solver.LookupVariable(dist(u, t))
                w_l = solver.LookupVariable(w(l))
                x_tl = solver.LookupVariable(x(t, l))
                solver.Add(d_ut <= d_vt + w_l)
                solver.Add(d_vt - d_ut + w_l <= M * (1 - x_tl))
                # print(1, '-', x_tl, '<= M * (', d_vt, '-', d_ut, '+', w_l, ')')
                solver.Add(1 - x_tl <= M * (d_vt - d_ut + w_l))

    # not part of the MIP, for possible future experiments
    # produces link weights that maximizes the worst congestion
    # don't forget to comment the original objective before using this
    def maximize_congestion():
        # not part of the main MIP
        nonlocal WP
        WP = 0
        for l in links:
            solver.NumVar(0.0, M, Z(l))

        for l in links:
            solver.BoolVar(I(l))

        Z_sum = []
        I_sum = []
        for l in links:
            I_l = solver.LookupVariable(I(l))
            I_sum.append(I_l)
            Z_l = solver.LookupVariable(Z(l))
            Z_sum.append(Z_l)
            solver.Add(Z_l <= I_l * M)
            f_pql = []
            for p, q in segments:
                f_pql.append(solver.LookupVariable(f(p, q, l)))
            solver.Add(Z_l <= solver.Sum(f_pql) / l[2])  # upper bound by congestion

        solver.Add(solver.Sum(I_sum) == 1)
        solver.Maximize(solver.Sum(Z_sum))

    ex_uniform_capacities()

    WP = 0
    variables()
    # draw()
    # maximize_congestion()
    flows()
    segments_paths()
    capacity()
    SP_tree()
    even_split()
    # weights()

    print('Number of constraints =', solver.NumConstraints())

    # [END constraints]

    # [START objective]
    solver.Minimize(solver.LookupVariable('L'))

    # [END objective]

    # [START solve]
    status = solver.Solve()
    # [END solve]
    print("status=", status)
    # [START print_solution]
    if status == pywraplp.Solver.OPTIMAL:
        print('Solution:')
        print('Objective value =', solver.Objective().Value())
        # print('variables =', solver.variables())
        variable_list = solver.variables()
        for variable in variable_list:
            if variable.solution_value() > 0:
                print(('%s = %f' % (variable.name(), variable.solution_value())))

    elif status == pywraplp.Solver.INFEASIBLE:
        print('Not feasible.')
    else:
        print('The problem does not have an optimal solution.')

    # [END print_solution]

    # [START advanced]
    print('\nAdvanced usage:')
    print('Problem solved in %f milliseconds' % solver.wall_time())
    print('Problem solved in %d iterations' % solver.iterations())
    print('Problem solved in %d branch-and-bound nodes' % solver.nodes())

    # [END advanced]


if __name__ == '__main__':
    main()
# [END program]
