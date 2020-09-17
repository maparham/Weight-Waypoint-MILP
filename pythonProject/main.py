from __future__ import print_function
from ortools.linear_solver import pywraplp


# from operator import itemgetter


def main():
    demands = [(1, 4, 1, 0), (1, 4, 2, 1)]  # s,t,d,i   Todo: check s != t
    nodes = [1, 2, 3, 4]
    links = [(1, 2, 1), (1, 3, 1), (2, 4, 1), (3, 4, 1)]  # u,v,capacity

    # st_pairs = list(map(itemgetter(0, 1), demands))
    segments = []  # pairs of nodes p,q where at least one of p or q is a source/terminal
    for s, t, d, i in demands:
        segments.append((s, t))
        for v in nodes:
            if v == s or v == t:
                continue  # to ensure p!=q for segment (p,q)
            segments += [(s, v), (v, t)]
    segments = set(segments)  # remove duplicates

    def x(v, l):
        l_str = str(l[0]) + str(l[1])
        return 'x^' + str(v) + '_' + l_str

    def f(p, q, l):
        l_str = str(l[0]) + ',' + str(l[1])
        pq = str(p) + ',' + str(q)
        return 'f^{' + pq + '}_{' + l_str + '}'

    # [START solver]
    # Create the mip solver with the CBC backend.
    solver = pywraplp.Solver.CreateSolver('simple_mip_program', 'CBC')
    # [END solver]

    # [Add variables]
    infinity = solver.infinity()

    # segments = [(p, q) for p in nodes for q in nodes if p != q]

    def w(s, t, i, v):
        return 'w^{' + str((s, t)) + '_' + str(i) + '}_' + str(v)

    def seg(s, t, i, p, q):
        return 'S^{' + str((s, t)) + '_' + str(i) + '}_{' + str((p, q)) + '}'

    for l in links:
        for v in nodes:
            solver.IntVar(0, 1, x(v, l))
        for p, q in segments:
            solver.NumVar(0.0, infinity, f(p, q, l))

    for s, t, d, i in demands:
        for v in nodes:
            solver.IntVar(0, 1, w(s, t, i, v))

    for s, t, d, i in demands:
        for p, q in segments:
            solver.IntVar(0, 1, seg(s, t, i, p, q))

    print('Number of variables =', solver.NumVariables())
    # [END variables]

    # [Add constraints]

    # for s, t in st_pairs:
    #     for v in nodes:
    #         if v == s: continue
    #         f_sv = []
    #         f_vt = []
    #         for l in links:
    #             if l[0] == v:  # outgoing link
    #                 f_vtl = solver.LookupVariable(f(v, t, l))
    #                 if f_vtl is None: assert 'wrong f_vtl!'
    #                 f_vt.append(f_vtl)
    #             if l[1] == v:  # incoming link
    #                 f_svl = solver.LookupVariable(f(s, v, l))
    #                 if f_svl is None: assert 'wrong f_svl!'
    #                 f_sv.append(f_svl)
    #             else:
    #                 continue
    #
    #         if len(f_sv) > 0:  # otherwise v cannot be a waypoint
    #             w_st = []
    #             for p, q, d, i in demands:
    #                 if p != s or q != t: continue
    #                 w_stiv = solver.LookupVariable(w(s, t, i, v))
    #                 w_st.append(w_stiv * d)
    #             solver.Add(solver.Sum(f_sv) == solver.Sum(w_st))
    #             if v != t and len(f_vt) > 0:  # if there v has outgoing links and not a terminal
    #                 solver.Add(solver.Sum(f_vt) == solver.Sum(f_sv))

    # segment constraints
    for p, q in segments:
        f_pql = []
        for l in links:
            if l[0] == p:
                f_pql.append(solver.LookupVariable(f(p, q, l)))
        w_stip = []
        w_stiq = []
        for s, t, d, i in demands:
            if p == s:
                w_stiq.append((solver.LookupVariable(w(s, t, i, q)) * d))
            elif q == t:
                w_stip.append((solver.LookupVariable(w(s, t, i, p)) * d))

        if len(w_stip) > 0 or len(w_stiq) > 0:
            solver.Add(solver.Sum(f_pql) == solver.Sum(w_stip) + solver.Sum(w_stiq))

    # flow conservation constraints
    for p, q in segments:
        for v in nodes:
            if v == p or v == q:
                continue
            f_pq_in = []
            f_pq_out = []
            for l in links:
                f_pql = solver.LookupVariable(f(p, q, l))
                # print(f_pql)
                # print(f(p, q, l))
                if l[0] == v:
                    f_pq_out.append(f_pql)
                if l[1] == v:
                    f_pq_in.append(f_pql)

            if len(f_pq_out) > 0 or len(f_pq_in) > 0:
                # print(solver.Sum(f_pq_out), '==', solver.Sum(f_pq_in))
                solver.Add(solver.Sum(f_pq_out) + 0 == solver.Sum(f_pq_in) + 0)

    # force single WP
    for s, t, d, i in demands:
        w_sti = []
        for v in nodes:
            if v != s:
                w_sti.append(solver.LookupVariable(w(s, t, i, v)))
        solver.Add(solver.Sum(w_sti) == 1)

    # # waypoint constraints
    # # solver.Add(w14_1_4 == 1)
    # # solver.Add(w14_2_4 == 1)
    # solver.Add(w14_1_2 + w14_1_3 + w14_1_4 == 1)
    # solver.Add(w14_2_2 + w14_2_3 + w14_2_4 == 1)
    #
    # # capacity constraints
    # solver.Add(f14_12 + f12_12 <= 1)
    # solver.Add(f14_24 + f24_24 <= 1)
    # solver.Add(f14_13 + f13_13 <= 2)
    # solver.Add(f14_34 + f34_34 <= 2)

    print('Number of constraints =', solver.NumConstraints())

    # [END constraints]

    # [START objective]
    # Maximize x + 10 * y.
    # solver.Maximize(1)
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
            print(('%s = %f' % (variable.name(), variable.solution_value())))

    elif status == pywraplp.Solver.INFEASIBLE:
        print('The problem is not feasible.')
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
