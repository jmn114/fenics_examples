#!/usr/bin/python2

__author__ = 'Maksim Vlasov'

"""
An example of Navier-Stoke coupled solver
with no slip and slip boundary conditions.
"""

import dolfin


# region Set parameters.

# A length of a rectangular prism.
height = 3.0

# Kinematic viscosity.
# nu = 1E-1
# nu = 1E-2
# nu = 1E-3
nu = 1E-5
# nu = 1E-7

# Force.
# f_exp = dolfin.Constant((0.0, 0.0, 0.0))
# f_exp = dolfin.Constant((0.0, 0.0, -1.0))
f_exp = dolfin.Constant((0.0, 0.0, 1.0))
# f_exp = dolfin.Constant((0.0, 0.0, -10.0))
# f_exp = dolfin.Constant((0.0, 0.0, 10.0))

# Start time.
t = 0.0

# End time.
#t_end = 1000.0
t_end = 1.0

# Time step.
# dt = 1E-1
dt = 1E-2

# Max iter.
max_iter = 10000

# Inlet velocity.
u_inlet = 1.0
#u_inlet = -1.0

# Outlet pressure,
p_outlet = 0.0
# p_outlet = 10.0

# Use slip bc or no slip bc.
#is_slip_bcs = True
is_slip_bcs = False

# Plot solution or not.
is_plot_solution = True
# is_plot_solution = False

# Plot solution interactive or not.
# is_plot_solution_interactive = True
is_plot_solution_interactive = False

# Save solution or not.
is_save_solution = True
# is_save_solution = False

# Define a mesh.
# mesh = dolfin.BoxMesh(0.0, 0.0, 0.0, 1.0, 1.0, height, 2, 2, 6)
#mesh = dolfin.BoxMesh(0.0, 0.0, 0.0, 1.0, 1.0, height, 4, 4, 10)
mesh = dolfin.BoxMesh(dolfin.Point(0.0, 0.0, 0.0), dolfin.Point(1., 1., height), 4, 4, 10)
# mesh = dolfin.BoxMesh(0.0, 0.0, 0.0, 1.0, 1.0, height, 8, 8, 20)

# Additional parameters.
max_inner_iter = 5
# max_inner_iter = 2

tol = 1E-2
# tol = 1E-6

# Define linear algebra solver settings.
#method = "bicgstab"
#preconditioner = "bjacobi"

method = "mumps"
preconditioner = "none"

krylov_solver_monitor_convergence = False
krylov_solver_relative_tolerance = 1E-16
krylov_solver_absolute_tolerance = 1E-15
# endregion

# region Define sub domains.

# region Define sub domain (walls).
if not is_slip_bcs:
    class Walls(dolfin.SubDomain):
        # noinspection PyUnusedLocal
        def inside(self, x, on_boundary):
            return on_boundary
else:
    class Wall1(dolfin.SubDomain):
        # noinspection PyUnusedLocal
        def inside(self, x, on_boundary):
            return on_boundary and \
                   (dolfin.near(x[0], 0.0) or dolfin.near(x[0], 1.0))


    class Wall2(dolfin.SubDomain):
        # noinspection PyUnusedLocal
        def inside(self, x, on_boundary):
            return on_boundary and \
                   (dolfin.near(x[1], 0.0) or dolfin.near(x[1], 1.0))


# endregion


# region Define sub domain (inlet).
class Inlet(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        # noinspection PyUnresolvedReferences
        return on_boundary and dolfin.near(x[2], 0.0)


# endregion


# region Define sub domain (outlet).
class Outlet(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        # noinspection PyUnresolvedReferences
        return on_boundary and dolfin.near(x[2], height)


# endregion


# Define mesh boundaries.
boundaries = dolfin.MeshFunction("size_t",
                                 mesh,
                                 mesh.topology().dim() - 1)

boundaries.set_all(0)

if not is_slip_bcs:
    walls = Walls()
else:
    wall_1 = Wall1()
    wall_2 = Wall2()

inlet = Inlet()
outlet = Outlet()

if not is_slip_bcs:
    walls.mark(boundaries, 3)
else:
    wall_1.mark(boundaries, 3)
    wall_2.mark(boundaries, 4)

inlet.mark(boundaries, 1)
outlet.mark(boundaries, 2)
# endregion

# region Plot mesh, boundaries.
# dolfin.plot(mesh)
# dolfin.plot(boundaries)
# dolfin.interactive()
# endregion

# region Define function spaces.
V = dolfin.VectorFunctionSpace(mesh, "CG", 2)
Q = dolfin.FunctionSpace(mesh, "CG", 1)
VQ = V * Q
# endregion

# region Define test and trial functions.
up = dolfin.TrialFunction(VQ)
u, p = dolfin.split(up)

v_u, v_p = dolfin.TestFunctions(VQ)
# endregion

# region Define functions.
up_1 = dolfin.Function(VQ)
up_s = dolfin.Function(VQ)

u_1, p_1 = dolfin.split(up_1)

up_expr = dolfin.Constant((0.0, 0.0, 0.0, 0.0))

up_0 = dolfin.interpolate(up_expr, VQ)
u_0, p_0 = dolfin.split(up_0)
# endregion

# region Define parameters.
nu_exp = dolfin.Constant(nu)
dt_exp = dolfin.Constant(dt)
n = dolfin.FacetNormal(mesh)
p_outlet_exp = dolfin.Constant(p_outlet)
# endregion

# region Define boundary conditions.
bcs_up = None

# region Inlet bc.
bc_u_0 = dolfin.DirichletBC(VQ.sub(0),
                            dolfin.Constant((0.0, 0.0, u_inlet)),
                            boundaries,
                            1)
# endregion

# region No-slip bc.
if not is_slip_bcs:
    bc_u_1 = dolfin.DirichletBC(VQ.sub(0),
                                dolfin.Constant((0.0, 0.0, 0.0)),
                                boundaries,
                                3)
else:
    bc_u_1 = dolfin.DirichletBC(VQ.sub(0).sub(0),
                                dolfin.Constant(0.0),
                                boundaries,
                                3)

    bc_u_2 = dolfin.DirichletBC(VQ.sub(0).sub(1),
                                dolfin.Constant(0.0),
                                boundaries,
                                4)
# endregion

# endregion

# region Collect boundary conditions.
if not is_slip_bcs:
    bcs_up = [bc_u_0, bc_u_1]
else:
    bcs_up = [bc_u_0, bc_u_1, bc_u_2]
# endregion


dolfin.ds = dolfin.Measure("ds")[boundaries]


# region Define equations.
def epsilon(u):
    """Return the symmetric gradient."""

    return 0.5 * (dolfin.grad(u) + dolfin.grad(u).T)


F = (1 / dt_exp) * dolfin.inner(u - u_1, v_u) * dolfin.dx \
    + dolfin.inner(dolfin.dot(dolfin.grad(u), u_0), v_u) * dolfin.dx \
    + nu_exp * 2.0 * dolfin.inner(epsilon(u), epsilon(v_u)) * dolfin.dx \
    - nu_exp * dolfin.inner(dolfin.grad(u).T * n, v_u) * dolfin.ds(2) \
    - p * dolfin.div(v_u) * dolfin.dx \
    + dolfin.inner(p_outlet_exp * n, v_u) * dolfin.ds(2) \
    + dolfin.div(u) * v_p * dolfin.dx \
    - dolfin.inner(f_exp, v_u) * dolfin.dx
# endregion

# region Create bilinear and linear form.
a = dolfin.lhs(F)
L = dolfin.rhs(F)
# endregion

# region Create files.
file_u = dolfin.File("n_s_3d/u.pvd")
file_p = dolfin.File("n_s_3d/p.pvd")
#file_u_div = dolfin.File("results/u_div.pvd")
# endregion


def solve_div(u, function_space):
    """
    Calculate divergence of a velocity.

    :rtype : dolfin.Function
    :param u: velocity
    :param function_space: function space
    :return: divergence of the velocity
    """

    u_div = dolfin.project(dolfin.div(u), function_space)

    return u_div


def calculate_error(func_s, func_0):
    """
    Calculate error.

    :rtype : double
    :param func_s: current solution
    :param func_0: previous solution
    :return: the relative error between current
    and previous solution
    """

    (func_s_u, func_s_p) = func_s.split()
    (func_0_u, func_0_p) = func_0.split()

    def error(u, u_e):
        """Return the error norm between two functions.

        :rtype : dolfin.Form
        :param u: the input function
        :param u_e: the input function
        """

        return pow((u - u_e), 2.0) * dolfin.dx

    def norm(u):
        """Return the norm of a function.

        :rtype : dolfin.Form
        :param u: the input function
        """

        return pow(u, 2.0) * dolfin.dx

    error_u = dolfin.sqrt(
        dolfin.assemble(error(func_s_u[2], func_0_u[2])))

    error_p = dolfin.sqrt(
        dolfin.assemble(error(func_s_p, func_0_p)))

    norm_u = dolfin.sqrt(dolfin.assemble(norm(func_s_u[2])))
    norm_p = dolfin.sqrt(dolfin.assemble(norm(func_s_p)))

    error_relative_u = error_u / norm_u
    error_relative_p = error_p / norm_p

    print "error_relative_u = {0}; error_relative_p = {1}".format(
        error_relative_u, error_relative_p)

    error_relative_max = max(error_relative_u, error_relative_p)
    return error_relative_max


iter_ = 0

# Transient loop.
while (t < t_end) and (iter_ < max_iter):
    t += dt
    iter_ += 1

    inner_iter_ = 0
    eps = 1.0

    print "t = {0}; iter_ = {1}".format(t, iter_)

    # Inner loop (Picard method).
    while (eps > tol) and (inner_iter_ < max_inner_iter):
        inner_iter_ += 1

        # for i in range(0, max_inner_iter):
        problem = dolfin.LinearVariationalProblem(a, L, up_s, bcs_up)
        solver = dolfin.LinearVariationalSolver(problem)

        solver.parameters["linear_solver"] = method
        solver.parameters["preconditioner"] = preconditioner

        solver.parameters["krylov_solver"]["monitor_convergence"] = \
            krylov_solver_monitor_convergence

        solver.parameters["krylov_solver"]["relative_tolerance"] = \
            krylov_solver_relative_tolerance

        solver.parameters["krylov_solver"]["absolute_tolerance"] = \
            krylov_solver_absolute_tolerance

        solver.solve()

        eps = calculate_error(up_s, up_0)

        up_0.assign(up_s)

        print "inner_iter = {0}; eps = {1}".format(inner_iter_, eps)

    up_1.assign(up_0)

    # Solve divergence.
    #u_div = solve_div(up_0, Q)

    # region Save solution results.
    if is_save_solution:
        file_u << up_0.split()[0]
        file_p << up_0.split()[1]

        #file_u_div << u_div
    # endregion

    # region Plot solution results.
    if is_plot_solution:
        dolfin.plot(up_0.split()[0], title="u", rescale=True, key="u")
        dolfin.plot(up_0.split()[1], title="p", rescale=True, key="p")

        #dolfin.plot(u_div, title="u_div", rescale=True, key="u_div")

        if is_plot_solution_interactive:
            dolfin.interactive()
            # endregion