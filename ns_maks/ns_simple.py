__author__ = 'maks'

"""
Solve rectangular geometry with Navier-Stokes Coupled solver.

Input parameters:
is_slip: bool: use Slip / NoSlip boundary condition for velocity
is_iterative_solver: bool: use / not use iterative solver
is_plot_solution: bool
is_save_solution: bool

f: volume function
folder_name: str: folder where solution will be saved

linear_solver: str
preconditioner: str
"""

from dolfin import *

import numpy

# Define a parameters
#is_problem_3d = True
is_problem_3d = False

is_slip = True
#is_slip = False

is_transient = True
# is_transient = False

#is_iterative_solver = True
is_iterative_solver = False

is_plot_solution = True
# is_plot_solution = False

is_save_solution = True
# is_save_solution = False

is_p_out = True
# is_p_out = False

# p_2 = 0.0
p_2 = 10.0
u_inlet = 1.0

# max_iter = 1
max_iter = 1000
tol = 1E-5

time_end = 1.0

#delta_t = 0.01
delta_t = 0.001
dt = Constant(delta_t)

if is_transient:
    omega = 1.0
else:
    omega = 0.4

linear_solver = "gmres"
#linear_solver = "tfqmr"

#preconditioner = "ilu"
#preconditioner = "amg"
preconditioner = "petsc_amg"

folder_name = "n_s_simple"

# Print the current parameters
print "is_slip = {}".format(is_slip)
print "is_iterative_solver = {}".format(is_iterative_solver)
print "is_plot_solution = {}".format(is_plot_solution)
print "is_save_solution = {}".format(is_save_solution)
print "folder_name = {}".format(folder_name)

if is_iterative_solver:
    print "linear_solver = {}".format(linear_solver)
    print "preconditioner = {}".format(preconditioner)

print "*" * 25

# Define a dolfin parameters
parameters["linear_algebra_backend"] = "PETSc"
#parameters["linear_algebra_backend"] = "uBLAS"

parameters["mesh_partitioner"] = "SCOTCH"
parameters["form_compiler"]["representation"] = "quadrature"
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
dolfin.parameters["form_compiler"]["quadrature_degree"] = 4

# Define a kinematic viscosity
nu = 0.1

# Define a force, mesh
if is_problem_3d:
    f = Constant((0.0, 0.0, 0.0))

    #mesh = dolfin.BoxMesh(0.0, 0.0, 0.0, 1.0, 1.0, 10.0, 4, 4, 8)    
    mesh = dolfin.BoxMesh(dolfin.Point(0.0, 0.0, 0.0), dolfin.Point(1., 1., height), 4, 4, 10)
    # mesh = dolfin.BoxMesh(0.0, 0.0, 0.0, 1.0, 1.0, 10.0, 8, 8, 16)
    # mesh = dolfin.BoxMesh(0.0, 0.0, 0.0, 1.0, 1.0, 10.0, 8, 8, 32)
else:
    f = Constant((0.0, 0.0))

    #mesh = RectangleMesh(0.0, 0.0, 10.0, 1.0, 10, 5, "crossed")
    mesh = RectangleMesh(Point(0, 0), Point(10, 1), 10, 5, 'crossed')
    # mesh = RectangleMesh(0.0, 0.0, 10.0, 1.0, 20, 10, "crossed")
    # mesh = RectangleMesh(0.0, 0.0, 10.0, 1.0, 40, 20, "crossed")

info(mesh)

n = FacetNormal(mesh)

# plot(mesh)
# interactive()

# Define function spaces
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
VQ = V * Q

# Define functions
up = TrialFunction(VQ)
u, p = split(up)

v_u, v_p = TestFunctions(VQ)

up_1 = Function(VQ)
up_0 = Function(VQ)
u_0, p_0 = split(up_0)


if is_problem_3d:
    up_expr = Constant((0.0, 0.0, 1.0, 0.0))
else:
    up_expr = Constant((1.0, 0.0, 0.0))

if is_transient:
    up_0 = interpolate(up_expr, VQ)
    u_0, p_0 = split(up_0)


# Define Sub domain (Walls)
class Walls(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


# Define Sub domain (inflow)
class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        if is_problem_3d:
            return on_boundary and dolfin.near(x[2], 0.0)
        else:
            return on_boundary and dolfin.near(x[0], 0.0)


# Define Sub domain (outflow)
class Outlet(SubDomain):
    def inside(self, x, on_boundary):
        if is_problem_3d:
            return on_boundary and dolfin.near(x[2], 10.0)
        else:
            return on_boundary and dolfin.near(x[0], 10.0)


# Define mesh boundaries
domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

domains.set_all(0)

walls = Walls()
inlet = Inlet()
outlet = Outlet()

walls.mark(domains, 3)
inlet.mark(domains, 1)
outlet.mark(domains, 2)

# plot(domains)
# interactive()

ds = dolfin.Measure("ds")[domains]

# Save sub domains to file
# file = File(folder_name + "/domains.pvd")
# file << domains

# Define boundary conditions
if is_problem_3d:
    u_inlet_const = Constant((0.0, 0.0, u_inlet))
    u_wall_const = Constant((0.0, 0.0, 0.0))
else:
    u_inlet_const = Constant((u_inlet, 0.0))
    u_wall_const = Constant((0.0, 0.0))

p_2_const = Constant(p_2)

# Inlet BC
bc_u_1 = DirichletBC(VQ.sub(0), u_inlet_const, domains, 1)
# Wall BC
bc_u_3 = DirichletBC(VQ.sub(0), u_wall_const, domains, 3)

# bcs_up = []

bcs_up = None
F = None

if not is_slip:
    bcs_up = [bc_u_1, bc_u_3]
else:
    bcs_up = [bc_u_1]

F = nu * inner(grad(u) + grad(u).T, grad(v_u)) * dx \
    - p * div(v_u) * dx \
    + div(u) * v_p * dx \
    - inner(f, v_u) * dx

# Add convective term
F += inner(dot(grad(u), u_0), v_u) * dx


def t(u, p):
    return dot(2.0 * nu * sym(grad(u)), n) - p * n

# Add Slip boundary condition
if is_slip:
    beta = Constant(10.0)

    cell_size = CellSize(mesh)

    F += - dot(n, t(u, p)) * dot(v_u, n) * ds(3) \
         - dot(u, n) * dot(n, t(v_u, v_p)) * ds(3) \
         + beta / cell_size * dot(u, n) * dot(v_u, n) * ds(3)

if is_p_out:
    # Weak form of p_out bc.
    F += inner(p_2_const * n, v_u) * ds(2)

if is_transient:
    F += (1 / dt) * inner(u - u_0, v_u) * dx

timer_solver_all = Timer("TimerSolveAll")
timer_solver_all.start()

# Create bilinear and linear form
a = lhs(F)
L = rhs(F)

# Define a parameters for a stationary loop
eps = 1.0
iter_ = 0
time = 0.0


def solve():
    problem = LinearVariationalProblem(a, L, up_1, bcs_up)
    solver = LinearVariationalSolver(problem)

    if is_iterative_solver:
        solver.parameters["linear_solver"] = linear_solver
        solver.parameters["preconditioner"] = preconditioner

    solver.solve()

    diff_up = up_1.vector().array() - up_0.vector().array()
    eps = numpy.linalg.norm(diff_up, ord=numpy.Inf)

    if omega == 1.0:
        up_0.assign(up_1)
    else:
        up_1.vector().axpy(-1, up_0.vector())
        up_0.vector().axpy(omega, up_1.vector())

    return eps

if is_transient:
    while (time < time_end and
           iter_ < max_iter and
           eps > tol):

        iter_ += 1
        time += delta_t

        eps = solve()

        print "iter = {:d}; eps_up = {:e}; time = {:.2f}\n".format(iter_, eps, time)

else:
    # Solve stationary Navier-Stokes problem with Picard method
    while (iter_ < max_iter and
            eps > tol):

        iter_ += 1

        eps = solve()

        print "iter = {:d}; eps_up = {:e}\n".format(iter_, eps)

timer_solver_all.stop()
info("TimerSolveAll = %g" % timer_solver_all.value())

u, p = up_0.split()

print "*" * 25
print "end"

# Plot a solution
if is_plot_solution:
    plot(u, title="u")
    plot(p, title="p")
    interactive()

# Save a solution to a file
if is_save_solution:
    File(folder_name + "/" + "u.pvd") << u
    File(folder_name + "/" + "p.pvd") << p