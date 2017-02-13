from dolfin import *

# Defining mesh, boundaries and regions
# This assumes some preprocessing (see explanation)
mesh = Mesh('unitsquare.xml')
boundaries = MeshFunction('size_t', mesh, 'unitsquare_facet_region.xml')
regions = MeshFunction('size_t', mesh, 'unitsquare_physical_region.xml')
dx = Measure('dx')[regions]
ds = Measure('ds')[boundaries]

# Define function space and basis functions
V = FunctionSpace(mesh, "Lagrange", 1)
u = TrialFunction(V)
v = TestFunction(V)

# Define Dirichlet boundary conditions
# The boundary labeled 3 can be uncommented to
# demand a discontinuous jump on the boundary
u0, u1 = 0.0, 0.05
bcs = [DirichletBC(V, Constant(u0), boundaries, 1),
       # DirichletBC(V, Constant(u0), boundaries, 3),
       DirichletBC(V, Constant(u1), boundaries, 2)]

# Define weak form and solve
# Although not used here, dx(i) and ds(i) could be used
# to specify an integral over a specific region/boundary
a = inner(grad(u), grad(v)) * dx
L = Constant(0.0) * v * dx
u = Function(V)
problem = LinearVariationalProblem(a, L, u, bcs)
solver = LinearVariationalSolver(problem)
solver.solve()
E = -project(grad(u), VectorFunctionSpace(mesh, "Lagrange", 1))

# Plotting
plot(mesh)
plot(boundaries)
plot(regions)
plot(u)
plot(E)
interactive()
