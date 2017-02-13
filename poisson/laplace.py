import dolfin as df

# The geometry
start_x, end_x, n_elem_x = 0.0, 2.0, 4
start_y, end_y, n_elem_y = 0.0, 1.0, 2

# Make a mesh
mesh = df.RectangleMesh(start_x, start_y, end_x, end_y, n_elem_x, n_elem_y)

def define_boundary(defining_func, boundary_parts, boundary_id):
    """
    Define a boundary

    - defining_func is used instead of the normal SubDomain inside method: defining_func(x, on_boundary)
    - boundary_parts is a MeshFunction used to distinguish parts of the boundary from each other
    - boundary_id is the id of the new boundary in the boundary_parts dictionary
    """
    class Boundary(df.SubDomain):
        def inside(self, x, on_boundary):
            return defining_func(x, on_boundary)
    boundary = Boundary()
    boundary.mark(boundary_parts, boundary_id)
    return boundary

# Each part of the mesh gets its own ID
ID_LEFT_SIDE, ID_RIGHT_SIDE, ID_BOTTOM, ID_TOP, ID_OTHER = 1, 2, 3, 4, 5
boundary_parts = df.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundary_parts.set_all(ID_OTHER)

# Define the boundaries
GEOM_TOL = 1e-6
define_boundary(lambda x, on_boundary: on_boundary and x[0] < start_x + GEOM_TOL, boundary_parts, ID_LEFT_SIDE)
define_boundary(lambda x, on_boundary: on_boundary and x[0] > end_x - GEOM_TOL, boundary_parts, ID_RIGHT_SIDE)
define_boundary(lambda x, on_boundary: on_boundary and x[1] < start_y + GEOM_TOL, boundary_parts, ID_BOTTOM)
define_boundary(lambda x, on_boundary: on_boundary and x[1] > end_y - GEOM_TOL, boundary_parts, ID_TOP)

# Define a function space of bi-linear elements on our mesh
V = df.FunctionSpace(mesh, 'Lagrange', 1)

# The trial and test functions are defined on the bi-linear function space
u = df.TrialFunction(V)
v = df.TestFunction(V)

# Our equation is defined as a=L. We only define what is inside the integrals, and the difference
# between integrals over the domain and the boundary is whether we use dx or ds:
#   - the term a contains both the trial (unknown) function u and the test function v
#   - the term L contains only test function v
a = df.inner(df.nabla_grad(u), df.nabla_grad(v))*df.dx
L = df.Constant(42)*v*df.ds(ID_BOTTOM) + df.Constant(42)*v*df.ds(ID_TOP)

# The Dirichlet boundary conditions
expr_left = df.Expression('-1')
expr_right = df.Expression('+1')
bc_left = df.DirichletBC(V, expr_left, boundary_parts, ID_LEFT_SIDE)
bc_right = df.DirichletBC(V, expr_right, boundary_parts, ID_RIGHT_SIDE)
boundary_conditions = [bc_left, bc_right]

# The function we want to calculate, unfortunately called u as well by convention
u = df.Function(V)

# Solve!
df.solve(a == L, u, boundary_conditions)

# Show the mesh and the result
df.plot(mesh)
df.plot(u)
df.interactive()