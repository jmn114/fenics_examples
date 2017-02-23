import numpy as np
import matplotlib.pyplot as plt
from dolfin import *

parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['optimize'] = True
parameters["ghost_mode"] = "shared_facet"

# Create mesh and define function space
#mesh = UnitSquareMesh(3, 3)
mesh = UnitCubeMesh(3, 3, 3)
#mesh = RectangleMesh(Point(0.0, 0.0), Point(1.0, 1.0), 10, 10)

boundary = BoundaryMesh(mesh, "exterior", True)

#print mesh.topology().dim()
#print boundary.topology().dim()
#print mesh.coordinates()
#print mesh.num_cells()
#print mesh.num_vertices()
#print mesh.cells()
#print mesh.hmin()
#print mesh.hmax()

# Structure sub domain
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        #return on_boundary
        return near(x[1], 0.0)
        #return x[1] <= DOLFIN_EPS
        #return x[0] == 0.0 and x[1] == 0.0
        #return x[0] > 0.4 - DOLFIN_EPS and x[0] < 0.6 + DOLFIN_EPS and x[1] < 0.5 + DOLFIN_EPS

# Initialize sub-domain instances
bottom = Bottom()

# Create sub domain markers and mark everything as 0
fluid = MeshFunction("size_t", mesh, mesh.topology().dim())
fluid.set_all(0)

boundaries = MeshFunction("size_t", boundary, boundary.topology().dim())
boundaries.set_all(0)
bottom.mark(boundaries, 1)

# Extract sub meshes
fluid_mesh = SubMesh(mesh, fluid, 0)
#print fluid_mesh.num_vertices()

bottom_mesh = SubMesh(boundary, boundaries, 1)
#print bottom_mesh.num_vertices()
#print bottom_mesh.coordinates()

# Plot meshes
#plot(fluid_mesh, title="Fluid")
#plot(bottom_mesh, title="Bed")
#interactive()

bottom_to_fluid = compute_vertex_map(bottom_mesh, fluid_mesh)
#print(bottom_to_fluid)

#################################################
#Fluid Mesh variational problem
#################################################

V = FunctionSpace(fluid_mesh, 'Lagrange', 1)

# Define boundary conditions
u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')
#u0 = Expression('5*(x[0]+x[1])')
#u0 = Constant(0.0)

def u0_boundary(x, on_boundary):
    return on_boundary

bc_fluid = DirichletBC(V, u0, u0_boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(nabla_grad(u), nabla_grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc_fluid)
#print u(0.5,0.5)

u_nodal_values = u.vector()
u_array = u_nodal_values.array()
#print u_array
#print len(u_array)

u_e = interpolate(u0, V)
u_e_array = u_e.vector().array()
#print 'Max error:', numpy.abs(u_e_array - u_array).max()

#center = (0.5, 0.5)
#print 'numerical u at the center point:',  u(center)
#print 'exact     u at the center point:', u0(center)

file1 = File('heat/fluid_mesh.pvd')
file2 = File('heat/bottom_mesh.pvd')

# Plot solutions and meshes
plot(u, interactive=False)
file1 << u
plot(fluid_mesh, interactive=False, title="Fluid")
plot(bottom_mesh, interactive=False, title="Bed")
file2 << bottom_mesh

################################################
#Bottom Mesh variational problem
################################################

#mydomains = CellFunction('size_t', bottom_mesh)
#mydomains.set_all(0)
#dx_subdomain = Measure('dx')[mydomains]

Pe = Constant(1e10) # K = 1/Pe
T = 2.0 # final time
num_steps = 30 # number of time steps
dt = T / num_steps # time step size
alpha = 3 # parameter alpha
beta = 1.2 # parameter beta

Q = FunctionSpace(bottom_mesh, 'CG', 1)
#Q = FunctionSpace(bottom_mesh, 'DG', 1)

#print bottom_mesh.topology().dim()

#Define initial conditions
#ic= Expression("((pow(x[0]-0.25,2)+pow(x[1]-0.25,2))<0.2*0.2)?(-25*((pow(x[0]-0.25,2)+pow(x[1]-0.25,2))-0.2*0.2)):(0.0)")
#ic= Expression("((pow(x[0]-0.3,2)+pow(x[1]-0.3,2))<0.2*0.2)?(1.0):(0.0)", domain=bottom_mesh)
#ic= Expression("1./(sqrt(2*pi)*0.1) * exp(-0.5*pow((x[0]-0.5),2)/0.01)", domain=bottom_mesh)
ic = Expression("1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t", degree=2, alpha=alpha, beta=beta, t=0)

#Define convection velocity function
#ac = Expression(("-(x[1]-0.5)","(x[0]-0.5)"), domain=bottom_mesh)
#ac = Constant((5.0,0.0))

# Define boundary conditions
def r0_boundary(x, on_boundary):
    return on_boundary

#bc_bottom = DirichletBC(Q, Constant(0.0), r0_boundary, method="geometric")
bc_bottom = DirichletBC(Q, ic, r0_boundary, method="geometric")

'''
r0 = Expression('-5*(x[0])')

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[0] - 1.0) < DOLFIN_EPS and on_boundary

bc_bottom = DirichletBC(Q, r0, r0_boundary)
#bc_bottom = DirichletBC(Q, r0, DomainBoundary())
#bc_bottom = DirichletBC(Q, r0, DirichletBoundary())
'''

# Define variational problem
r = TrialFunction(Q)
s = TestFunction(Q)

# Define initial value
#r_n = interpolate(ic, Q)
r_n = Function(Q)

### HEAT EQUATION ###
g = Constant(beta - 2 - 2*alpha)
F = r*s*dx + dt*dot(grad(r), grad(s))*dx - (r_n + dt*g)*s*dx
a, L = lhs(F), rhs(F)

# Time-stepping
r = Function(Q)
t = 0.0

out_file = File ( "heat/temperature.pvd" )

for n in range(num_steps):
# Update current time
    t += dt
    ic.t = t # update for bc
    # Compute solution
    solve(a == L, r, bc_bottom)
    #solve(lhs(F), rhs(F), r, bc_bottom)
    plot(r, interactive=False)
    out_file << (r, t)
    # Compute error at vertices
    #r_e = interpolate(ic, Q)
    #error = np.abs(r_e.vector().array() - r.vector().array()).max()
    #print("t = %.2f: error = %.3g" % (t, error))
    # Update previous solution
    r_n.assign(r)

#r_nodal_values = r.vector()
#r_array = r_nodal_values.array()
#print r_array
#print len(r_array)

# Plot solutions and meshes
#plot(u, interactive=True)
#plot(fluid_mesh, interactive=True, title="Fluid")
#plot(r, interactive=True)
#plot(bottom_mesh, interactive=True, title="Bed")