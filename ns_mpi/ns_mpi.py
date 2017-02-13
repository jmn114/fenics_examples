from dolfin import *
#import mshr

comm = mpi_comm_world()
rank = MPI.rank(comm)
set_log_level(INFO if rank==0 else INFO+1)
parameters["std_out_all_processes"] = False



# Define domain
center = Point(1.0, 0.15)
radius = 0.05
L = 3.0
W = 0.4

'''

geometry = mshr.Rectangle(Point(0.0, 0.0), Point(L, W)) \
         - mshr.Circle(center, radius, 10)

# Build mesh
mesh = mshr.generate_mesh(geometry, 50)

'''

mesh = Mesh("../pipe/pipe.xml")

# Construct facet markers
bndry = FacetFunction("size_t", mesh)
for f in facets(mesh):
    mp = f.midpoint()
    if near(mp[0], 0.0): # inflow
        bndry[f] = 1
    elif near(mp[0], L): # outflow
        bndry[f] = 2
    elif near(mp[1], 0.0) or near(mp[1], W): # walls
        bndry[f] = 3
    elif mp.distance(center) <= radius: # cylinder
        bndry[f] = 5



#mesh = Mesh("../pipe/pipe.xml")
#bndry = MeshFunction('size_t', mesh, "../pipe/pipe_facet_region.xml")
#plot(bndry, interactive=True)

# Build function spaces (Taylor-Hood)
V = VectorFunctionSpace(mesh, "CG", 2)
P = FunctionSpace(mesh, "CG", 1)
E = FunctionSpace(mesh, "CG", 1)
W = MixedFunctionSpace([V, P, E])

# No-slip boundary condition for velocity on walls and cylinder - boundary id 3
noslip = Constant((0, 0))
bcv_walls = DirichletBC(W.sub(0), noslip, bndry, 3)

vc= Expression(("-0.5*t*cos(atan2(x[0]-0.2,x[1]-0.2))","0.5*t*sin(atan2(x[0]-0.2,x[1]-0.2))"),t=0)
bcv_cylinder = DirichletBC(W.sub(0), vc, bndry, 5)

bce_cylinder = DirichletBC(W.sub(2), Constant(1.0), bndry, 5)
# Inflow boundary condition for velocity - boundary id 1
v_in = Expression(("1.5 * 4.0 * x[1] * (0.41 - x[1]) / ( 0.41 * 0.41 )", "0.0"))
bcv_in = DirichletBC(W.sub(0), v_in, bndry, 1)

# Collect boundary conditions
bcs = [bcv_cylinder, bcv_walls, bcv_in, bce_cylinder]

# Facet normal, identity tensor and boundary measure
n = FacetNormal(mesh)
I = Identity(mesh.geometry().dim())
ds = Measure("ds", subdomain_data=bndry)
nu = Constant(0.001)

dt = 0.1
t_end = 10
theta=0.5   # Crank-Nicholson timestepping
k=0.01

# Define unknown and test function(s)
(v_, p_, e_) = TestFunctions(W)

# current unknown time step
w = Function(W)
(v, p, e) = split(w)

# previous known time step
w0 = Function(W)
(v0, p0, e0) = split(w0)

def a(v,u) :
    D = sym(grad(v))
    return (inner(grad(v)*v, u) + inner(2*nu*D, grad(u)))*dx

def b(q,v) :
    return inner(div(v),q)*dx

def c(v,e,g) :
    return ( inner(k*grad(e),grad(g)) + inner(v,grad(e))*g )*dx

# Define variational forms without time derivative in previous time
F0_eq1 = a(v0,v_) + b(p,v_)
F0_eq2 = b(p_,v)
F0_eq3 = c(v0,e0,e_)
F0 = F0_eq1 + F0_eq2 + F0_eq3

# variational form without time derivative in current time
F1_eq1 = a(v,v_) + b(p,v_)
F1_eq2 = b(p_,v)
F1_eq3 = c(v,e,e_)
F1 = F1_eq1 + F1_eq2 + F1_eq3

#combine variational forms with time derivative
#
#  dw/dt + F(t) = 0 is approximated as
#  (w-w0)/dt + (1-theta)*F(t0) + theta*F(t) = 0
#
F = (inner((v-v0),v_)/dt + inner((e-e0),e_)/dt)*dx + (1.0-theta)*F0 + theta*F1

# Create files for storing solution
name="a"
vfile = XDMFFile(mpi_comm_world(),"results_%s/v.xdmf" % name)
pfile = XDMFFile(mpi_comm_world(),"results_%s/p.xdmf" % name)
efile = XDMFFile(mpi_comm_world(),"results_%s/e.xdmf" % name)
vfile.parameters["flush_output"] = True
pfile.parameters["flush_output"] = True
efile.parameters["flush_output"] = True


J = derivative(F, w)
problem=NonlinearVariationalProblem(F,w,bcs,J)
solver=NonlinearVariationalSolver(problem)

prm = solver.parameters
#info(prm,True)  #get full info on the parameters
prm['nonlinear_solver'] = 'newton'
prm['newton_solver']['absolute_tolerance'] = 1E-12
prm['newton_solver']['relative_tolerance'] = 1e-12
prm['newton_solver']['maximum_iterations'] = 20
prm['newton_solver']['linear_solver'] = 'mumps'


# Time-stepping
t = dt
while t < t_end:

    print "t =", t
    #vc.t=t

    # Compute
    begin("Solving ....")
    solver.solve()
    end()

    # Extract solutions:
    (v, p, e) = w.split()

    v.rename("v", "velocity")
    p.rename("p", "pressure")
    e.rename("e", "temperature")
    # Save to file
    vfile << v
    pfile << p
    efile << e

    # Move to next time step
    w0.assign(w)
    t += dt