from dolfin import *

#H =  0.41
#Lstart = -0.2
#Lend = 2.0
H =  0.40
Lstart = 0.0
Lend = 3.0
# Define mesh
#M = 100
#N = 20
#mesh = Rectangle(0,0,L,H,M,N)

#mesh =  Mesh("cylinder_coarse.xml.gz")
mesh = Mesh("../pipe/pipe.xml")

# Define function spaces
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)

# Define boundary domains
def noslip_boundary(x,on_boundary):
    return (on_boundary and not inflow_boundary(x) and not outflow_boundary(x) and not near(x[1],0))
    #return (near(x[1],0) or near(x[1],H)  
	#    or on_boundary and not inflow_boundary(x) and not outflow_boundary(x))

def inflow_boundary(x):
    return near(x[0],Lstart)

def outflow_boundary(x):
    return near(x[0],Lend)

# Define boundary values
noslip = Constant((0.0,0.0))

#bc_walls = DirichletBC(V, noslip, bndry, 3)
#bc_cylinder = DirichletBC(V, noslip, bndry, 5)

#inflow = Expression(("1.5 * 4.0 * x[1] * (0.41 - x[1]) / ( 0.41 * 0.41 )", "0.0"))
U_m  = 1.5
inflow = Expression(("4*U_m*x[1]*(H-x[1])*sin(DOLFIN_PI*t/8.0)/(H*H)","0"),t=0.0,U_m = U_m,H=H)

bcu1 = DirichletBC(V, inflow, inflow_boundary)
#bcu1  = DirichletBC(V, inflow, bndry, 1)
bcu2 = DirichletBC(V, noslip, noslip_boundary)
bcu = [bcu1, bcu2]
#bcu = [bc_walls, bc_cylinder, bcu1]

bcp1 = DirichletBC(Q, 0.0, outflow_boundary)
#bcp1 = DirichletBC(Q, 0.0, bndry, 4)
bcp = [bcp1]

nu = 0.001
f =  Constant((0.0,0.0))

dt = 0.001
k =  Constant(dt)

# Define functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define test and trial functions
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Tentative velocity step
#F1 = ( (1/k)*inner(u - u0,v) + 0.5*inner(grad(u)*(u0+u),v)
F1 = ( (1/k)*inner(u - u0,v) + inner(grad(u0)*u0,v)
	+ nu*inner(grad(u),grad(v)) - inner(f,v) )*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
F2 = (inner(grad(p),grad(q)) + (1/k)*div(u1)*q)*dx
a2 = lhs(F2)
L2 = rhs(F2)

# Velocity update
F3 = (inner(u - u1,v) + k*inner(grad(p1),v))*dx
a3 = lhs(F3)
L3 = rhs(F3)

# Assemble matrices
A1 = assemble(a1)
b1 = None

A2 = assemble(a2)
b2 = None

A3 = assemble(a3)
b3 = None

# Create files 
u_file = File("results/velocity.pvd")
p_file = File("results/pressure.pvd")

#Solver tweaks
prm = parameters["krylov_solver"]
prm["absolute_tolerance"]= 1.0e-14
prm["relative_tolerance"]= 1.0e-10

# Time stepping
T = 8
t = dt

# Use amg preconditioner for the pressure poisson equation
# if available
precond = "amg" if has_krylov_solver_preconditioner("amg") else "default"

while t <= T + DOLFIN_EPS:
    
    print "Computing time step t = %e" % t

    # Update inflow boundary condition
    inflow.t = t

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1,tensor=b1)
    [ bc.apply(A1,b1) for bc in bcu ]
    solve(A1,u1.vector(),b1,"gmres","default")
    end()

    # Update pressure
    begin("Computing pressure correction")
    b2 = assemble(L2,tensor=b2)
    [ bc.apply(A2,b2) for bc in bcp ]
    solve(A2,p1.vector(),b2,"cg",precond)
    end()

    # Update velocity
    begin("Computing velocity correction")
    b3 = assemble(L3,tensor=b3)
    [ bc.apply(A3,b3) for bc in bcu ]
    solve(A3,u1.vector(),b3,"cg","default")
    end()

    # Plot solution
    plot(p1, title="Pressure", rescale=True)
    plot(u1, title="Velocity", rescale=True)

    # Save to file
    u_file << u1
    p_file << p1

    t += dt
    u0.assign(u1)

interactive()