import time
import os
import math
import numpy as np
from dolfin import *

def iter(c):

    # get file name
    #fileName = os.path.splitext(__file__)[0]
    
    # Create mesh and define function space
    #mesh = RectangleMesh(Point(0, 0), Point(1, 1), 40, 40, 'crossed')
    #mesh = UnitSquareMesh(10, 10)
    #mesh = UnitCubeMesh(10, 10, 10)
    
    #mesh = Mesh('unitsquaremesh.xml')
    mesh = Mesh('../mesh_pipe/pipe_coarse.xml')
    ##print mesh
    ##print mesh.topology().dim()
    
    ##print mesh.num_cells()
    ##print mesh.num_vertices()
    #print mesh.coordinates()
    #print mesh.cells()
    #print mesh.hmin()
    #print mesh.hmax()
    
    boundary = BoundaryMesh(mesh, "exterior", True)
    ##print boundary
    ##print boundary.topology().dim()
    #File("boundarymesh.xml") << boundary
    
    # bc labels for unitsquare
    #bottom_bdry = 1
    #top_bdry = 2
    #inlet_bdry = 3
    #outlet_bdry = 4
    #internal = 5
    
    # bc labes for pipe 2d
    bottom_bdry = 3
    top_bdry = 2
    inlet_bdry = 1
    outlet_bdry = 4
    cylinder = 5
    #internal = 14
    
    #boundaries = MeshFunction('size_t', mesh, 'unitsquaremesh_facet_region.xml')
    boundaries = MeshFunction('size_t', mesh, '../mesh_pipe/pipe_coarse_facet_region.xml')
    ##print boundaries
    ##print boundaries.array()
    ds_p = Measure('ds')[boundaries]
    
    #regions = MeshFunction('size_t', mesh, 'unitsquaremesh_physical_region.xml')
    regions = MeshFunction('size_t', mesh, '../mesh_pipe/pipe_coarse_physical_region.xml')
    dx_p = Measure('dx')[regions]
    
    '''
    boundary_ids = boundaries.array()    
    facet_nodes = np.array([f.entities(0) for f in facets(mesh)])
    facet_indices = [f.index() for f in facets(mesh)]
        
    # create list of nodes on pointwise boundary
    ptws_bndy_nds = []
    for j in range(len(facet_nodes)):
        if boundary_ids[j] == bottom_bdry:
            ptws_bndy_nds.extend(list(facet_nodes[j]))
    ptws_bndy_nds = sorted(set(ptws_bndy_nds))
     
    print ptws_bndy_nds
    print len(ptws_bndy_nds)
    
    # define a submesh, composed of bottom boundary
    bottom = FacetFunction('size_t', boundary)
    print bottom
    print bottom.array()
    bottom.set_all(0)
    bottom.array()[ptws_bndy_nds] = 1
    #bottom.array()[boundaries.array() == 1] = 1
    print bottom.array()
    bottom_mesh = SubMesh(boundary, bottom, 1)
    print bottom_mesh
    '''
    
    
    # copy the meshfunction, but defined on the boundarymesh
    bdim = boundary.topology().dim()
    boundary_boundaries = MeshFunction('size_t', boundary, bdim)
    boundary_boundaries.set_all(0)
    for i, facet in enumerate(entities(boundary, bdim)):
        parent_meshentity = boundary.entity_map(bdim)[i]
        parent_boundarynumber = boundaries.array()[parent_meshentity]
        boundary_boundaries.array()[i] = parent_boundarynumber
    
    bottom_mesh = SubMesh(boundary, boundary_boundaries, bottom_bdry)
    ##print bottom_mesh
    
    ##fluid_mesh = SubMesh(mesh, regions, internal)
    ##print fluid_mesh
    
    
    '''
    # Structure sub domain
    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            #return on_boundary
            return near(x[1], 0.0)
            #return on_boundary and not (near(x[0], 0.0) or near(x[0], 1.0) or near(x[1], 1.0))
            #return x[1] <= DOLFIN_EPS
            #return x[0] == 0.0 and x[1] == 0.0
            #return x[0] > 0.4 - DOLFIN_EPS and x[0] < 0.6 + DOLFIN_EPS and x[1] < 0.5 + DOLFIN_EPS
    
    # Initialize sub-domain instances
    bottom = Bottom()
    
    # Create sub domain markers and mark everything as 0
    fluid = MeshFunction("size_t", mesh, mesh.topology().dim())
    fluid.set_all(0)
    bottom.mark(fluid, 1)
    
    boundaries = MeshFunction("size_t", boundary, boundary.topology().dim())
    boundaries.set_all(0)
    bottom.mark(boundaries, 1)
    
    # Extract sub meshes
    fluid_mesh = SubMesh(mesh, fluid, 0)
    bottom_mesh = SubMesh(boundary, boundaries, 1)
    '''
    
    
    ##bottom_to_fluid = compute_vertex_map(bottom_mesh, fluid_mesh)
    ##print(bottom_to_fluid)
    
    #plot(mesh, interactive = False, title="Mesh")
    #plot(boundary, interactive = False, title="Boundary")
    #plot(fluid_mesh, interactive = False, title="Fluid")
    #plot(bottom_mesh, interactive = False, title="Bed")
    
    #File("mesh.pvd") << mesh
    #File("bottom_mesh.pvd") << bottom_mesh 
    
    # Hold plot
    #interactive()
    
    #exit()
    
    ### NAVIER-STOKES IN PARENT MESH###
    
    ### N-S CHORIN'S SPLITTING METHOD ###
    '''
    T = 10.0 # final time
    num_steps = 20 # number of time steps
    dt = T / num_steps # time step size
    mu = 1 # kinematic viscosity
    rho = 1 # density
    
    # Create mesh and define function spaces
    V = VectorFunctionSpace(mesh, "CG", 1)
    Q = FunctionSpace(mesh, "CG", 1)
    
    # Define boundaries
    inflow = "near(x[0], 0)"
    outflow = "near(x[0], 1)"
    walls = "near(x[1], 0) || near(x[1], 1)"
    
    # Define boundary conditions
    # 2D #
    bcu_noslip = DirichletBC(V, Constant((0, 0)), walls)
    # 3D #
    #bcu_noslip = DirichletBC(V, Constant((0, 0, 0)), walls)
    bcp_inflow = DirichletBC(Q, Constant(8), inflow)
    bcp_outflow = DirichletBC(Q, Constant(0), outflow)
    bcu = [bcu_noslip]
    bcp = [bcp_inflow, bcp_outflow]
    
    # Define trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)
    p = TrialFunction(Q)
    q = TestFunction(Q)
    
    # Define functions for solutions at previous and current time steps
    u_n = Function(V)
    u_ = Function(V)
    p_n = Function(Q)
    p_ = Function(Q)
    
    # Define expressions used in variational forms
    U = 0.5*(u_n + u)
    n = FacetNormal(mesh)
    # 2D #
    f = Constant((0, 0))
    # 3D #
    #f = Constant((0, 0, 0))
    k = Constant(dt)
    mu = Constant(mu)
    rho = Constant(rho)
    
    # Define strain-rate tensor
    def epsilon(u):
        return sym(nabla_grad(u))
    
    # Define stress tensor
    def sigma(u, p):
        return 2*mu*epsilon(u) - p*Identity(len(u))
    
    # Define variational problem for step 1
    F1 = rho*dot((u - u_n) / k, v)*dx + \
    rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
    + inner(sigma(U, p_n), epsilon(v))*dx \
    + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
    - rho*dot(f, v)*dx
    a1 = lhs(F1)
    L1 = rhs(F1)
    
    # Define variational problem for step 2
    a2 = dot(nabla_grad(p), nabla_grad(q))*dx
    L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx
    
    # Define variational problem for step 3
    a3 = dot(u, v)*dx
    L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx
    
    # Assemble matrices
    A1 = assemble(a1)
    A2 = assemble(a2)
    A3 = assemble(a3)
    
    # Apply boundary conditions to matrices
    [bc.apply(A1) for bc in bcu]
    [bc.apply(A2) for bc in bcp]
    
    # Create VTK files for visualization output
    vtkfile_u = File("ns_chorin/velocity.pvd")
    vtkfile_p = File("ns_chorin/pressure.pvd")
    
    # Create time series for saving solution for later
    #timeseries_u = TimeSeries("n_s_chorin/velocity")
    #timeseries_p = TimeSeries("n_s_chorin/pressure")
    
    # Save mesh to file for later
    File("newmesh.xml") << mesh
    
    # Time-stepping
    t = 0
    for n in range(num_steps):
        # Update current time
        t += dt
        # Step 1: Tentative velocity step
        b1 = assemble(L1)
        [bc.apply(b1) for bc in bcu]
        solve(A1, u_.vector(), b1)
        # Step 2: Pressure correction step
        b2 = assemble(L2)
        [bc.apply(b2) for bc in bcp]
        solve(A2, p_.vector(), b2)
        # Step 3: Velocity correction step
        b3 = assemble(L3)
        solve(A3, u_.vector(), b3)
        # Plot solution
        plot(u_)
        # Save solution to file (VTK)
        vtkfile_u << (u_, t)
        vtkfile_p << (p_, t)
        # Save solution to file (HDF5)
        #timeseries_u.store(u_.vector(), t)
        #timeseries_p.store(p_.vector(), t)   
        # Compute error
        #u_e = Expression(("4*x[1]*(1.0 - x[1])", "0"), degree=2)
        #u_e = interpolate(u_e, V)
        #error = np.abs(u_e.vector().array() - u_.vector().array()).max()
        #print("t = %.2f: error = %.3g" % (t, error))
        print("max u:", u_.vector().array().max())
        # Update previous solution
        u_n.assign(u_)
        p_n.assign(p_)
    
    # Hold plot
    interactive()
    '''
    ### N-S NON-LINEAR VARIATIONAL PROBLEM ###
    '''
    # Build function spaces (Taylor-Hood)
    V = VectorFunctionSpace(mesh, "CG", 1)
    P = FunctionSpace(mesh, "CG", 1)
    W = MixedFunctionSpace([V, P])
    
    # Define boundaries
    inflow = "near(x[0], 0)"
    outflow = "near(x[0], 1)"
    walls = "near(x[1], 0) || near(x[1], 1)"
    
    # Define boundary conditions
    # 2D #
    bc_walls = DirichletBC(W.sub(0), Constant((0, 0)), walls)
    # 3D #
    #bc_walls = DirichletBC(W.sub(0), Constant((0, 0, 0)), walls)
    
    # Inflow boundary condition for velocity
    #v_in = Expression(("0.3 * 4.0 * x[1] * (0.41 - x[1]) / ( 0.41 * 0.41 )", "0.0"))
    v_in = Expression(("1.0", "0.0"))
    bc_in = DirichletBC(W.sub(0), v_in, inflow)
    
    # Pressure BCs
    #bcp_inflow = DirichletBC(W.sub(1), Constant(8), inflow)
    #bcp_outflow = DirichletBC(W.sub(1), Constant(0), outflow)
    
    # Collect boundary conditions
    bcs = [bc_walls, bc_in]
    
    # Facet normal, identity tensor and boundary measure
    #n = FacetNormal(mesh)
    #ds = Measure("ds", subdomain_data=bndry)
    I = Identity(mesh.geometry().dim())
    nu = Constant(0.001)
    
    dt = 0.1
    t_end = 10
    theta = 0.5   # Crank-Nicholson timestepping
    
    # Define unknown and test function(s)
    (v_, p_) = TestFunctions(W)
    
    # current unknown time step
    w = Function(W)
    (v, p) = split(w)
    
    # previous known time step
    w0 = Function(W)
    (v0, p0) = split(w0)
    
    def a(v,u) :
        #D = 0.5*(grad(v)+grad(v).T)
        D = sym(grad(v))
        #T = -p*I + 2*nu*D    
        #return (inner(grad(v)*v, u) + inner(T, grad(u)))*dx
        return (inner(grad(v)*v, u) + inner(2*nu*D, grad(u)))*dx
    
    def b(q,v) :
        return inner(div(v),q)*dx
    
    # Define variational forms without time derivative in previous time
    F0_eq1 = a(v0,v_) + b(p,v_)
    F0_eq2 = b(p_,v)
    F0 = F0_eq1 + F0_eq2
    
    # Define variational forms without time derivative in current time
    F1_eq1 = a(v,v_) + b(p,v_)
    F1_eq2 = b(p_,v)
    F1 = F1_eq1 + F1_eq2
    
    #combine variational forms with time derivative
    #
    #  dw/dt + F(t) = 0 is approximated as
    #  (w-w0)/dt + (1-theta)*F(t0) + theta*F(t) = 0
    #
    F = (inner((v-v0),v_)/dt)*dx + (1.0-theta)*F0 + theta*F1
    
    # Create files for storing solution
    #name="a"
    #vfile = XDMFFile(mpi_comm_world(),"results_%s/v.xdmf" % name)
    #pfile = XDMFFile(mpi_comm_world(),"results_%s/p.xdmf" % name)
    #efile = XDMFFile(mpi_comm_world(),"results_%s/e.xdmf" % name)
    #vfile.parameters["flush_output"] = True
    #pfile.parameters["flush_output"] = True
    #efile.parameters["flush_output"] = True
    
    vtkfile_u = File("ns_nlvs/velocity.pvd")
    vtkfile_p = File("ns_nlvs/pressure.pvd")
    
    J = derivative(F, w)
    ffc_options = {"optimize": True, "quadrature_degree": 8}
    problem = NonlinearVariationalProblem(F, w, bcs, J, form_compiler_parameters=ffc_options)
    solver = NonlinearVariationalSolver(problem)
    
    prm = solver.parameters
    #info(prm,True)  #get full info on the parameters
    prm['nonlinear_solver'] = 'newton'
    prm['newton_solver']['absolute_tolerance'] = 1e-12
    prm['newton_solver']['relative_tolerance'] = 1e-12
    prm['newton_solver']['maximum_iterations'] = 50
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
        (v, p) = w.split()
    
        v.rename("v", "velocity")
        p.rename("p", "pressure")
        
        # Plot solution
        plot(v, interactive = True)
        
        # Save to file
        vtkfile_u << v
        vtkfile_p << p
    
        # Move to next time step
        w0.assign(w)
        t += dt
    
    # Hold plot
    interactive()
    '''
    ### N-S LINEAR VARIATIONAL PROBLEM ###
    
    # Define a parameters
    #is_problem_3d = True
    is_problem_3d = False
    
    is_slip = True
    #is_slip = False
    
    is_transient = True
    #is_transient = False
    
    #is_iterative_solver = True
    is_iterative_solver = False
    
    is_plot_solution = True
    #is_plot_solution = False
    
    is_save_solution = True
    #is_save_solution = False
    
    is_p_out = True
    #is_p_out = False
    
    p_2 = 0.0
    #p_2 = 10.0
    u_inlet = 1.0
    
    # max_iter = 1
    max_iter = 1000
    tol = 1E-5
    
    time_end = 1.0
    
    delta_t = 0.1
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
    
    folder_name = "ns"
    
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
    parameters["std_out_all_processes"] = False
    parameters["form_compiler"]["quadrature_degree"] = 8
    
    # Define a kinematic viscosity
    #nu = 0.1
    nu = 1E-4
    
    # Define a force, mesh
    if is_problem_3d:
        f = Constant((0.0, 0.0, 0.0))
    
        #mesh = dolfin.BoxMesh(0.0, 0.0, 0.0, 1.0, 1.0, 10.0, 4, 4, 8)    
        #mesh = dolfin.BoxMesh(dolfin.Point(0.0, 0.0, 0.0), dolfin.Point(1., 1., height), 4, 4, 10)
        #mesh = dolfin.BoxMesh(0.0, 0.0, 0.0, 1.0, 1.0, 10.0, 8, 8, 16)
        #mesh = dolfin.BoxMesh(0.0, 0.0, 0.0, 1.0, 1.0, 10.0, 8, 8, 32)
    else:
        f = Constant((0.0, 0.0))
    
        #mesh = RectangleMesh(0.0, 0.0, 10.0, 1.0, 10, 5, "crossed")
        #mesh = RectangleMesh(Point(0, 0), Point(10, 1), 10, 5, 'crossed')
        #mesh = RectangleMesh(0.0, 0.0, 10.0, 1.0, 20, 10, "crossed")
        #mesh = RectangleMesh(0.0, 0.0, 10.0, 1.0, 40, 20, "crossed")
    
    #info(mesh)
    #plot(mesh)
    #interactive()
    
    # Define function spaces
    #V = VectorFunctionSpace(mesh, "CG", 2)
    V = VectorFunctionSpace(mesh, "CG", 1)
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
    #class Walls(SubDomain):
    #    def inside(self, x, on_boundary):
    #        return on_boundary
    
    # Define Sub domain (inflow)
    #class Inlet(SubDomain):
    #    def inside(self, x, on_boundary):
    #        if is_problem_3d:
    #            return on_boundary and dolfin.near(x[2], 0.0)
    #        else:
    #            return on_boundary and dolfin.near(x[0], 0.0)
    
    # Define Sub domain (outflow)
    #class Outlet(SubDomain):
    #    def inside(self, x, on_boundary):
    #        if is_problem_3d:
    #            return on_boundary and dolfin.near(x[2], 10.0)
    #        else:
    #            return on_boundary and dolfin.near(x[0], 10.0)
    
    # Define mesh boundaries
    #domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    
    #domains.set_all(0)
    
    #walls = Walls()
    #inlet = Inlet()
    #outlet = Outlet()
    
    #walls.mark(domains, 3)
    #inlet.mark(domains, 1)
    #outlet.mark(domains, 2)
    
    # plot(domains)
    # interactive()
    
    #ds = dolfin.Measure("ds")[domains]
    
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
    #bc_u_1 = DirichletBC(VQ.sub(0), u_inlet_const, domains, 1)
    bc_u_in = DirichletBC(VQ.sub(0), u_inlet_const, boundaries, 1)
    # Wall BC
    #bc_u_3 = DirichletBC(VQ.sub(0), u_wall_const, domains, 3)
    bc_u_bottom = DirichletBC(VQ.sub(0), u_wall_const, boundaries, 3)
    bc_u_top = DirichletBC(VQ.sub(0), u_wall_const, boundaries, 2)
    bc_u_cylinder = DirichletBC(VQ.sub(0), u_wall_const, boundaries, 5)
    
    # bcs_up = []
    
    bcs_up = None
    F = None
    
    n = FacetNormal(mesh)
    
    if not is_slip:
        #bcs_up = [bc_u_1, bc_u_3]
        bcs_up = [bc_u_in, bc_u_bottom, bc_u_top, bc_u_cylinder]
    else:
        bcs_up = [bc_u_in, bc_u_top, bc_u_cylinder]
    
    F = nu * inner(grad(u) + grad(u).T, grad(v_u)) * dx_p \
        - p * div(v_u) * dx_p \
        + div(u) * v_p * dx_p \
        - inner(f, v_u) * dx_p
    
    # Add convective term
    F += inner(dot(grad(u), u_0), v_u) * dx_p
    
    def t(u, p):
        return dot(2.0 * nu * sym(grad(u)), n) - p * n
    
    # Add Slip boundary condition
    if is_slip:
        beta = Constant(10.0)
    
        cell_size = CellSize(mesh)
    
        F += - dot(n, t(u, p)) * dot(v_u, n) * ds_p(3) \
             - dot(u, n) * dot(n, t(v_u, v_p)) * ds_p(3) \
             + beta / cell_size * dot(u, n) * dot(v_u, n) * ds_p(3)
    
    if is_p_out:
        # Weak form of p_out bc.
        F += inner(p_2_const * n, v_u) * ds_p(4)
    
    if is_transient:
        F += (1 / dt) * inner(u - u_0, v_u) * dx_p
    
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
        eps = np.linalg.norm(diff_up, ord=np.Inf)
    
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
        while (iter_ < max_iter and eps > tol):
    
            iter_ += 1
    
            eps = solve()
    
            print "iter = {:d}; eps_up = {:e}\n".format(iter_, eps)
    
    timer_solver_all.stop()
    info("TimerSolveAll = %g" % timer_solver_all.value())
    
    u, p = up_0.split()
    
    print "*" * 25
    
    # Hold plot
    #interactive()
    
    #vel, pre = up_0.split(deepcopy=True)
    
    vel = Function(u)
    pre = Function(p)
    
    ##coordinates = mesh.coordinates()
    ##print coordinates
    ##print len(coordinates)
    ##print coordinates[0]
    ##print coordinates[0][1]
    ##print vel(coordinates[0])
    
    nodal_values_vel = vel.vector().array()
    nodal_values_pre = pre.vector().array()
    ##print nodal_values_vel
    ##print len(nodal_values_vel)
    #print nodal_values_pre
    #print len(nodal_values_pre)
    
    ##vertex_values = vel.compute_vertex_values()
    ##print vertex_values
    ##print len(vertex_values)
    
    ##for i, x in enumerate(coordinates):        
    ##    print("vertex %d: \t u(%s) = %s" % (i, x, vel(x)))    
    
    ##v2d = vertex_to_dof_map(V)
    ##print v2d
    ##print nodal_values_vel[v2d]
    
    ##d2v = dof_to_vertex_map(V)
    ##print d2v
    #print coordinates[d2v]
    
    #element = V.element()
    #print element
    #dofmap = V.dofmap()
    #print dofmap
    #for cell in cells(mesh):
    #    print(element.tabulate_dof_coordinates(cell))
    #    print(dofmap.cell_dofs(cell.index()))
    
    print "*" * 25
    
    # Plot a solution
    ##if is_plot_solution:
        ##plot(u, title="u")
        ##plot(p, title="p")
        ##interactive()
    
    # Save a solution to a file
    if is_save_solution:
        File(folder_name + "/" + "u.pvd") << u
        File(folder_name + "/" + "p.pvd") << p
        
    # Save for the iterations        
    
    u.rename("u", "velocity")
    p.rename("p", "pressure")
    
    File("./output/u" + str(c) + ".pvd") << u
    File("./output/p" + str(c) + ".pvd") << p    
    
    print "END of N-S"
    
    ### INTERIM DEFINING NEW MEASURES ###
    #sub_bottom = CellFunction("size_t", bottom_mesh)
    #sub_bottom = MeshFunction("size_t", bottom_mesh)
    #sub_bottom.set_all(0)
    #dx_b = Measure('dx')[sub_bottom]
    
    ### CONVECTION-DIFFUSION IN CHILD MESH OR SUBMESH###
    
    # Parameters
    Pe = Constant(1.0e10)
    #Pe = Constant(5.0e2)
    
    #t_end = 10.0
    t_end = 1.0 
    
    dt = 1.0 # only ONE time step
    
    smesh = bottom_mesh
    #smesh = UnitSquareMesh(5, 5)
    ##print smesh.topology().dim()
    ##print smesh.geometry().dim()
    
    def z0_boundary(x, on_boundary):
        #return on_boundary
        return on_boundary and not near(x[0], 3.0)
    
    # Define function spaces
    #Z = FunctionSpace(smesh, "DG", 1)
    Z = FunctionSpace(smesh, "CG", 1)
    ##print Z
    
    # Vector function space for N-S Velocity
    W = VectorFunctionSpace(smesh, "CG", 1)
    ##print W
    
    #ic = Expression("((pow(x[0]-0.25,2)+pow(x[1]-0.25,2))<0.2*0.2)?(-25*((pow(x[0]-0.25,2)+pow(x[1]-0.25,2))-0.2*0.2)):(0.0)")
    #ic = Expression("((pow(x[0]-0.3,2)+pow(x[1]-0.3,2))<0.5*0.5)?(1.0):(0.0)")
    #ic = Expression(('x[0]'), domain = smesh)
    #ic = Constant(0.005)
    
    class InitialCondition(Expression):
        def eval_cell(self, value, x, ufc_cell):
            if x[0] <= 0.95 and x[0] >= 1.05:
            #if x[0] >= 0.4:    
                value[0] = 0.001
            else:
                value[0] = 0.001
    
    # 2D advection velocity#
    #b = Expression(("-(x[1]-0.5)","(x[0]-0.5)"), domain=smesh)
    #b = Expression(("0.5","0.0"), domain=smesh)
    # 3D advection velocity#
    #b = Expression(("-(x[1]-0.5)","(x[0]-0.5)","0.0"), domain=smesh)
    #b = Expression(("0.5","0.0","0.0"), domain=smesh)
    
    #bc = DirichletBC(Z,Constant(0.0),DomainBoundary())
    bc = DirichletBC(Z, Constant(0.0), z0_boundary)
    	
    # Define unknown and test function(s)
    y = TestFunction(Z)
    z = TrialFunction(Z)
    ##print "about z; grad(z):"
    ##print z.rank()
    ##print z.geometric_dimension()
    ##print grad(z).rank()
    #print grad(z).geometric_dimension()
    ##print div(z).rank()
    #print div(z).geometric_dimension()
    
    z1 = Function(Z)
    z0 = Function(Z)
    #z0 = interpolate(ic, Z)
    z0.interpolate(InitialCondition())
    #print "about z0; grad(z0):"
    #print z0.rank()
    #print z0.geometric_dimension()
    #print grad(z0).rank()
    #print grad(z0).geometric_dimension()
    #print div(z0).rank()
    #print div(z0).geometric_dimension()
    
    #U = b
    U = Function(W)
    ##print U.vector().array()
    ##print "about advection velocity U:"
    ##print U.rank()
    ##print U.geometric_dimension()
    ##print grad(U).rank()
    ##print nabla_grad(U).rank()
    #print grad(U).geometric_dimension()
    ##print div(U).rank()
    ##print nabla_div(U).rank()
    #print div(U).geometric_dimension()
           
    # make a map of dofs from the bottom mesh to the global/parent mesh
    gdim = smesh.geometry().dim()
    
    smsh_dof_coordinates = W.dofmap().tabulate_all_coordinates(smesh).reshape(-1, gdim)
    ##print smsh_dof_coordinates
    
    mesh_dof_coordinates = V.dofmap().tabulate_all_coordinates(mesh).reshape(-1, gdim)
    ##print mesh_dof_coordinates
    
    sub_to_glob_map = {}
    for sub_dof_nr, sub_dof_coords in enumerate(smsh_dof_coordinates):
        corresponding_dofs = [i for i, coords in enumerate(mesh_dof_coordinates) if np.array_equal(coords, sub_dof_coords)]
        
        #print corresponding_dofs
        
        if len(corresponding_dofs) == 2:
            if sub_dof_nr % 2 == 0:
                sub_to_glob_map[sub_dof_nr] = corresponding_dofs[0]
            else:
                sub_to_glob_map[sub_dof_nr] = corresponding_dofs[1]
        else:
            raise NameError("Degrees of freedom not matching.")
    ##print sub_to_glob_map
    
    for sub_dof, glob_dof in sub_to_glob_map.iteritems():    
        U.vector()[sub_dof] = nodal_values_vel[glob_dof]
        
    ##print U.vector().array()
    
    ### SEDIMENT TRANSPORT FLUX FORMULAE ###
    
    #a = 1E-4
    a = 0.0002391
    
    #q = a * U * U**2 # simple Grass model
    q = a * U        # simple linear model
    
    #U.vector()[:] = np.array([ 1.0, 0.0, 1.0, 0.0, 1.0, 0.0])
    #U.vector()[:] = np.array([ 1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,
    #  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,
    #  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,
    #  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,  0.])
    #print U.rank()
    #print U.geometric_dimension()
    
    # Stabilization
    h = CellSize(smesh)
    #print h
    n = FacetNormal(smesh)
    #print n
    theta = Constant(0.5)
    alpha = Constant(0.1)
    #f = Constant(0.0)
    
    # Define variational forms
    
    ### CG INTERIOR PENALTY ###
    
    #a0 = (1.0/Pe)*inner(grad(z0), grad(y))*dx + inner(q, grad(z0))*y*dx
    a0 = div(q)*y*dx
    
    #a1 = (1.0/Pe)*inner(grad(z), grad(y))*dx + inner(q, grad(z))*y*dx
    a1 = div(q)*y*dx
    
    r = avg(alpha)*avg(h)**2*inner(jump(grad(z),n), jump(grad(y),n))*dS
    
    # other alternatives for stabilisation:
    #r = avg(alpha)*avg(h)**2*inner(jump(grad(div(q)),n), jump(grad(y),n))*dS
    #r = avg(alpha)*avg(h)**2*inner(jump(div(q)), jump(grad(y),n))*dS # not good
    
    ### DG ###
    '''
    # ( dot(v, n) + |dot(v, n)| )/2.0
    qn = (dot(q, n) + abs(dot(q, n)))/2.0
    
    def a(u,v) :
            # Bilinear form
            a_int = dot(grad(v), (1.0/Pe)*grad(u) - q*u)*dx
            
            a_fac = (1.0/Pe)*(alpha/avg(h))*dot(jump(u, n), jump(v, n))*dS \
                    - (1.0/Pe)*dot(avg(grad(u)), jump(v, n))*dS \
                    - (1.0/Pe)*dot(jump(u, n), avg(grad(v)))*dS
            
            a_vel = dot(jump(v), qn('+')*u('+') - qn('-')*u('-') )*dS  + dot(v, qn*u)*ds
            
            a = a_int + a_fac + a_vel
            return a
    
    a0 = a(z0, y)
    
    a1 = a(z, y)
    '''
    
    A = inner((z - z0)/dt, y)*dx + theta*a1 + (1-theta)*a0
    
    #F = A
    F = A + r
    
    # Create files for storing results
    #file = File("results_%s/u.pvd" % (fileName))
    file = File ( "advection/z.pvd" )
    
    ffc_options = {"optimize": True, "quadrature_degree": 8}
    problem = LinearVariationalProblem(lhs(F), rhs(F), z1, [bc], form_compiler_parameters=ffc_options)
    solver = LinearVariationalSolver(problem)
    
    #J = derivative(F, z1)
    
    #problem = NonlinearVariationalProblem(F, z1, [bc], J, form_compiler_parameters=ffc_options)
    #solver = NonlinearVariationalSolver(problem)
    
    # Time-stepping
    t = 0.0
    
    z1.assign(z0)
    
    while t < t_end:
    
    	print "t =", t, "end t=", t_end
    
    	# Compute
    	solver.solve()
    	##plot(z1, interactive = False)
    	# Save to file
    	file << z1    
    	# Move to next time step
    	z0.assign(z1)
    	t += dt
    
    # Hold plot
    ##interactive()
    
    # Save for the iterations        
    z0.rename("z", "height")
    File("./output/z" + str(c) + ".pvd") << z0
    
    print "*" * 25
    
    ##scoordinates = smesh.coordinates()
    ##print scoordinates
    ##print len(scoordinates)
    ##print scoordinates[0]
    ##print scoordinates[0][1]
    ##print z0(scoordinates[0])
    
    # make a map of dofs from the bottom mesh to the boundarymesh
    Zsmsh_dof_coordinates = Z.dofmap().tabulate_all_coordinates(smesh).reshape(-1, gdim)
    ##print Zsmsh_dof_coordinates
    
    sub_to_bmesh_map = {}
    for sub_dof_nr, sub_dof_coords in enumerate(Zsmsh_dof_coordinates):
        corresponding_dofs = [i for i, coords in enumerate(boundary.coordinates()) if np.array_equal(coords, sub_dof_coords)]
        
        #print corresponding_dofs
        
        if len(corresponding_dofs) == 1:
            sub_to_bmesh_map[sub_dof_nr] = corresponding_dofs[0]
        else:
            raise NameError("Degrees of freedom not matching.")
    ##print sub_to_bmesh_map
    
    nodal_values_z = z0.vector().array()
    ##print nodal_values_z
    ##print len(nodal_values_z)
    
    ##vertex_values = z0.compute_vertex_values()
    ##print vertex_values
    ##print len(vertex_values)
    
    ##for i, x in enumerate(scoordinates):        
    ##    print("vertex %d: \t u(%s) = %s" % (i, x, z0(x)))    
    
    ##v2d = vertex_to_dof_map(Z)
    ##print v2d
    ##print nodal_values_z[v2d]
    
    ### move boundarymesh coordinates according to advection solution
    
    ##print boundary.coordinates()
    
    for sub_dof, bnd_dof in sub_to_bmesh_map.iteritems():    
        boundary.coordinates()[bnd_dof][1] += nodal_values_z[sub_dof] * -1 ## vertical displacement of vertices
    
    print "*" * 25
    ##print boundary.coordinates()
    
    ### ALE ###
    
    #plot(mesh, interactive = True)
    #plot(boundary, interactive = True)
    #plot(smesh, interactive = True)
    
    '''
    bcoordinates = boundary.coordinates()
    
    for x in bcoordinates:
        print x
        if x[1] <= DOLFIN_EPS:
            x[0] += 0.0            
            x[1] -= 5.0
    
    print bcoordinates
    
    for x in scoordinates:
        print x
        x[0] += 0.0            
        x[1] -= 5.0
     
    print scoordinates
    '''
       
    # Move mesh
    #disp = Expression(("0.0", "-5*x[1]"), degree=1)
    #ALE.move(mesh, disp)
    
    ALE.move(mesh, boundary) 
    
    #plot(mesh, interactive = True)
    
    mesh.smooth()
    
    ##plot(mesh, interactive = False)
    
    #interactive()
    
    #File("pipe_mesh.pvd") << mesh 
    
    File("../mesh_pipe/pipe_coarse.xml") << mesh
