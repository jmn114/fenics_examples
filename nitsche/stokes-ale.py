#! /usr/bin/env python
"""
Demo for Nitsche-type free-slip boundary conditions
"""

from dolfin import *
import numpy as np
#from vtktools import *
#vtk_writer = VTK_XML_Serial_Unstructured()

np.set_printoptions(threshold=np.inf)
#parameters['allow_extrapolation'] = True
#from __future__ import print_function
set_log_level(ERROR) # avoid outputs like 'solving variational problem' by FEniCS

#f = open("./log.txt", "w")

def diff_fwd(x, y): 
    return np.diff(y)/np.diff(x)

def diff_central(x, y):
    x0 = x[:-2]
    x1 = x[1:-1]
    x2 = x[2:]
    y0 = y[:-2]
    y1 = y[1:-1]
    y2 = y[2:]
    f = (np.array(x2) - np.array(x1))/(np.array(x2) - np.array(x0))
    return (1-f)*(np.array(y2) - np.array(y1))/(np.array(x2) - np.array(x1)) + f*(np.array(y1) - np.array(y0))/(np.array(x1) - np.array(x0))

#mesh = refine(Mesh("cylinderbump.xml"))

def ale_stokes(c):

    dt = 0.00001
    
    #mesh = Mesh("./pipe_coarse.xml")
    mesh = Mesh("./ale_mesh.xml")
    #plot(mesh, interactive=True)               
        
    ### STOKES SECTION ###
    
    V = VectorFunctionSpace(mesh, 'CG', 1)
    Q = FunctionSpace(mesh, 'CG', 1)
    
    mu = Constant(1.0)
    F = Constant((0.0,0.0))
    
    #W = V*Q
    W = MixedFunctionSpace([V, Q])    
    
    # variational problem
    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)
    
    # boundary conditions
    #leftexp = Expression(("x[1]*(3-x[1])/10", "0.0"))
    #left = DirichletBC(W.sub(0), leftexp, "near(x[0], -5)")
    #top = DirichletBC(W.sub(0), (0,0), "near(x[1], 3)")
    
    #leftexp = Expression(("x[1]*(0.4-x[1])/10", "0.0"))
    leftexp = Expression(("0.5 * 4.0 * x[1] * (0.4 - x[1]) / ( 0.4 * 0.4 )", "0.0"))
    left = DirichletBC(W.sub(0), leftexp, "near(x[0], 0.0)")
    top = DirichletBC(W.sub(0), (0,0), "near(x[1], 0.4)")
    bcs = [left, top]
    
    class Bottom(SubDomain):
      def inside(self, x, on_boundary):
        #return on_boundary and not (near(x[0], -5) or near(x[0], 5) or near(x[1], 3))
        return on_boundary and not (near(x[0], 0.0) or near(x[0], 3.0) or near(x[1], 3))
    bottom = Bottom()    
        
    bounds = MeshFunction('size_t', mesh, "../pipe/pipe_coarse_facet_region.xml")
    boundary_ids = bounds.array()    
    facet_nodes = np.array([f.entities(0) for f in facets(mesh)])
    facet_indices = [f.index() for f in facets(mesh)]
    
    # create list of nodes on pointwise boundary
    ptws_bndy_nds = []
    for j in range(len(facet_nodes)):
      if boundary_ids[j] == 3:
        ptws_bndy_nds.extend(list(facet_nodes[j]))
    ptws_bndy_nds = sorted(set(ptws_bndy_nds))    
    print len(ptws_bndy_nds)
    #print ptws_bndy_nds
    #print min(ptws_bndy_nds)
    #print max(ptws_bndy_nds)    
    
    boundaries = FacetFunction("size_t", mesh)
    #boundaries = FacetFunction("uint", mesh)   
    
    boundaries.set_all(0)
    boundaries.set_all(0)
    bottom.mark(boundaries, 1)
    ds = Measure("ds")[boundaries]
    
    alpha = Constant(1./10)
    beta = Constant(10)
    h = CellSize(mesh)
    n = FacetNormal(mesh)
    
    # (bi)linear forms
    def a(u,v): return inner(mu*grad(u),grad(v))*dx
    def b(v,q): return - div(v)*q*dx
    def f(v):   return dot(F, v)*dx
    def t(u,p): return dot(2*mu*sym(grad(u)),n) - p*n
    
    stokes = a(u,v) + b(v,p) + b(u,q) - f(v) \
           - dot(n,t(u,p))*dot(v,n)*ds(1) - dot(u,n)*dot(n,t(v,q))*ds(1) \
           + beta/h*dot(u,n)*dot(v,n)*ds(1) \
           + alpha*h**2*dot(F - grad(p), grad(q))*dx
    
    # solve variational problem
    wh = Function(W)
    #print 'size:', wh.vector().size()
    #print mesh.num_vertices()
    #print mesh.num_cells()
    
    solve(lhs(stokes) == rhs(stokes), wh, bcs)
    #solve(lhs(stokes) == rhs(stokes), wh, bcs, solver_parameters = { "linear_solver": "tfqmr", "preconditioner": "amg"} )
    
    uh, ph = wh.split()
    #plot(uh, interactive = True)
    #plot(ph, interactive = True)
    File("velocity.pvd") << uh
    File("pressure.pvd") << ph
    
    uhc = Function(uh)
    u_nodal_values = uhc.vector()
    u_array = u_nodal_values.array()
    #print len(u_array)
    #print u_array
    
    #print uhc(0.0,0.0)[0]
    #print uhc(0.00993792544783,0.0)[0]
    #print uhc(0.00993792544783,-0.000755005170072)[0]
    #print uhc(3.0,0.0)[0]
    
    coor = mesh.coordinates()
    #print len(coor)
    #print coor
    
    #v_d = W.dofmap().vertex_to_dof_map(mesh)
    #d_v = W.dofmap().dof_to_vertex_map(mesh)
    
    #dofs_at_vertices = uhc[d_v]
    
    # let's iterate over vertices
    #for v in vertices(mesh):
    #    print 'vertex index', v.index()
    #    print 'at point', v.point().str()
    #    print 'at coordinates', coor[v.index()]
    #    print 'dof', dofs_at_vertices[v.index()]
    
    d = mesh.geometry().dim()
    #print d
    
    # Coordinates of all dofs in the mixed space
    Wdofs_x = W.dofmap().tabulate_all_coordinates(mesh).reshape((-1, d))
    # Dofs of second subspace of the mixed space
    Q_dofs = W.sub(1).dofmap().dofs()
    # Coordinates of dofs of second subspace of the mixed space
    Q_dofs_x = Wdofs_x[Q_dofs, :]
    di_dx = Q_dofs_x.tolist()
    #print len(di_dx)
    #print di_dx
    #print di_dx[4]
    #print di_dx[7]
    
    #if mesh.num_vertices() == len(Q_dofs_x):
    #    for i in range(mesh.num_vertices()):
    #        if Q_dofs_x[i,1] == 0.0:
    #            print uhc(Q_dofs_x[i,0],0.0)
    #            print 'u(%8g,%8g) = %g' % (Q_dofs_x[i][0], Q_dofs_x[i][1], uhc(Q_dofs_x[i,0],Q_dofs_x[i,1])[0]) 
    
    # Coordinates of all vertices
    vertex_x = coor.reshape((-1, d))
    # vi_vx[vertex_index] = coordinates of index
    vi_vx = vertex_x.tolist()
    #print len(vi_vx)
    #print vi_vx    
    #print vi_vx[4]
    #print vi_vx[5]
    #print vi_vx[6]
    #print vi_vx[7]
    #print vi_vx[99]
    
    #print di_dx.index(vi_vx[4])
    #print di_dx.index(vi_vx[5])
    #print di_dx.index(vi_vx[6])
    #print di_dx.index(vi_vx[7])
    #print di_dx.index(vi_vx[99])
    
    G = []
    if mesh.num_vertices() == len(vi_vx):
        for j in ptws_bndy_nds:
            G.append(di_dx.index(vi_vx[j]))
    
    #print G
    
    X = []
    if mesh.num_vertices() == len(vi_vx):
        for j in G:
            #X.append(vi_vx[j][0])
            X.append(di_dx[j][0])
                
    print len(X)
    #print X
    
    Y = []
    if mesh.num_vertices() == len(vi_vx):
        for j in G:
            #Y.append(vi_vx[j][1])
            Y.append(di_dx[j][1])
    
    print len(Y)
    #print Y
    
    X_sort = sorted(X, key=float)
    print len(X_sort)
    #print X_sort
    
    #print X_sort[1]
    #print X.index(X_sort[1])
    #print Y[X.index(X_sort[1])]
                   
    U_x = []
    for i in range(len(X_sort)):
        #print X_sort[i], Y[X.index(X_sort[i])]    
        U_x.append(uhc(X_sort[i],Y[X.index(X_sort[i])])[0])
      
    print len(U_x)
    #print U_x
    
    Q = []
    for i in range(len(U_x)):
        Q.append(0.4*0.5*U_x[i]*abs(U_x[i])**0.5)
    
    print len(Q)
    #print Q
    
    #Grad Q
    grad_Q = diff_central(X_sort, Q)
    grad_Q = np.insert(grad_Q,0,0.0)
    grad_Q = np.insert(grad_Q,len(grad_Q),0.0)
    #print len(grad_Q)
    #print grad_Q
    
    #if mesh.num_vertices() == len(Q_dofs_x):
    #	for i in range(mesh.num_vertices()):        
    #         print 'u(%8g,%8g) = %g' % (Q_dofs_x[i][0], Q_dofs_x[i][1], u_array[i])		
    '''
    dof_coordinates = Q.dofmap().tabulate_all_coordinates(mesh)
    n = Q.dim()
    d = mesh.geometry().dim()
    dof_coordinates.resize((n, d))
    
    x_dofs = Q.sub(0).dofmap().dofs()
    #print x_dofs
    y_dofs = Q.sub(1).dofmap().dofs()
    #print y_dofs
    
    for x_dof, y_dof in zip(x_dofs, y_dofs):
      print dof_coordinates[x_dof], dof_coordinates[y_dof], uh_c[x_dof], uh_c[y_dof]
      #print dof_coordinates[x_dof], dof_coordinates[y_dof], uh_c[x_dof], uh_c[y_dof]
    '''
    
    ### ALE SECTION ###
          
    # Create boundary mesh
    boundary = BoundaryMesh(mesh, "exterior", True)
    print boundary
    #print boundary.num_vertices()
    #print boundary.coordinates()    
    #print len(boundary.coordinates())
    #print >> f, set(boundary.cells().flat)
    
    # Create sub-domain for mesh and mark everything as 0
    sub_domain_fluid = MeshFunction("size_t", mesh, mesh.topology().dim())
    sub_domain_fluid.set_all(0)
    
    #Bottom sub-domain
    class Seabed(SubDomain):
        def inside(self, x, on_boundary):
            #return x[1] == 0.0
            #return x[0] > 0.0 and x[0] < 3.0 and x[1] < 0.005
            return x[1] <= DOLFIN_EPS             
    seabed = Seabed()      
    
    # Create sub-domain for boundary and mark bottom sub-domain as 1
    sub_domain_seabed = MeshFunction("size_t", boundary, boundary.topology().dim())    
    sub_domain_seabed.set_all(0)
    seabed.mark(sub_domain_seabed, 1)
    #bottom.mark(sub_domain_seabed, 1)
    
    # Extract sub meshes
    fluid_mesh = SubMesh(mesh, sub_domain_fluid, 0)
    print fluid_mesh
    #print fluid_mesh.num_vertices()
    
    seabed_mesh = SubMesh(boundary, sub_domain_seabed, 1)
    print seabed_mesh
    #print seabed_mesh.num_vertices()
    #print seabed_mesh.coordinates()
    #print len(seabed_mesh.coordinates())
    
    # Move mesh
    #disp = Expression(("0.0", "-5*x[1]"), degree=1)
    #ALE.move(boundary, disp)    
    
    count = 0    
    for x in sorted(boundary.coordinates(), key=lambda x: x[0]):
    #for x in seabed_mesh.coordinates():        
        #if x[1] <= DOLFIN_EPS:
        if x[0] == X_sort[count] and x[1] == Y[X.index(X_sort[count])]:    
            x[0] *= 1.0            
            x[1] -= grad_Q[count] * dt
            count += 1
            #print x
    
    ALE.move(mesh, boundary)
    #ALE.move(mesh, seabed_mesh)
    
    # Plot mesh
    #plot(mesh, interactive=False)    
    
    # Write mesh to file
    if c == 1 or c % 50 == 0:
        File("./output/deformed_mesh" + str(c) + ".pvd") << mesh
    
    File("ale_mesh.xml") << mesh      
      
    #W = project(mesh, W)
    #vtk_writer.writePVD("deformed_mesh.pvd")    
        
    #exit()