from dolfin import *



mesh = UnitSquareMesh(10, 10)
boundary = BoundaryMesh(mesh, 'exterior')
#plot (mesh)
#plot (boundary)
#interactive()

disp = Expression(("0.3*x[0]*x[1]", "0.5*(1.0-x[1])"), degree=2)
ALE.move(boundary, disp)

ALE.move(mesh, boundary)

boundary_new = BoundaryMesh(mesh, 'exterior')
assert boundary.topology().hash() == boundary_new.topology().hash()

plot (mesh)
interactive()

'''

mesh = UnitSquareMesh(4, 5)

plot(mesh)
interactive()

# Make some cell function
# FIXME: Initialization by array indexing is probably
#        not a good way for parallel test
cellfunc = CellFunction('size_t', mesh)
cellfunc.array()[0:24] = 0
cellfunc.array()[24:]  = 1

# Create submeshes - this does not work in parallel
submesh0 = SubMesh(mesh, cellfunc, 0)
submesh1 = SubMesh(mesh, cellfunc, 1)

plot(submesh0)
plot(submesh1)
interactive()

# Move submesh0
disp = Constant(("0.1", "-0.1"))
ALE.move(submesh0, disp)

# Move and smooth submesh1 accordignly
ALE.move(submesh1, submesh0)

# Move mesh accordingly
parent_vertex_indices_0 = submesh0.data().array('parent_vertex_indices', 0)
parent_vertex_indices_1 = submesh1.data().array('parent_vertex_indices', 0)
mesh.coordinates()[parent_vertex_indices_0[:]] = submesh0.coordinates()[:]
mesh.coordinates()[parent_vertex_indices_1[:]] = submesh1.coordinates()[:]

plot (mesh)
interactive()

'''