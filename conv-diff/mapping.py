from dolfin import *
import numpy as np

degree = 1
n = 2

mesh = UnitSquareMesh(n, n)
plot(mesh, interactive = True)
print mesh
print mesh.coordinates()

V = FunctionSpace(mesh, 'CG', degree)
class BottomBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0)
bottomboundary = BottomBoundary()
boundaries = FacetFunction('size_t', mesh)
boundaries.set_all(0)
bottomboundary.mark(boundaries, 1)
ds = Measure('ds')[boundaries]

boundarymesh = BoundaryMesh(mesh, 'exterior')
plot(boundarymesh, interactive = True)
print boundarymesh
print boundarymesh.coordinates()

Vb = FunctionSpace(boundarymesh, 'CG', degree)

bottom_mesh = SubMesh(boundarymesh, bottomboundary)
plot(bottom_mesh, interactive = True)
L = FunctionSpace(bottom_mesh, 'CG', degree)
bf = Function(L)
# here goes some more code to find the correct bf
# but here a dummy function will do
bf.vector()[:] = np.linspace(0.0, 10.0, n * degree + 1)

# node ids of bmesh in Bmesh
bottom_mesh_vertex = dof_to_vertex_map(L)
print bottom_mesh_vertex

map_n_sm_bm = bottom_mesh.data().array('parent_vertex_indices', 0)
print map_n_sm_bm
# node ids of Bmesh in mesh
map_n_bm_m = boundarymesh.entity_map(0)
#print map_n_bm_m
# node ids of bmesh in mesh
map_n_sm_m = map_n_bm_m.array()[map_n_sm_bm]
print map_n_sm_m

print "*" * 25

print vertex_to_dof_map(Vb)
print dof_to_vertex_map(Vb)
print boundarymesh.coordinates()[dof_to_vertex_map(Vb)]

print "*" * 25

print vertex_to_dof_map(L)
print dof_to_vertex_map(L)
print bottom_mesh.coordinates()[dof_to_vertex_map(L)]

print "*" * 25

print vertex_to_dof_map(V)
print dof_to_vertex_map(V)
print mesh.coordinates()[dof_to_vertex_map(V)]

#bottom_to_mesh = compute_vertex_map(bottom_mesh, mesh)
#print(bottom_to_mesh)

print "*" * 25

# make a map of dofs from the submesh to boundarymesh and the original mesh
gdim = bottom_mesh.geometry().dim()
#print gdim
#print bottom_mesh.topology().dim()

bmesh_dof_coordinates = Vb.dofmap().tabulate_all_coordinates(boundarymesh).reshape(-1, gdim)
print bmesh_dof_coordinates

print "*" * 25

smsh_dof_coordinates = L.dofmap().tabulate_all_coordinates(bottom_mesh).reshape(-1, gdim)
print smsh_dof_coordinates

print "*" * 25

mesh_dof_coordinates = V.dofmap().tabulate_all_coordinates(mesh).reshape(-1, gdim)
print mesh_dof_coordinates

### MAPPINGS ###

sub_to_glob_map = {}
for sub_dof_nr, sub_dof_coords in enumerate(smsh_dof_coordinates):
    corresponding_dofs = [i for i, coords in enumerate(mesh_dof_coordinates) if np.array_equal(coords, sub_dof_coords)]
    
    #print corresponding_dofs
    
    if len(corresponding_dofs) == 1:
        sub_to_glob_map[sub_dof_nr] = corresponding_dofs[0]
    else:
        raise NameError("Degrees of freedom not matching.")
print sub_to_glob_map

sub_to_bnd_map = {}
for sub_dof_nr, sub_dof_coords in enumerate(smsh_dof_coordinates):
    corresponding_dofs = [i for i, coords in enumerate(bmesh_dof_coordinates) if np.array_equal(coords, sub_dof_coords)]
    
    #print corresponding_dofs
    
    if len(corresponding_dofs) == 1:
        sub_to_bnd_map[sub_dof_nr] = corresponding_dofs[0]
    else:
        raise NameError("Degrees of freedom not matching.")
print sub_to_bnd_map

bmesh_to_bnd_map = {}
for bmesh_dof_nr, bmesh_dof_coords in enumerate(boundarymesh.coordinates()):
    corresponding_dofs = [i for i, coords in enumerate(bmesh_dof_coordinates) if np.array_equal(coords, bmesh_dof_coords)]
    
    #print corresponding_dofs
    
    if len(corresponding_dofs) == 1:
        bmesh_to_bnd_map[bmesh_dof_nr] = corresponding_dofs[0]
    else:
        raise NameError("Degrees of freedom not matching.")
print bmesh_to_bnd_map

sub_to_bmesh_map = {}
for sub_dof_nr, sub_dof_coords in enumerate(smsh_dof_coordinates):
    corresponding_dofs = [i for i, coords in enumerate(boundarymesh.coordinates()) if np.array_equal(coords, sub_dof_coords)]
    
    #print corresponding_dofs
    
    if len(corresponding_dofs) == 1:
        sub_to_bmesh_map[sub_dof_nr] = corresponding_dofs[0]
    else:
        raise NameError("Degrees of freedom not matching.")
print sub_to_bmesh_map

print "*" * 25

'''
un = Function(V)
un.vector().array().fill(0.0)
for Vs_dof, val in enumerate(us.vector().array()):
    submesh_vertex = dof_to_vertex_map(L)[Vs_dof]
    boundary_vertex = bottom_mesh.data().array('parent_vertex_indices', 0)[submesh_vertex]
    mesh_vertex = boundarymesh.entity_map(0)[int(boundary_vertex)] # np.uint not accepted
    V_dof = vertex_to_dof_map(V)[mesh_vertex]
    un.vector()[V_dof] = val
'''

print "*" * 25

'''
coordinates = mesh.coordinates()
bcoordinates = boundarymesh.coordinates()

for x in bcoordinates:
    print x
    if x[1] <= DOLFIN_EPS:
        x[0] += 0.0            
        x[1] -= 1.0

print bcoordinates
'''

z = np.array([ 0., -0.05, 0. ])

print bmesh_dof_coordinates
print bmesh_dof_coordinates[3]
print bmesh_dof_coordinates[3][1]
print bmesh_dof_coordinates[3][0]
print boundarymesh.coordinates()
print boundarymesh.coordinates()[3]
print boundarymesh.coordinates()[3][0]
print boundarymesh.coordinates()[3][1]

for sub_dof, bnd_dof in sub_to_bmesh_map.iteritems():
    #bmesh_dof_coordinates[bnd_dof][1] = z[sub_dof]
    boundarymesh.coordinates()[bnd_dof][1] = z[sub_dof]

print boundarymesh.coordinates()

print "***"

ALE.move(mesh, boundarymesh) 

plot(mesh, interactive = True)

mesh.smooth()

plot(mesh, interactive = True)

# warp bf to V
bf_to_V = Function(V)
for sub_dof, glob_dof in sub_to_glob_map.iteritems():
    bf_to_V.vector()[glob_dof] = bf.vector().array()[sub_dof]
plot(bf_to_V, interactive=True)