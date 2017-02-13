from dolfin import *
import numpy as np
#from __future__ import print_function

f = open("./log.txt", "w")

#Bottom subdomain
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] == 0.0
        #return x[0] > 0.0 and x[0] < 3.0 and x[1] < 0.005
bottom = Bottom()

# Create mesh
#mesh = dolfin.UnitSquareMesh(4, 4)
mesh = Mesh("../pipe/pipe.xml")
print mesh
print mesh.num_vertices()
#print >> f, set(mesh.cells().flat)

# Create boundary mesh
boundary = BoundaryMesh(mesh, "exterior", True)
print boundary
print boundary.num_vertices()
#print boundary.coordinates()
#print >> f, set(boundary.cells().flat)

# Create subdomain for mesh and mark everything as 0
sub_domain_fluid = MeshFunction("size_t", mesh, mesh.topology().dim())
sub_domain_fluid.set_all(0)

# Create subdomain for boundary and mark bottom subdomain as 1
sub_domain_seabed = MeshFunction("size_t", boundary, boundary.topology().dim())
sub_domain_seabed.set_all(0)
bottom.mark(sub_domain_seabed, 1)

# Extract sub meshes
fluid_mesh = SubMesh(mesh, sub_domain_fluid, 0)
print fluid_mesh
print fluid_mesh.num_vertices()
#print fluid_mesh.coordinates()

seabed_mesh = SubMesh(boundary, sub_domain_seabed, 1)
print seabed_mesh
print seabed_mesh.num_vertices()
#print seabed_mesh.coordinates()

print "\n"

"""
boundaries = dolfin.MeshFunction('size_t', mesh, "../pipe/pipe_facet_region.xml")
print boundaries
#print >> f, set(boundaries.array().flat)

boundary_ids = boundaries.array()
facet_nodes = np.array([f.entities(0) for f in dolfin.facets(mesh)])

# create list of nodes on pointwise boundary
ptws_bndy_nds = []
for j in range(len(facet_nodes)):
  if boundary_ids[j] == 3:
    ptws_bndy_nds.extend(list(facet_nodes[j]))
ptws_bndy_nds = sorted(set(ptws_bndy_nds))
print len(ptws_bndy_nds)
"""
#print seabed_mesh.data().array('parent_vertex_indices', 0)
#print seabed_mesh.data().array('parent_cell_indices', 1)
#mesh.data().create_array("parent_vertex_indices", len(bottom.coordinates()))
#boundary.data().create_array("parent_vertex_indices", len(bottom.coordinates()))
#print mesh.data()
#print boundary.data()

# Move vertices in boundary
for x in boundary.coordinates():
	if x[1] == 0.0:
    		x[0] *= 1.0    
    		x[1] += -0.05*sin(1.0*x[0])    

# Move mesh
ALE.move(mesh, boundary)
#fluid_mesh.smooth()

# Plot mesh
plot(mesh, interactive=True)
plot(fluid_mesh, title="Fluid")
plot(seabed_mesh, title="Seabed")
interactive()

# Write mesh to file
fid = File("deformed_mesh.pvd") << mesh
#for i in range(0, 20):
#    u[i].rename("u[i]", "u[i]")
#    fid << u[i], i