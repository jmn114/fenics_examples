"""
2D sloshing in rectangular tank
===============================

- Potential theory
- Linear free surface conditions
- 2D geometry
- Pure horizontal motion
- Rectangular tank
"""
from __future__ import division
from math import pi, tanh
import time
import numpy
import matplotlib.pyplot as plt
import dolfin as df
from scipy.integrate import ode

#############################################################
# Input definition

wall_motion_amplitude = 0.1  # Max amplitude of the horizontal tank wall motion 
wall_motion_period = 5.3     # seconds
dt = wall_motion_period/30   # we use a fixed number of time steps per oscilation
tmax = wall_motion_period*10  # How long will our simulation last
tramp = wall_motion_period*3 # How long will we spend ramping up the motion from zero

Nelemx = 40
Nelemz = 20

# Dimension in meters
width = 20
still_water_height = 10
GEOM_TOL = width/1e6

# Acceleration of gravity
g = 9.81

# Order of the elements
order = 1

# Make snapshots for animation
animate = False

# Print expected natural periods
print 'Natural periods for a rectangular tank (Faltinsen):'
print 'Mode Frequency    Period'
print ' [-]   [rad/s]       [s]'
for n in range(1, 5):
    wn = (g * n*pi/width * tanh(n*pi/width * still_water_height))**0.5
    print '%4d %9.6f %9.6f' % (n, wn, 2*pi/wn)
print

#############################################################
# Create mesh and define function space

mesh = df.RectangleMesh(df.Point(0, 0), df.Point(width, still_water_height), Nelemx, Nelemz, 'left/right')
coords = mesh.coordinates()

V = df.FunctionSpace(mesh, 'Lagrange', order)

#############################################################
# Get the coordinates of the upper vertices and define the
# free surface functions that will be integrated in time

# Indices of the top vertices
idx_top = numpy.where(coords[:,1] > still_water_height - GEOM_TOL)[0]

# Sort indices by x coordinate of the corresponding vertex
idx_top = list(idx_top)
idx_top.sort(key=lambda i: coords[i,0])
idx_top = numpy.array(idx_top, int)

# The corresponding x-coordinates
xcoords = coords[idx_top,0]

#############################################################
# Define boundary conditions

def define_boundary(defining_func, boundary_parts, boundary_id):
    """
    Define a boundary
    
    - defining_func is used instead of the normal SubDomain eval method: defining_func(x, on_boundary)
    - boundary_parts is a MeshFunction used to distinguish parts of the boundary from each other
    - boundary_id is the id of the new boundary in the boundary_parts dictionary
    """
    class Boundary(df.SubDomain):
        def inside(self, x, on_boundary):
            return defining_func(x, on_boundary)
    boundary = Boundary()
    boundary.mark(boundary_parts, boundary_id)
    return boundary

# Each boundary part gets its own ID
ID_LEFT_SIDE, ID_RIGHT_SIDE, ID_BOTTOM, ID_TOP, ID_OTHER = 1, 2, 3, 4, 5
boundary_parts = df.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundary_parts.set_all(ID_OTHER)

# Define the boundaries (tank side walls, bottom and the free surface)
define_boundary(lambda x, on_boundary: on_boundary and x[0] < GEOM_TOL, boundary_parts, ID_LEFT_SIDE)
define_boundary(lambda x, on_boundary: on_boundary and x[0] > width-GEOM_TOL, boundary_parts, ID_RIGHT_SIDE)
define_boundary(lambda x, on_boundary: on_boundary and x[1] < GEOM_TOL, boundary_parts, ID_BOTTOM)
define_boundary(lambda x, on_boundary: on_boundary and x[1] > still_water_height - GEOM_TOL, boundary_parts, ID_TOP)

class DynamicFreeSurfaceCondition(df.Expression):
    """
    This function will apply the value of the velocity potential 
    calculated in the previous time step to the free surface
    
    The phi_on_surface attribute will be set by the time stepping
    code and changed for each time step
    """
    def eval(self, values, x):
        values[0] = numpy.interp(x[0], xp=xcoords, fp=self.phi_on_surface)
bc_exp_top = DynamicFreeSurfaceCondition()
bc_top = df.DirichletBC(V, bc_exp_top, boundary_parts, ID_TOP)

#############################################################
# Define variational problems

# The velocity potential
u = df.TrialFunction(V)
v = df.TestFunction(V)
a = df.inner(df.nabla_grad(u), df.nabla_grad(v))*df.dx
phi = df.Function(V)

# The z-derivative or the velocity potential
u2 = df.TrialFunction(V)
v2 = df.TestFunction(V)
a2 = u2*v2*df.dx
L2 = df.Dx(phi, 1)*v2*df.dx
phi_z = df.Function(V)

#############################################################
# Compute solution

# Pre-compute stiffness  matrix A, the left hand side of the a=L equation system
#A = df.assemble(a, exterior_facet_domains=boundary_parts)
A = df.assemble(a)

# The angular velocity of the motion
omega = 2*pi/wall_motion_period

# The time vector
t_vec = numpy.arange(0, tmax+dt/10, dt)

#print len(t_vec)

# Somewhere to store the solution after each time step
etas = numpy.zeros((len(t_vec), len(idx_top)), float)

#print numpy.size(etas,0)
#print numpy.size(etas,1)

# Create ramp function for the tank motion and velocity
ramp = 0.5 - 0.5*numpy.cos(t_vec/tramp*pi)
ramp[t_vec > tramp] = 1.0

# Pre compute the velocity and motion at each time step
motion = wall_motion_amplitude * numpy.cos(omega * t_vec) * ramp
velocity = - omega * wall_motion_amplitude * numpy.sin(omega * t_vec) * ramp

def compute_vertical_velocity_on_free_surface(wall_velocity, phi_on_surface):
    """
    Compute the vertical velocity on the surface from the velocity potential
    The input is the velocity potential on the free surface to use in Dirichlet
    boundary condition and the horizontal wall velocity for Neumann conditions
    """
    
    #ds = df.Measure("ds", domain=mesh, subdomain_data=boundary_parts)
    
    # Calculate right hand side of the a=L equation system
    L = (-df.Constant(-wall_velocity)*v*df.ds(ID_LEFT_SIDE, subdomain_data=boundary_parts)
         -df.Constant(wall_velocity)*v*df.ds(ID_RIGHT_SIDE, subdomain_data=boundary_parts))
    
    # Update the top boundary condition
    bc_exp_top.phi_on_surface = phi_on_surface
    
    # Assemble right hand side and apply boundary conditions
    #b = df.assemble(L, exterior_facet_domains=boundary_parts)
    b = df.assemble(L)
    bc_top.apply(A, b)     
    
    # Solve for the velocity potential
    df.solve(A, phi.vector(), b)
    
    #print phi.vector().array()
    
    #boundary_conditions = [bc_top]
    #df.solve(a == L, phi, boundary_conditions)    
    
    # Solve for the z-derivative of the velocity potential
    A2 = df.assemble(a2)
    b2 = df.assemble(L2)
    df.solve(A2, phi_z.vector(), b2)
    
    # Return the vertical velocity at the free surface
    dphidz_vert = phi_z.compute_vertex_values()
    #print dphidz_vert
    phi_vert = phi.compute_vertex_values()
    #print phi_vert
    return dphidz_vert, phi_vert

def ode_deriv(t, y):
    """
    We define an ode for N*2 equations where the first N equations
    are d/dt(eta) and the last N equations are d/dz(phi) on the
    free surface
    """
    global phi_vert # this is global just to allow plotting it in the end
    
    # Height of free surface (displacement from mean position)
    eta = y[:N]
    
    # The value of the velocity potential at the surface    
    phi_surf = y[N:]
    
    # Calculate the vertical velocity on the free surface from the velocity potential
    wall_velocity = numpy.interp(x=t, xp=t_vec, fp=velocity)
    #print wall_velocity
    
    dphidz_vert, phi_vert = compute_vertical_velocity_on_free_surface(wall_velocity, phi_surf)
    #print dphidz_vert
    #print phi_vert
    
    dphidz_surface = dphidz_vert[idx_top]    
        
    # Assign the time derivatives of eta and phi on the surface
    dydt = numpy.zeros(N*2, float)    
    dydt[:N] = dphidz_surface
    dydt[N:] = -g*eta
    
    return dydt

# Define initial conditions
# The first N are the starting surface elevations,
# The last N are the stsrting velocity potential on the free surface
N = len(idx_top)
y0 = numpy.zeros(N*2, float)

# Use scipy to integrate the boundary conditions at the free surface 
ode_integrator = ode(ode_deriv)
ode_integrator.set_initial_value(y0)
ode_integrator.set_integrator('dop853')

# Loop through each time step

t1 = time.time()
for it, t in enumerate(t_vec):
    if it == 0:
        continue
    
    # Integrate to the next time step
    y = ode_integrator.integrate(t)    
            
    #print len(y)
    
    assert ode_integrator.successful()
    
    # Store the solution
    etas[it] = y[:N]
    
    print '%10.5f %15.2e %15.2e %7.2f' % (t, etas[it,0], etas[it,-1], time.time()-t1)
    

#while ode_integrator.successful() and ode_integrator.t < t_vec[300]:
#    print(ode_integrator.t+dt, ode_integrator.integrate(ode_integrator.t+dt))

#df.plot(mesh)

plt.figure('velocity potential, phi', figsize=(6,4))
X = mesh.coordinates()[:,0].reshape((Nelemz+1, Nelemx+1)).T
Y = mesh.coordinates()[:,1].reshape((Nelemz+1, Nelemx+1)).T
Z = phi_vert.reshape((Nelemz+1, Nelemx+1)).T
plt.contourf(X, Y, Z)
plt.tight_layout()

plt.figure('Eta at ends', figsize=(6,4))
plt.plot
plt.plot(t_vec, etas[:,0], label='End 1')
plt.plot(t_vec, etas[:,-1], label='End 2')
plt.plot(t_vec, motion[:], label='Wall motion')
plt.legend(loc='upper left')
plt.tight_layout()

#####################################################################
# Animation
if animate:
    metadata = dict(title='Free surface', artist='Matplotlib')
    fps = int(round(1/dt))
    print 'fps', fps
    
    fig = plt.figure()
    ax = plt.subplot(111)
    
    xcoords = coords[idx_top][:,0]
    maxampl = abs(etas).max().max()
    for i, (t, eta) in enumerate(zip(t_vec, etas)):
        ax.plot(xcoords, eta + still_water_height)
        plt.xlim(0, width)
        plt.ylim(still_water_height-maxampl*1.2, still_water_height+maxampl*1.2)
        plt.title('t  = %5.2f' % t)
        #fig.savefig('fig/fs%07d.png' % i)
        ax.clear()
    plt.close(fig)

last = t_vec > tmax-wall_motion_period*1.5
print etas[last,0].max(), etas[last,0].min()
print etas[last,-1].max(), etas[last,-1].min()
 
plt.show()