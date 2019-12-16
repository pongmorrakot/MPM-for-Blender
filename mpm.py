import bpy
import copy
import math
import numpy as np

# Grid resolution (cells)
n = 80;

dt = 1e-4;
frame_dt = 1e-3;
dx = 1.0 / n;
inv_dx = 1.0 / dx;

particle_mass = 1.0
vol = 1.0       # Particle Volume
hardening = 10.0; # Snow hardening factor
E = 1e4;          # Young's Modulus
nu = 0.2;         # Poisson ratio
plastic = True;

mu_0 = E / (2 * (1 + nu));
lambda_0 = E * nu / ((1+nu) * (1 - 2 * nu));

cur_frame = 0
particles = None
py_particles = []
totalParticles = 0
objects = []
gravity = None
domain = None
grid = []
young = None
compression = None
stretch = None
time = None
step = None

class particle:
	x = np.array([0.,0.,0.]) # location
	v = np.array([0.,0.,0.]) # velocity
	F = np.matrix([1.,1.,1.], # deformation Gradient
				 [1.,1.,1.],
				 [1.,1.,1.])
	C = np.matrix([0.,0.,0.], # Affine momentum from APIC
				 [0.,0.,0.],
				 [0.,0.,0.])
	Jp = 1.
	
	def __init__(self, x,v):
		self.x = x
		self.v = v

def parse(b_objects, b_domain, start_frame, b_gridRes, b_force, b_young, b_compression, b_stretch, b_hardening, b_time, b_step):
	# parse particle into particle array
	# parse other object to be objects
	for obj in b_objects:
		try:
			degp = bpy.context.evaluated_depsgraph_get()
			# Emitter Object
			cur = obj
			# Evaluate the depsgraph (Important step)
			particle_systems = cur.evaluated_get(degp).particle_systems
			# All particles of first particle-system which has index "0"	
			particles = particle_systems[0].particles		
		except:
			# this object doesn't have a particle system
			objects.append(obj)
		else:
			# this object has a particle system
			totalParticles += len(particles)
			locList = [0]*(3*totalParticles)
			vList = [0]*(3*totalParticles)
			particles.foreach_get("location", locList)
			particles.foreach_get("velocity", vList)
			i = 0
			while i < 3*len(particles):
				x = np.array([locList[i],locList[i+1],locList[i+2]])
				v = np.array([vList[i],vList[i+1],vList[i+2]])
				particles.append(particle(x,v))
				i += 3			
	cur_frame = start_frame
	domain = b_domain
	u,v,w = b_domain.dimensions
	x,y,z = b_domain.location
	n = b_gridRes
	dx = 1.0 / n;
	inv_dx = 1.0 / dx;
	gravity = np.array([b_force[0],b_force[1],b_force[2]])
	E = b_young
	compression = b_compression
	stretch = b_stretch
	hardening = b_hardening
	time = b_time
	step = b_step
	
def main():
	# set up
	setup()
	t = 0
	while t < time:
		cur_frame += 1
		s = 0
		while s < step:
			#do the simulation
			simulate()
			s += 1
		#output the info to keyframe
		draw(cur_frame)
		t += 1
		
def setup():
	# initialize grid node
	for i in range(n+1):
		row = []
		for j in range(n+1):
			col = []
			for k in range(n+1):
				# initialize grid cell: v_x, v_y, v_z, mass
				col.append(np.array([0.,0.,0.,0.]))
			row.append(col)
		grid.append(row)
	# relate the position of the domain to the grid node
	# add scripted object in the domain
	
def simulate():
	cur_grid = copy.deepcopy(grid)
	# particle to grid
	for p in py_particles:
		# do something
		base_coord = p.x * inv_dx - vec(0.5)
		fx = p.x * inv_dx - base_coord
		w[3] = [
			vec(0.5) * Math.pow(vec(1.5) - fx, 2),
			vec(0.75) - Math.pow(fx - vec(1.0), 2),
			vec(0.5) * Math.pow(fx - vec(0.5), 2)
		]
		e = hardening * (1.0 - p.Jp)
		mu = mu_0 * e
		da = lambda_0 * e
		J = np.linalg.det(p.F)
		r,s = np.polar(p.F)
		Dinv = 4 * inv_dx * inv_dx;
		PF = (2 * mu * (p.F-r) * transposed(p.F) + da * (J-1) * J)
		stress = - (dt * vol) * (Dinv * PF)
		affine = stress + particle_mass * p.C
		# mapping to grid
		for i in range(n+1):
			for j in range(n+1):
				for k in range(n+1):
					# do something
					dpos = (np.array([i,j,k]) - fx) * dx
					mass_x_velocity = np.array([p.v[0]*particle_mass,p.v[1]*particle_mass,p.v[2]*particle_mass,particle_mass])
					# probably need to fix this line as don't know if np.array can be accessed as .x .y
					cur_grid[base_coord[0] + i][base_coord[1] + j][base_coord[2] + k] += w[i][0]*w[j][1]*w[k][2] * (mass_x_velocity + np.array([affine[0]*dpos,affine[1]*dpos,affine[2]*dpos,0]))
	# for all grid node
	for i in range(n+1):
		for j in range(n+1):
			for k in range(n+1):
				# do something
				g = cur_grid[i][j][k]
				if g[3] > 0:
					g /= g[2]
					g += dt * np.array([0,-200,0])
					boundary = 0.05
					# change the coordinate here?
					x = i / n
					y = j / n
					z = k / n
					if x < boundary or x > 1-boundary or y < boundary or y > 1-boundary or z < boundary or z > 1-boundary:
						g = vec(0.)
					
	# grid to particle
	for p in py_particles:
		# do something
		base_coord = p.x * inv_dx - vec(0.5)
		fx = p.x * inv_dx - base_coord
		w[3] = [
			vec(0.5) * Math.pow(vec(1.5) - fx, 2),
			vec(0.75) - Math.pow(fx - vec(1.0), 2),
			vec(0.5) * Math.pow(fx - vec(0.5), 2)
		]
		p.C = mat(0.)
		p.v = vec(0.)
		for i in range(n+1):
			for j in range(n+1):
				for k in range(n+1):
					dpos = np.array([i,j,k]) - fx
					grid_v = cur_grid[base_coord[0] + i][base_coord[1] + j][base_coord[2] + k]
					weight = w[i][0] * w[j][1]
					p.v += weight * grid_v
					p.C += 4 * inv_dx * np.outer(weight * grid_v, dpos)
		p.x += dt * p.v
		F = (mat(1.) + dt * p.C) * p.F;
		svd_u, sig, svd_v = np.linalg.svd(F)
		for i in range(2*plastic):
			sig[i][i][i] = clamp(sig[i][i], 1.0 - 2.5e-2, 1.0f + 7.5e-3)
		oldJ = np.linalg.det(F)
		F = svd_u * sig * np.transpose(svd_v)

		Jp_new = clamp(p.Jp * oldJ / np.linalg.det(F), 0.6f, 20.0f)

		p.Jp = Jp_new
		p.F = F
		

def draw(frame):
	for p in particles:
		p.location = (py_particles.x[0],py_particles.x[1],py_particles.x[2])
	bpy.context.scene.frame_set(frame)
	bpy.ops.anim.keyframe_insert()
	
def vec(val):
	return np.array([val,val,val])
	
def mat(val):
	return np.array([val,val,val],[val,val,val],[val,val,val])
	
def clamp(a, min, max):
	if a < min:
		return min
	if a > max:
		return max
	return a
	
	