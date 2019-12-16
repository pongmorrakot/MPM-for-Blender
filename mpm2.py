import bpy
import numpy as np
import graphyt
import pyck  # if using pyck for geometry
import VTKWriter # for writing output files

def parse(b_objects, b_domain, b_force, b_gridRes, b_young, b_compression, b_stretch, b_hardening, b_time, b_step):
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
	gravity = np.array(b_force)
	domain = b_domain
	u,v,w = b_domain.dimensions
	x,y,z = b_domain.lication
	gridRes = b_gridRes
	young = b_young
	compression = b_compression
	stretch = b_stretch
	hardening = b_hardening
	time = b_time
	step = b_step