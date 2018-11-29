"""
============
3D animation
============

A simple example of an animated plot... In 3D!
"""
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

class TrajVideo:
	"""
	Plot the input trajectory.
	Shape of the trajectory should be (nFrames, nParticles, dim)
	nParticles<Natoms is the number of particles that need to be plotted
	"""

	def __init__(self, traj, L, times=None):
		self.SetTraj(traj, times)
		self.SetL(L)
		return

	def SetTimes(self, times):
		if np.all(times)==None:
			self.times=np.arange(self.nFrames)
		else: 
			self.times=times
		if len(self.times)!=self.nFrames:
			raise ValueError('TrajVideo: len(times)!=len(traj) is not acceptable')
		return

	def SetTraj(self, traj, times):
		self.traj=traj
		self.nFrames, self.nParticles, self.dim = np.shape(traj)
		self.SetTimes(times)
		return

	def SetL(self, L):
		self.L=L
		self.Lm=0.5*L
		return

	def update_traj_lines(self,iframe, mytraj, lines, old_points, time_text):
		for iparticle in range(self.nParticles):
			if old_points: #Whether to plot also the old points, or only the new ones
				lines[iparticle].set_data( mytraj[ :iframe, iparticle, 0:2].transpose() )
				lines[iparticle].set_3d_properties( mytraj[:iframe, iparticle, 2].transpose() )
			else:
				lines[iparticle].set_data( mytraj[ iframe, iparticle, 0:2].transpose() )
				lines[iparticle].set_3d_properties( mytraj[iframe, iparticle, 2].transpose() )

		time_text.set_text('t = %g' % self.times[iframe])
		return lines

	def PlotLines3D(self, marker='o', markersize=10, linestyle='', old_points=True, color=None, savetitle=None, show=True):
		# Attaching 3D axis to the figure
		self.fig = plt.figure()

		ax = p3.Axes3D(self.fig)

		lines = [ax.plot(self.traj[0:1, iparticle,0], self.traj[0:1, iparticle,1], self.traj[0:1, iparticle,2], marker=marker, markersize=markersize, linestyle=linestyle, color=color)[0] for iparticle in range(self.nParticles)]

		# Setting the axes properties
		ax.set_xlim3d([-self.Lm, self.Lm])
		ax.set_xlabel('X')
		ax.set_ylim3d([-self.Lm, self.Lm])
		ax.set_ylabel('Y')
		ax.set_zlim3d([-self.Lm, self.Lm])
		ax.set_zlabel('Z')
		if not savetitle==None: ax.set_title(savetitle)

		time_text = ax.text(x=0.02, y=0.95, z=0.95, s='',fontdict=None, withdash=False, transform=ax.transAxes)

		# Creating the Animation object
		line_ani = animation.FuncAnimation(self.fig, self.update_traj_lines, self.nFrames, fargs=(self.traj, lines, old_points, time_text), interval=500, blit=False, repeat=False)
		
		# Save animation
		if not savetitle==None: 
			# Set up formatting for the movie files
			Writer = animation.writers['ffmpeg']
			writer = Writer(fps=2, metadata=dict(artist='MBJ'), bitrate=1800)
			line_ani.save(savetitle+'.mp4', writer=writer)
		if show: plt.show()



def Gen_RandLine(length, dims=2):
	"""
	Create a line using a random walk algorithm

	length is the number of points for the line.
	dims is the number of dimensions the line has.
	"""
	lineData = np.empty((dims, length))
	lineData[:, 0] = np.random.rand(dims)
	for index in range(1, length):
		# scaling the random numbers by 0.1 so
		# movement is small compared to position.
		# subtraction by 0.5 is to change the range to [-0.5, 0.5]
		# to allow a line to move backwards.
		step = ((np.random.rand(dims) - 0.5) * 0.1)
		lineData[:, index] = lineData[:, index - 1] + step

	return lineData


def update_lines(num, dataLines, lines):
	for line, data in zip(lines, dataLines):
		# NOTE: there is no .set_data() for 3 dim data...
		line.set_data(data[0:2, :num])
		line.set_3d_properties(data[2, :num])
	return lines



#--------------#
#Run an example#
#--------------#
if __name__=="__main__":

	nParticles=10
	dim=3
	nFrames=20

	# Attaching 3D axis to the figure
	fig = plt.figure()
	ax = p3.Axes3D(fig)

	# Fifty lines of random 3-D lines
	data = [Gen_RandLine(nFrames, dim) for index in range(nParticles)]

	# Creating fifty line objects.
	# NOTE: Can't pass empty arrays into 3d version of plot()
	lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

	# Setting the axes properties
	ax.set_xlim3d([0.0, 1.0])
	ax.set_xlabel('X')

	ax.set_ylim3d([0.0, 1.0])
	ax.set_ylabel('Y')

	ax.set_zlim3d([0.0, 1.0])
	ax.set_zlabel('Z')

	ax.set_title('3D Test')

	# Creating the Animation object
	line_ani = animation.FuncAnimation(fig, update_lines, nFrames, fargs=(data, lines),
									   interval=50, blit=False, repeat=False)

	plt.show()
