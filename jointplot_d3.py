from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import seaborn as sns
from scipy import stats
import mpld3
from matplotlib import gridspec
import matplotlib.cm as cm
from matplotlib import animation
from functools import partial
import multiprocessing as mp

def beautify(ax):
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(0)
        line.set_color("gray")
        line.set_markeredgewidth(1.4)

    for line in ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True):
        line.set_markersize(0)
    for line in ax.xaxis.get_ticklines(minor=False) + ax.yaxis.get_ticklines(minor=False):
        line.set_markersize(0)

def jointplots(data, save_dir, titles = None, main = "", refcoords_file = None, axes=None, reshape=True, data_j=None, titles_j=None):

	if(refcoords_file is not None):
		refcoords = load_file(refcoords_file)
		if reshape:
			refcoords = np.transpose(np.vstack(refcoords))
	else:
		refcoords = None
	print((np.shape(refcoords)))
	print(refcoords)

	num_columns = np.shape(data)[1]
	jointplot_partial = partial(jointplot, main = main, save_dir = save_dir, titles = titles, all_data = data, axes=axes, data_j=data_j, titles_j=titles_j)
	pool = mp.Pool(mp.cpu_count())
	pool.map(jointplot_partial, range(0,num_columns))
	pool.terminate()
	#for i in range(0,num_columns):
	#	jointplot_partial(i)

	print("Done plotting columns")
	return

def jointplot(i, all_data, save_dir, make_animation=False, trajectory=None, video_file=None, titles=None,
			  main=None, include_1d_kde=False, custom_xlim=None, custom_ylim=None, axes=None, data_j=None, titles_j=None):
	if data_j is None:
		j_range = range(i+1, all_data.shape[1])
	else:
		j_range = range(0, data_j.shape[1])

	for j in j_range:
		if data_j is not None: 
			data = np.column_stack([all_data[:,i], data_j[:,j]])
			partial_titles = [titles[i], titles_j[j]]
		else:
			data = all_data[:,[i,j]]
			partial_titles = [titles[i], titles[j]]

		values = data[:,].T
		print(np.shape(values))

		kde1 = stats.gaussian_kde(data[:,0])
		kde2 = stats.gaussian_kde(data[:,1])
		kde = stats.gaussian_kde(values, bw_method='silverman')

		print("Computed 2D KDE")

		# Create a regular 2D grid with 50 points in each dimension
		xmin, ymin = data.min(axis=0)# + np.array([-.5,.5])
		xmax, ymax = data.max(axis=0)# + np.array([-.5,.5])

		if custom_xlim is not None:
			xmin, xmax = custom_xlim
			ymin, ymax = custom_ylim

		xi, yi = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]

		# Create two 1D grid with 50 points in each dimension
		x = np.linspace(0.8*xmin, 1.2*xmax,100)		
		y = np.linspace(0.8*ymin,1.2*ymax,100)
		dx = kde1(x)
		dy = kde2(y)

		# Evaluate the KDE on a regular grid...
		coords = np.vstack([item.ravel() for item in [xi, yi]])
		density = kde(coords).reshape(xi.shape)
		density = -0.6 * np.log(density)
		density -= density.min()
		print("Computed Density. Now plotting.")

		#Set grid
		gs = gridspec.GridSpec(2, 2, width_ratios=[3,1], height_ratios=[1,3])

		#Contour Plot
		fig = plt.figure()
		
		#fig.suptitle('My Contour Plot')
		#ax = plt.subplot()
		#ax = plt.axes(xlim=(0,15), ylim=(0,3))
		if include_1d_kde:
			ax = plt.subplot(gs[1,0])
		else:
			ax = fig.add_subplot(111)
		
		"""
		if custom_xlim is None:
			print(density)
			print(density.shape)
			xi[np.argmin(np.abs(2.4-density))], 
			xlim = [0.75*xmin, 1.25*xmax]
			ylim = [0.75*ymin, 1.25*ymax]
		else:
			xlim = custom_xlim
			ylim = custom_ylim
		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
		"""
		if titles is not None:
			ax.set_xlabel(partial_titles[0])
			ax.set_ylabel(partial_titles[1])
		if main is not None:
			ax.set_title(main)
		
		cax = ax.contourf(density.T, origin='lower', aspect='auto',cmap=cm.coolwarm, extent=(xmin,xmax,ymin,ymax), levels=np.linspace(0.0, 2.5,10), vmax=2.5)
		ax.contour(density.T, origin='lower', aspect='auto',cmap=cm.bone, extent=(xmin,xmax,ymin,ymax), levels=np.linspace(0.0, 2.5,10), vmax=2.5)
		#ax.contour(density.T, extent = (xmin,xmax,ymin,ymax))#, origin='lower', aspect='auto', cmap=cm.bone)	
		
		#Marginals
		
		if include_1d_kde:
			axr = plt.subplot(gs[1,1], sharey=ax,xticks=[],xlim=(0,dy.max()),ylim=(ymin,ymax))
			cbar = plt.colorbar(cax)#, ticks=[0.001, .006,.011])
			#cbar.ax.set_yticklabels(['low', 'med', 'high'])
			axt = plt.subplot(gs[0,0], sharex=ax,frameon=False,yticks=[],xlim=(xmin,xmax),ylim=(0,dx.max()))
			axr.plot(dy,y,color='black')
			axt.plot(x,dx,color='black')
			#axr.fill(dy, y, alpha=.75,color='#5673E0')
			#axt.fill(x,dx, alpha=.75,color='#5673E0')

			#Clean Up
			beautify(axr)
			beautify(axt)

		#if make_animation:
		#	ax.scatter([data[0,0]], [data[0,1]], marker = 's', c='g',s=15)

		fig.savefig("%s/%s_%s_%s.svg" %(save_dir, main, partial_titles[0], partial_titles[1]), format='svg', dpi=1200)
		#pp = PdfPages("%s/%s_%s_%s.pdf" %(save_dir, main, partial_titles[0], partial_titles[1]))
		#pp = PdfPages(save_file)
		#plt.xlabel(titles[i])
		#plt.ylabel(titles[j])
		#plt.title(main)
		#pp.savefig()
		#pp.close()
		
		if make_animation:
			line, = ax.plot([], [], lw=0.6)
			scatter = ax.scatter([], [], marker = 's', c='g',s=15)
			def init():
				line.set_data([], [])
				scatter.set_offsets([])
				return line, scatter,

			def animate(i):
	 			#line = ax.scatter([trajectory[1,0]], [trajectory[1,1]], marker = 's', c='g',s=15)
	 			line.set_data(trajectory[:i,0], trajectory[:i,1])
	 			scatter.set_offsets(trajectory[i,:])
	 			return line, scatter,

	 		anim = animation.FuncAnimation(fig, animate, 
	 									   init_func=init, 
			                               frames=trajectory.shape[0]-1, 
			                               interval=10,
			                               blit=False)
	 		anim.save(video_file, fps=30, 
	          extra_args=['-vcodec', 'h264', 
	                      '-pix_fmt', 'yuv420p'])
		plt.clf()



def test(filename):
	fig = plt.figure()
	ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
	line, = ax.plot([], [], lw=2)

	# initialization function: plot the background of each frame
	def init():
	    line.set_data([], [])
	    return line,

	# animation function.  This is called sequentially
	def animate(i):
	    x = np.linspace(0, 2, 1000)
	    y = np.sin(2 * np.pi * (x - 0.01 * i))
	    line.set_data(x, y)
	    return line,

	# call the animator.  blit=True means only re-draw the parts that have changed.
	anim = animation.FuncAnimation(fig, animate, init_func=init,
	                               frames=200, interval=20, blit=True)

	# save the animation as an mp4.  This requires ffmpeg or mencoder to be
	# installed.  The extra_args ensure that the x264 codec is used, so that
	# the video can be embedded in html5.  You may need to adjust this for
	# your system: for more information, see
	# http://matplotlib.sourceforge.net/api/animation_api.html
	anim.save(filename)#, fps=30, extra_args=['-vcodec', 'libx264'])


