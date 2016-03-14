from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import seaborn as sns
from scipy import stats
import mpld3
from matplotlib import gridspec
import matplotlib.cm as cm

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

def jointplot(data, save_file):
	values = data[1:,:].T

	kde1 = stats.gaussian_kde(data[:,0])
	kde2 = stats.gaussian_kde(data[:,1])
	kde = stats.gaussian_kde(values, bw_method='silverman')

	# Create a regular 2D grid with 50 points in each dimension
	xmin, ymin = data.min(axis=0)# + np.array([-.5,.5])
	xmax, ymax = data.max(axis=0)# + np.array([-.5,.5])
	xi, yi = np.mgrid[xmin:xmax:500j, ymin:ymax:500j]

	# Create two 1D grid with 50 points in each dimension
	x = np.linspace(xmin,xmax,100)
	y = np.linspace(ymin,ymax,100)
	dx = kde1(x)
	dy = kde2(y)

	# Evaluate the KDE on a regular grid...
	coords = np.vstack([item.ravel() for item in [xi, yi]])
	density = kde(coords).reshape(xi.shape)

	#Set grid
	gs = gridspec.GridSpec(2, 2, width_ratios=[3,1], height_ratios=[1,3])

	#Contour Plot
	fig = plt.figure()
	fig.suptitle('My Contour Plot')
	#ax = plt.subplot()
	ax = plt.subplot(gs[1,0])
	xlim = (np.percentile(data[:,0], 0.5), np.percentile(data[:,0], 99.5))
	ylim = (np.percentile(data[:,1], 0.5), np.percentile(data[:,1], 99.5))
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	cax = ax.contourf(density.T, origin='lower', aspect='auto',cmap=cm.coolwarm, extent=(xmin,xmax,ymin,ymax))
	ax.contour(density.T, origin='lower', aspect='auto',cmap=cm.bone, extent=(xmin,xmax,ymin,ymax))
	#ax.contour(density.T, extent = (xmin,xmax,ymin,ymax))#, origin='lower', aspect='auto', cmap=cm.bone)

	#Marginals
	
	axr = plt.subplot(gs[1,1], sharey=ax,xticks=[],xlim=(0,dy.max()),ylim=(ymin,ymax))
	cbar = plt.colorbar(cax)#, ticks=[0.001, .006,.011])
	cbar.ax.set_yticklabels(['low', 'med', 'high'])
	axt = plt.subplot(gs[0,0], sharex=ax,frameon=False,yticks=[],xlim=(xmin,xmax),ylim=(0,dx.max()))
	axr.plot(dy,y,color='black')
	axt.plot(x,dx,color='black')
	#axr.fill(dy, y, alpha=.75,color='#5673E0')
	#axt.fill(x,dx, alpha=.75,color='#5673E0')

	#Clean Up
	beautify(axr)
	beautify(axt)
	
	#pp = PdfPages("%s/%s_%s_%s.pdf" %(save_dir, main, titles[i], titles[j]))
	pp = PdfPages(save_file)
	#plt.xlabel(titles[i])
	#plt.ylabel(titles[j])
	#plt.title(main)
	pp.savefig()
	pp.close()
	plt.clf()