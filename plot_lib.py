import numpy as np

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

from matplotlib import pyplot as plt


def colorline(x, y, z=None, cmap=plt.get_cmap('copper'),
              norm=True, linewidth=5, alpha=1.0):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''
    
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = y
           
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
        
    z = np.asarray(z)
    
    if norm:
        norm = plt.Normalize(np.amin(y), np.amin(y))
    
    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)
    
    #ax = plt.gca()
    #ax.add_collection(lc)
    
    return lc


def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments


def elevation_plot(x, elevs, xlim):

    
    colored_line = colorline(x, elevs)
    
    plt.gca().add_collection(colored_line)

    bbox = {'fc':'w', 'pad':0, 'ec' :'none', 'alpha' : 0.5}  
    props = {'ha':'left', 'va':'center', 'fontsize': 6, 'bbox':bbox}
    
    plt.ylim(0,3*np.amax(elevs))  
    plt.xlim(xlim)  
    plt.gca().set_frame_on(False)
    plt.yticks([])
    plt.xticks([])
    plt.gca().patch.set_alpha(0.1)
    
    for spine in plt.gca().spines.values():
        spine.set_edgecolor('#F3F3F3')

    plt.grid(True)
    plt.text(0.0, 2*(np.ptp(elevs)), "Elevation [m]",
             props, rotation=0)

