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
    # to check for numerical input -- this is a hack
    if not hasattr(z, "__iter__"):  
        z = np.array([z])
        
    z = np.asarray(z)
    
    if norm:
        norm = plt.Normalize(np.amin(y), np.amin(y))
    
    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm,
                        linewidth=linewidth, alpha=alpha)
    
    return lc


def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates,
    in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y)
    array
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


def plot_striplog(ax, striplog, legend, width=1, ladder=False,
                  summaries=False, minthick=1, alpha=0.75):
    
    for t, b, l in zip(striplog['tops'],striplog['bases'],
                       striplog['liths']):
        origin = (0, t)
        colour = '#' + l['map colour'].strip()
        thick = b - t
        
        if ladder:
            w = legend[colour[1:]]['width']
        else:
            w = width
            
        rect = Rectangle(origin, w, thick, color=colour,
                         edgecolor='k', linewidth=1.0, alpha=alpha)
        ax.add_patch(rect)

        if summaries:
            if thick >= minthick:
                ax.text(6, t+thick/2, summarize(l),
                        verticalalignment='center')

    return ax


def get_curve_params(abbrev, fname):
    """
    returns a dictionary of petrophysical parameters for 
    plotting purposes
    """
    params = {}
    with open(fname, 'rU') as csvfile:
        reader = csv.DictReader(csvfile) 
        for row in reader:
            if row['acronymn'] == abbrev:
                params['acronymn'] = row['acronymn']
                params['track'] = int(row['track'])
                params['units'] = row['units']
                params['xleft'] = float(row['xleft'])
                params['xright'] = float(row['xright'])
                params['logarithmic'] = row['logarithmic']
                params['hexcolor']= row['hexcolor']
                params['fill_left_cond']= bool(row['fill_left_cond'])
                params['fill_left']= row['fill_left']
                params['fill_right_cond']= bool(row['fill_right_cond'])
                params['fill_right']= row['fill_right']
                params['xticks'] = row['xticks'].split(',')
    return params

def rolling_window(a, window):
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        rolled = np.lib.stride_tricks.as_strided(a, 
                                                 shape=shape, 
                                                 strides=strides)
        rolled = np.median(rolled, -1)
        rolled = np.pad(rolled, window/2, mode='edge')
        return rolled
    
def despike(curve, curve_sm, max_clip): 
    spikes = np.where(curve - curve_sm > max_clip)[0]
    spukes = np.where(curve_sm - curve > max_clip)[0]
    out = np.copy(curve)
    out[spikes] = curve_sm[spikes] + max_clip  # Clip at the max allowed diff
    out[spukes] = curve_sm[spukes] - max_clip  # Clip at the min allowed diff
    return out

