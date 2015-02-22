import numpy as np

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

from matplotlib import pyplot as plt

from matplotlib.colors import hsv_to_rgb, hex2color

import matplotlib.transforms as transforms
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle

from lithology.legends import legend

import csv
import StringIO
import re 


def elevation_plot(elevation, bedrock, xlim):


    ele_x = elevation.coords
    ele_y = elevation.data

    bed_x = bedrock.coords

    # Could be done with LineCollection?
    for i, geo in enumerate(bedrock.data[:-1]):

        lim1 = bed_x[i]
        lim2 = bed_x[i+1]

        idx = np.where(np.logical_and(ele_x >= lim1,
                                      ele_x <= lim2))[0] # list 

        if idx[-1] < ele_x.size-1:
            idx = np.append(idx, (idx[-1] + 1))

        hsv = np.array([[geo["AV_HUE"],geo["AV_SAT"],
                         geo["AV_VAL"]]]).reshape(1,1,3)
        
        color = hsv_to_rgb(hsv/255.)
        plt.plot(ele_x[idx], ele_y[idx], linewidth=3,
                 color=color.flatten())
        

    bbox = {'fc':'w', 'pad':0, 'ec' :'none', 'alpha' : 0.5}  
    props = {'ha':'left', 'va':'center', 'fontsize': 6, 'bbox':bbox}
    
    plt.ylim(0,3*np.amax(ele_y))  
    plt.xlim(xlim)  
    plt.yticks([])
    plt.xticks([])
    plt.gca().patch.set_alpha(0.1)


    plt.grid(True)
    plt.text(0.0, 1.2*(np.ptp(ele_y)), "Surface Geology",
             props, rotation=0)
    plt.gca().set_frame_on(False)



def plot_striplog(ax, striplog, width=1,
                  ladder=False, summaries=False, minthick=1,
                  alpha=0.75):


    
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
                         edgecolor='k', 
                         linewidth=1.0, alpha=alpha)
        ax.add_patch(rect)

        if summaries:
            if thick >= minthick:
                ax.text(6, t+thick/2, summarize(l),
                        verticalalignment='center')

        ax.set_ylim([striplog['bases'][-1],0])


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
                params['hexcolor']= (row['hexcolor'])
                params['fill_left_cond']= bool(row['fill_left_cond'])
                params['fill_left']= (row['fill_left'])
                params['fill_right_cond']= bool(row['fill_right_cond'])
                params['fill_right']=(row['fill_right'])
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


def composite_plot(fig, gs,striplog, logs):

    fname = 'templates/Petrophysics_display_template.csv'
    params = get_curve_params('RT_HRLT', fname)


    Z = logs.data['DEPT']
    GR = logs.data['GR']
    DT = logs.data['DT'] 
    DPHISS = logs.data['DPHI_SAN']
    NPHISS = logs.data['NPHI_SAN']
    DTS = logs.data['DTS']
    RT_HRLT = logs.data['RT_HRLT']
    RHOB = logs.data['RHOB']
    DRHO = logs.data['DRHO']

    log_dict = {
            'GR':GR, 
            'DT':DT, 
            'DPHI_SAN':DPHISS,
            'NPHI_SAN':NPHISS,
            'DTS':DTS, 
            'RT_HRLT':RT_HRLT,
            'RHOB':RHOB,
            'DRHO':DRHO
            }
    
    left  = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.1   # the amount of width reserved for blank space between subplots
    hspace = 0.5   # the amount of height reserved for white space between subplots

    window = 51 # window length for smoothing must be an odd integer
    frac = 0.05
    ntracks = 5
    lw = 1.0
    smooth = True
    has_striplog = True
    height = 2.5*ntracks  # in inches
    width = 1.5*ntracks # in inches
    fs = 12  #font size for curve labels


    naxes = 0
    ncurv_per_track = np.zeros(ntracks)

    if has_striplog:
        ncurv_per_track[0]=1

    for curve, values in log_dict.iteritems(): 
        naxes += 1
        params =  get_curve_params(curve, fname)
        ncurv_per_track[params['track']] += 1

    
    #fig, axs = plt.subplots(1, ntracks, figsize=(width, height))
    #fig.subplots_adjust(left, bottom, right, top, wspace, hspace)
    #fig.set_facecolor('w')
    
    axss = plt.subplot(gs[2:12,0])
    axs0 = [axss, axss.twiny()]
    axs1 = [plt.subplot(gs[2:12,1])]
    axs2 = [plt.subplot(gs[2:12,2])]
    axs3 = [plt.subplot(gs[2:12,3])]
    axs4 = [plt.subplot(gs[2:12,4])]

    axs = [axs0, axs1, axs2, axs3, axs4]

    plot_striplog(axs0[0], striplog,width = 5, alpha = 0.75,
                  ladder=True)

    axs0[0].set_ylim([striplog['bases'][-1], 0])

    # Plot each curve with a white fill to fake the curve fill.

    label_shift = np.zeros(len(axs))

    for curve, values in log_dict.iteritems():

        params = get_curve_params(curve, fname)

        i = params['track']

        if i == 0:
            j = 1

        j = 0  # default number of tracks to index into

        ncurves = ncurv_per_track[i]

        label_shift[i] += 1

        units = '$%s$' % params['units']

        #crvtxt = params['acronymn']+ '\n' + '\n' + units

        linOrlog = params['logarithmic']

        sxticks = np.array(params['xticks'])
        xticks = np.array(sxticks, dtype = float)
        whichticks = 'major'

        if linOrlog == 'log':
            midline = np.log(np.mean(xticks))
            xpos = midline 
            whichticks = 'minor'
        else:
            midline = np.mean(xticks)
            xpos = midline 

        if smooth == True:
            values = rolling_window(values,window)

        if params['acronymn']=='GR':
            j = 1  # second axis in first track
            label_shift[i] = 1
        #if params['acronymn']=='RHOB':
        #    j = 2
        #    label_shift[i] = 1

        #fill_left
        if params['fill_left_cond']:
            print params['fill_left_cond']
            if params['acronymn']=='GR':
                # do the fill for the lithology track
                axs[i][j].fill_betweenx(Z, params['xleft'], values,
                                facecolor=params['fill_left'], alpha = 1.0, 
                                zorder = 11)

            if params['acronymn']=='DPHI_SAN':
                # do the fill for the neutron porosity track
                nphi = rolling_window(log_dict['NPHI_SAN'],window)
                axs[i][j].fill_betweenx(Z, 
                                        nphi,
                                        values,
                                        where=nphi>=values,
                                        facecolor=params['fill_left'], 
                                        alpha = 1.0, 
                                        zorder = 11)

                axs[i][j].fill_betweenx(Z, 
                                        nphi,
                                        values,
                                        where=nphi<=values,
                                        facecolor='#8C1717', 
                                        alpha = 0.5, 
                                        zorder = 12)

        if params['acronymn']=='DRHO':
            blk_drho = 3.2
            values += blk_drho   # this is a hack to get DRHO on RHOB scale
            axs[i][j].fill_betweenx(Z, 
                                    blk_drho,
                                    values,
                                    where=nphi<=values,
                                    facecolor='#CCCCCC', 
                                    alpha = 0.5, 
                                    zorder = 12)

        # fill right        
        if params['fill_right_cond']:
            axs[i][j].fill_betweenx(Z, values, params['xright'],
                                facecolor=params['fill_right'], 
                                alpha = 1.0, zorder = 12)
        # plot curve        
        axs[i][j].plot(values, Z, color = params['hexcolor'],
                       lw = lw, zorder = 13)

        # set scale of curve
        axs[i][j].set_xlim([params['xleft'],params['xright']])
        print 'xticks', xticks
        #print 'strings', sxticks
        axs[i][j].xaxis.set_ticks([xticks])
        axs[i][j].set_xticklabels([sxticks])

        for label in axs[i][j].axes.xaxis.get_ticklabels():
            label.set_rotation(90)

        axs[i][j].tick_params(axis='x', direction='out')

        # set curve label and units beneath it
        # if this is the first curve on the axis

        # curve label

        trans = transforms.blended_transform_factory(axs[i][j].transData, axs[i][j].transData)
        axs[i][j].text(xpos, -120 - 40*(label_shift[i] - 1), params['acronymn'],  
                       horizontalalignment='center',                    
                       verticalalignment='bottom',
                       fontsize=12, color=params['hexcolor'],
                       transform=trans)
        # curve units
        if label_shift[i] <= 1:
            axs[i][j].text(xpos, -100, units,
                       horizontalalignment='center',
                       verticalalignment='top',
                       fontsize=12, color='k',
                       transform=trans)

        axs[i][j].set_xscale(linOrlog)
        axs[i][j].set_ylim([striplog['bases'][-1], 0]) 

        axs[i][j].axes.xaxis.set_ticklabels(xticks)

        axs[i][j].xaxis.tick_top()
        axs[i][j].xaxis.set_label_position('top') 

        axs[i][j].xaxis.grid(True, which=whichticks,
                             linewidth=0.5, linestyle='-', 
                             color='0.85', zorder = 100)

        axs[i][j].yaxis.grid(True, which=whichticks,
                             linewidth=0.5, linestyle='-',
                             color='0.85', zorder = 100)

        axs[i][j].yaxis.set_ticks(np.arange(0,max(Z),100))
        if i != 0:
            axs[i][j].set_yticklabels("")

        # add Depth label

    #axs[0][0].set_ylabel('MD \n $m$',color=hgrey, top)
    axs[0][0].text(-0.25, 50, 'MD \n $m$ ', fontsize = '9',
                   horizontalalignment='left',
                   verticalalignment='center',)
    axs[0][0].axes.yaxis.get_ticklabels()
    axs[0][0].axes.xaxis.set_ticklabels('')
    for label in axs[0][0].axes.yaxis.get_ticklabels():
            label.set_rotation(90)
            label.set_fontsize(10)
    
