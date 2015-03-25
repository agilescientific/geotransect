#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defines a hillshade function that's a bit nicer than the one in matplotlib.

Not currently being used in geotransect, but it will be.

http://rnovitsky.blogspot.ca/2010/04/using-hillshade-image-as-intensity.html

:copyright: 2015 Ran Novitsky Nof
:license: Unclear, but expressed as open.
"""
from pylab import gradient, pi, sin, hypot, arctan, arctan2, cos, cm


def set_shade(a,
              intensity=None,
              cmap=cm.jet,
              scale=10.0,
              azdeg=165.0,
              altdeg=45.0):
    '''
    Sets shading for data array based on intensity layer
    or the data's value itself.

    Args:
    a - a 2-d array or masked array
    intensity - a 2-d array of same size as a (no chack on that)
                  representing the intensity layer. if none is given
                  the data itself is used after getting the hillshade values
                  see hillshade for more details.
    cmap - a colormap (e.g matplotlib.colors.LinearSegmentedColormap
            instance)
    scale,azdeg,altdeg - parameters for hilshade function see there for
            more details
    output:
    rgb - an rgb set of the Pegtop soft light composition of the data and
         intensity can be used as input for imshow()
    based on ImageMagick's Pegtop_light:
    http://www.imagemagick.org/Usage/compose/#pegtoplight
    '''

    if intensity is None:
        # hilshading the data
        intensity = hillshade(a, scale=10.0, azdeg=165.0, altdeg=45.0)
    else:
        # or normalize the intensity
        mi, ma = intensity.min(), intensity.max()
        intensity = (intensity - mi)/(ma - mi)
    # get rgb of normalized data based on cmap
    rgb = cmap((a-a.min())/float(a.max()-a.min()))[:, :, :3]
    # form an rgb eqvivalent of intensity
    d = intensity.repeat(3).reshape(rgb.shape)
    # simulate illumination based on pegtop algorithm.
    rgb = 2*d*rgb+(rgb**2)*(1-2*d)
    return rgb


def hillshade(data, scale=10.0, azdeg=165.0, altdeg=45.0):
    '''
    Convert data to hillshade based on matplotlib.colors.LightSource class.

    Args:
        data - a 2-d array of data
        scale - scaling value of the data. higher number = lower gradient
        azdeg - where the light comes from: 0 south ; 90 east ; 180 north ;
                        270 west
        altdeg - where the light comes from: 0 horison ; 90 zenith

    Returns:
        a 2-d array of normalized hilshade
    '''
    # convert alt, az to radians
    az = azdeg*pi/180.0
    alt = altdeg*pi/180.0
    # gradient in x and y directions
    dx, dy = gradient(data/float(scale))
    slope = 0.5*pi - arctan(hypot(dx, dy))
    aspect = arctan2(dx, dy)
    az = -az - aspect - 0.5*pi
    intensity = sin(alt)*sin(slope) + cos(alt)*cos(slope)*cos(az)
    mi, ma = intensity.min(), intensity.max()
    intensity = (intensity - mi)/(ma - mi)
    return intensity
