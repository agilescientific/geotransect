# -*- coding: utf 8 -*-
"""
Functions to read and write an LAS file containing
lithology data as a set of tops, and to process the
data in some simple ways.

For example:

Top, Base, Lithology
240, 260,  Sandstone
260, 290,  Shale
300, 320,  Limestone
320,    ,  Shale

The LAS file will ideally be LAS 2.0 and 3.0
compatible, if possible.
"""

import logging
import re
import StringIO
import csv
import numpy as np
from PIL import Image

import templates
import utils

# Set up logging.
log = logging.getLogger('lithlog')
fstring = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

fh = logging.FileHandler('/tmp/lithlog.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(logging.Formatter(fstring))
log.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
ch.setFormatter(logging.Formatter(fstring))
log.addHandler(ch)


def read_image(filename, hex=True):
    """
    Read an image and return an array of unique floats.

    Args:
       filename (str): The image to read, with path.
       hex (boolean): Do you want hex strings, or RGB triples?

    Returns:
       ndarray. An array of floats representing lithologies.

    """
    im = np.array(Image.open(filename))

    col = im.shape[1]/10.
    rgb = im[:, col, :3]

    if hex:
        return np.array([utils.rgb_to_hex(t) for t in rgb])
    else:
        return rgb


def intervals_from_loglike(loglike, offset=2):
    """
    Take a log-like stream of numbers or strings representing lithologies.

    Args:
       loglike (array-like): The input stream of loglike data.
       offset (int): Offset (down) from top at which to get lithology,
       to be sure of getting 'clean' pixels.

    Returns:
       ndarray. Two ndarrays: tops and lithologies.

    """
    loglike = np.array(loglike)
    all_edges = loglike[1:] == loglike[:-1]
    edges = all_edges[1:] & (all_edges[:-1] == 0)

    tops = np.where(edges)[0]
    tops = np.append(0, tops)

    liths = loglike[tops + offset]

    return tops, liths


def intervals_from_image(filename, start, stop, offset=2):
    """
    Read an image and generate a 'tops' dictionary.

    """
    loglike = read_image(filename)
    pixels, liths = intervals_from_loglike(loglike, offset)

    length = float(loglike.size)

    tops = [start + (p/length) * (stop-start) for p in pixels]

    # This is a horrible hack. Turn it into a CSV so I can use my existing
    # function for reading from an LAS file.
    s = '{0}, , {1}'
    as_list = [s.format(t, liths[i]) for i, t in enumerate(tops)]
    csv_string = ""
    for row in as_list:
        csv_string += row + '\n'

    return intervals_from_string(csv_string)


def find_best_match(colour, legend, tolerance=0.0):
    """
    Given a colour and a list of colours, find the closest colour in the
    legend to the given colour.
    """

    r1, g1, b1 = utils.hex_to_rgb(colour)

    # Start with a best match of black.
    best_match = '000000'
    best_match_dist = np.sqrt(r1**2. + g1**2. + b1**2.)

    # Now compare to each colour in the legend.
    for l in legend:
        r2, g2, b2 = utils.hex_to_rgb(l)

        distance = np.sqrt((r2-r1)**2. + (g2-g1)**2. + (b2-b1)**2.)

        if distance < best_match_dist:
            best_match_dist = distance
            best_match = l

    if best_match_dist <= tolerance:
        return best_match
    else:
        return None


def intervals_from_las3_string(string, dlm=','):
    # Conditioning.
    string = re.sub(r'(\n+|\r\n|\r)', '\n', string.strip())
    string = re.sub(r'None', '', string.strip())

    # Parse... not doing anything special with sections yet.
    as_strings = []
    f = StringIO.StringIO(string)
    reader = csv.reader(f, delimiter=dlm, skipinitialspace=True)
    section = None
    for row in reader:
        if re.search(r'~.+?parameter', row[0].lower()):
            section = "param"
        if re.search(r'~.+?definition', row[0].lower()):
            section = "def"
        if re.search(r'~.+?data ?\|.+?definition', row[0].lower()):
            section = 'data'
        else:
            if section == 'data':
                as_strings.append(row)

    result = {'tops': [], 'bases': [], 'liths': []}

    for i, row in enumerate(as_strings):
        # TOP
        this_top = float(row[0])

        # BASE
        # Base is null: use next top if this isn't the end.
        if not row[1]:
            if i < len(as_strings)-1:
                this_base = float(as_strings[i+1][0])  # Next top.
            else:
                this_base = this_top + 1  # Default to 1 m thick at end.
        else:
            this_base = float(row[1])

        # LITH
        this_lith = {}
        this_lith['lithology'] = row[2].strip()
        this_lith['colour'] = row[3].strip()
        this_lith['grainsize'] = row[4].strip()
        this_lith['map colour'] = row[5].strip()

        # If this top is not the same as the last base, add a layer first.
        if (i > 0) and (this_top != result['bases'][-1]):
            l = "Interpolating gap: {0} to {1}"
            log.info(l.format(result['bases'][-1], this_top))
            result['tops'].append(result['bases'][-1])  # Last base.
            result['bases'].append(this_top)
            result['liths'].append('')

        # ASSIGN
        result['tops'].append(this_top)
        result['bases'].append(this_base)
        result['liths'].append(this_lith)

    return result  # A striplog dict for this well


def intervals_from_string(string, dlm=','):
    """
    Convert an OTHER field containing interval data to an interval dictionary.
    Handles an arbitrary number of fields; it's up to you to know what they
    are.

    Args:
       string (str): The input string, given by ``well.other``.
       dlm (str): The delimiter, given by ``well.dlm``. Default: CSV

    Returns:
       dict. A dictionary with keys 'tops', 'bases', and 'liths'.

    Example:
        # Lithology interval data
        ~OTHER
        # TOP       BOT        LITH
          312.34,   459.61,    225
          459.71,   589.61,    200
          589.71,   827.50,    225
          827.60,   1010.84,   225

    TODO:
       Space delimiter is hard, b/c it's hard to count the null entries.

    """
    string = re.sub(r'(\n+|\r\n|\r)', '\n', string.strip())

    as_strings = []
    f = StringIO.StringIO(string)
    reader = csv.reader(f, delimiter=dlm, skipinitialspace=True)
    for row in reader:
        as_strings.append(row)

    result = {'tops': [], 'bases': [], 'liths': []}

    for i, row in enumerate(as_strings):
        # TOP
        this_top = float(row[0])

        # BASE
        # Base is null: use next top if this isn't the end.
        if not row[1]:
            if i < len(as_strings)-1:
                this_base = float(as_strings[i+1][0])  # Next top.
            else:
                this_base = this_top + 1  # Default to 1 m thick at end.
        else:
            this_base = float(row[1])

        # LITH
        this_lith = row[2]

        # If this top is not the same as the last base, add a layer first.
        if (i > 0) and (this_top != result['bases'][-1]):
            l = "Interpolating gap: {0} to {1}"
            log.info(l.format(result['bases'][-1], this_top))
            result['tops'].append(result['bases'][-1])  # Last base.
            result['bases'].append(this_top)
            result['liths'].append('')

        # ASSIGN
        result['tops'].append(this_top)
        result['bases'].append(this_base)
        result['liths'].append(this_lith)

    return result


def write_other(lith, dlm=",", source=""):
    """
    Accept a dict of intervals and write them out as an OTHER section.

    Args:
       lith (ndarray): The input intervals dictionary.
       dlm (str): The delimiter.

    Returns:
       str. A string forming the OTHER section of an LAS file.

    """
    result = "# Lithology data: {0}\n~OTHER\n#    TOP          BOT    LITH\n"
    result = result.format(source)

    for i, t in enumerate(lith['tops']):
        r = '{0:12.3f}'.format(t)
        r += '{0}{1:12.3f}'.format(dlm, lith['bases'][i])
        r += '{0}   {1:s}'.format(dlm, lith['liths'][i])
        result += "{0}\n".format(r.strip())

    return result


def write_lithology(lith, dlm=",", source=""):
    """
    Accept a dict of intervals and write them out as an LAS 3.0 section.

    Args:
       lith (ndarray): The input intervals dictionary.
       dlm (str): The delimiter.
       source (str): The sourse of the data.

    Returns:
       str. A string forming the OTHER section of an LAS file.

    """
    data = ''

    for i, t in enumerate(lith['tops']):
        data += '{0:12.3f}'.format(t)
        data += '{0}{1:12.3f}'.format(dlm, lith['bases'][i])
        data += '{0}   {1:16s}'.format(dlm, lith['liths'][i]['Lithology']+dlm)
        data += '{0:16s}'.format(lith['liths'][i]['Colour']+dlm)
        data += '{0:6s}'.format(lith['liths'][i]['Grainsize']+dlm)
        data += '{0:8s}'.format(lith['liths'][i]['Map colour'])
        data += '\n'

    result = templates.lithology

    return result.format(source=source, data=data)


def discretize_values(a, no_bins, scale=None):
    """
    Reduce the number of values in an array and optionally map to some new
    (min, max) scale. The new values are in the centres of the bins.

    Args:
       a (ndarray): The input log array.
       no_bins (int): Number of bins.

    Kwargs:
       scale (tuple): The min and max of the blocky log.

    Returns:
       ndarray. The blocky version of the log.

    TODO:
       Allow 'smart' scale that provides some sensible round-number divisions.
       E.g. for a log with min 171 and max 388, I'd like::

       no_bins = 22
       scale = (170, 390)

    """
    if not scale:
        scale = np.nanmin(a), np.nanmax(a)

    # Define the intervals the data will be binned into.
    # Note that this is the not the range of the output,
    # we will have to rescale for that.
    bins = np.linspace(scale[0], scale[1], no_bins)
    bin_half_width = 0.5 * (bins[1] - bins[0])

    blocky = np.digitize(a, bins)                           # Unitless.
    blocky = blocky * (scale[1] - scale[0]) / (no_bins-1)   # Rescale.
    blocky = blocky + scale[0] - bin_half_width             # Shift.

    blocky *= a/a   # Workaround to restore the nans in the original.

    return blocky


def compress_log(a, depth, shift=True):
    """
    Takes a regularly sampled log and returns a set of tops where the log
    changes. Essentially this is a kind of data compression. Intended for
    logs with discrete values, not ordinary logs.

    You can optionally choose not to shift the tops by half the depth interval.

    Args:
       a (ndarray): The input log array.
       depth (ndaray): The depth basis of the log.

    Kwargs:
       shift (bool): Whether or not to shift by 0.5*dz.

    Returns:
       list. A list of (depth, value) pairs ('tops').

    """
    edges = a[1:] - a[:-1]
    values = a[1:][edges != 0]
    values = np.append(a[0], values)

    if shift:
        s = 0.5 * (depth[1] - depth[0])
    else:
        s = 0
    tops = depth[:-1][edges != 0] + s
    tops = np.append([depth[0]], tops)

    bots = np.append(tops[1:], depth[-1])

    return np.array(zip(tops, bots, values))


def undiscretize(a, depth):
    """
    Return a fully sampled log from a set of interval data, given a start
    depth and depth interval.
    """

    result = np.zeros_like(depth)     # Make a container for the result
    init = depth[0]                   # The top depth of the logs, to shift z
    dz = depth[1] - depth[0]          # The increment

    for i in a:
        top_index = np.ceil((i[0]-init)/dz)
        base_index = np.ceil((i[1]-init)/dz)+1
        result[top_index:base_index] = i[2]

    return result

