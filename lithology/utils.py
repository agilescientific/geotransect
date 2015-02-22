# -*- coding: utf 8 -*-
"""
Helper functions for the library.

"""


def rgb_to_hex(rgb):
    """
    Utility function to convert (r,g,b) triples to hex.
    http://stackoverflow.com/questions/214359/converting-hex-color-to-rgb-and-vice-versa

    """
    h = '%02x%02x%02x' % tuple(rgb)
    return h.upper()


def hex_to_rgb(hex):
    """
    Utility function to convert hex to (r,g,b) triples.
    http://stackoverflow.com/questions/214359/converting-hex-color-to-rgb-and-vice-versa

    """
    h = hex.strip('#')
    l = len(h)

    return tuple(int(h[i:i+l//3], 16) for i in range(0, l, l//3))


def summarize_rock(rock, fmt='%C %g %l'):
    """
    Given a rock dict and a format string,
    return a summary description of a rock.

    Args:
    rock (dict): A rock dictionary.

    Kwargs:
    fmt (str): Describes the format with
        c = colour
        g = grainsize
        l = lithology
        m = map colour

    Returns:
    str. A summary string.

    Example:
    r = {'code': 4,
         'colour': 'Red',
         'grainsize': 'VF-F',
         'lithology': 'Sandstone',
         'map colour': 'F2FF42',
         'summary': 'Sandstone, F-M, Grey',
         'width': 4}

    summarize(r)  -->  'Red vf-f sandstone'
    """

    fields = {'c': 'colour',
              'g': 'grainsize',
              'l': 'lithology',
              'm': 'map colour'
              }

    flist = fmt.split('%')

    fmt_items = []
    fmt_string = flist.pop(0)

    skip = 0

    for i, item in enumerate(flist):
        if rock[fields[item[0].lower()]]:
            this_item = rock[fields[item[0].lower()]].lower()
            if item.isupper():
                this_item = this_item.capitalize()
            fmt_items.append(this_item)
            fmt_string += '{' + str(i-skip) + '}' + item[1:]
        else:
            skip += 1

    # The star trick lets us unpack from a list.
    summary = fmt_string.format(*fmt_items)

    return summary
