import numpy as N
import scipy as S
import scipy.ndimage

temperature = <whatever> 
# This is the data in your polar grid.
# The 0th and 1st axes correspond to r and ., respectively.
# For the sake of simplicity, . goes from 0 to 2., 
# and r's units are just its indices.

def polar2cartesian(outcoords, inputshape, origin):
    """Coordinate transform for converting a polar array to Cartesian coordinates. 
    inputshape is a tuple containing the shape of the polar array. origin is a
    tuple containing the x and y indices of where the origin should be in the
    output array."""

    xindex, yindex = outcoords
    x0, y0 = origin
    x = xindex - x0
    y = yindex - y0

    r = N.sqrt(x**2 + y**2)
    theta = N.arctan2(y, x)
    theta_index = N.round((theta + N.pi) * inputshape[1] / (2 * N.pi))

    return (r,theta_index)

temperature_cartesian = S.ndimage.geometric_transform(temperature, polar2cartesian, 
    order=0,
    output_shape = (temperature.shape[0] * 2, temperature.shape[0] * 2),
    extra_keywords = {'inputshape':temperature.shape,
        'center':(temperature.shape[0], temperature.shape[0])})
