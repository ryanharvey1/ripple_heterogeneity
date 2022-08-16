import numpy as np
from skimage import measure
import scipy.ndimage as ndimage

# detect peaks and countour coordinates in 2D histogram from scratch
def detect_peaks(image):
    """
    Takes an image and detect the peaks using the local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an 8-connected neighborhood
    neighborhood = generate_binary_structure(2,2)
    # apply the local maximum filter; all pixel of maximal value
    # in their neighborhood are set to 1
    local_max = maximum_filter(image, footprint=neighborhood)==image
    # (note: the mask is in fact a numpy array)
    # find the pixels that are the maximum in their neighborhood
    # but are not connected to a high pixel
    # (i.e. the peaks)
    return local_max

def generate_binary_structure(n, m):
    """
    Generates a binary structure.
    """
    return np.array([[0]*m]*n)

def maximum_filter(image, footprint):
    """
    Returns the maximum of the image within the footprint.
    """
    # compute the maximum of the image within the footprint
    return ndimage.maximum_filter(image, footprint)

def detect_contours(image, level):
    """
    Takes an image and a level and returns the contours at that level.
    """
    # compute the contour at the desired level
    contours = find_contours(image, level)
    # simplify the contours
    contours = [simplify_contour(c) for c in contours]
    # return the contours
    return contours

def find_contours(image, level):
    """
    Finds the contours at the given level.
    """
    # find the contours at the given level
    contours = measure.find_contours(image, level)
    # return the contours
    return contours

def simplify_contour(contour):
    """
    Takes a contour and returns a simplified version.
    """
    # compute the Douglas-Peucker approximation of the contour
    contour = douglas_peucker(contour, epsilon=1)
    # return the simplified contour
    return contour

def douglas_peucker(contour, epsilon):
    """
    Returns the simplified version of the contour using the Douglas-Peucker algorithm.
    """
    # compute the perpendicular distance to the line between every two points
    # and select the point that is farthest from the line
    d = perpendicular_distance(contour)
    index = np.argmax(d)
    # if the farthest point is too far, return the contour unchanged
    if d[index] > epsilon:
        return contour
    # otherwise, recursively simplify the contour
    contour1 = douglas_peucker(contour[:index+1], epsilon)
    contour2 = douglas_peucker(contour[index:], epsilon)
    return contour1 + contour2

def perpendicular_distance(contour):
    """
    Returns the perpendicular distance to the line between every two points.
    """
    # compute the perpendicular distance to the line between every two points
    # and select the point that is farthest from the line
    d = np.empty(len(contour))
    for i in range(len(contour)):
        if i == 0:
            d[i] = np.linalg.norm(contour[i] - contour[-1])
        else:
            d[i] = np.linalg.norm(contour[i] - contour[i-1])
    return d

