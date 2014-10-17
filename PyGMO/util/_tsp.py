#!/usr/bin/python
"""
The tsplib.py module contains helper routines instantiate TSP problems.
"""

def read_tsplib(file_name):
    """
    This function parses an XML file defining a TSP (from TSPLIB
    http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/)
    and returns the Adjacency Matrix that cna be then used to construct a PyGMO.problem.tsp

    Args:
            file_name (string): The XML file to be opened for parsing.
    Returns:
            adj_mat (double): Adjacency Matrix, 0 per diagonal.
    Raises:
    IOError: 
            The input file was not found.
    TypeError: 
            At least one of the attributes in an edge 
            of the XML file is missing or of the wrong type.
    xml.etree.ElementTreeParseError:
            There was an error parsing the file.
            See: https://docs.python.org/2.7/library/xml.etree.elementtree.html
            
    """
    import xml.etree.ElementTree as ET
    try:
            tree = ET.parse(file_name)
    except ET.ParseError as e:
            print('There was a problem parsing', fileName, ':\n', e)
            return
    except IOError as e:
            print('There was a problem opening the file:\n', e)
            return
    
    # get root
    root = tree.getroot()

    # graph data (list of list: [[]])
    adj_mat = []
    symmetric = False
    try: # in case vertex.cost attribute is not set or incorrect type
            for idx_from, vertice in enumerate(root.iter('vertex')):
                    tmp = []
                    for idx_to, edge in enumerate(vertice):
                            # symmetric problems don't have values for main diagonal
                            if idx_from == idx_to != int(edge.text): # insert diagonal 0's
                                   tmp.append(0)
                                   symmetric = True
                            tmp.append(float(edge.get('cost')))
                    adj_mat.append(tmp)
            if symmetric:
                    adj_mat[idx_to+1].append(0) # set last diagonal element to 0
    except TypeError:
            print('One of the values of the graph attributes is not valid.')
            print('Hint:', idx_from, '->', idx_to, '=', edge.get('cost'))
            return

    return adj_mat


def _symmetric_tril(mat):
    """
    If a matrix is symmetric, returns a copy with elements above the main diagonal zeroed.

    Args:
            mat ([[]]]): A square matrix.
    Returns:
    tril ([[]]): The matrix with only lower diagonal items.
    """
    import numpy
    mat = numpy.array(mat)
    
    if (mat.transpose() == mat).all():
            mat = numpy.tril(mat)
    
    return mat.tolist()


def _print_matrix(mat, show_all = False):
    import numpy
    numpy.set_printoptions(linewidth=100)
    numpy.set_printoptions(precision=3)
    # this forces to print all elements on a long row, on the next line
    # otherwise, center elements are snipped '...,' to fit line of 100
    if show_all: 
            numpy.set_printoptions(threshold='nan')

    print(numpy.array(mat))

from PyGMO import *



def _three_impulse(pl1,pl2):
    from math import cos,sqrt,sin
    MU=pl1.mu_central_body

    a1 = pl1.orbital_elements[0]
    i1 = pl1.orbital_elements[2]
    W1 = pl1.orbital_elements[3] 
    e1 = pl1.orbital_elements[1]

    a2 = pl2.orbital_elements[0]
    i2 = pl2.orbital_elements[2]
    W2 = pl2.orbital_elements[3] 
    e2 = pl2.orbital_elements[1]
        
    ra1 = a1 * (1 + e1); # radius of apocenter starting orbit (km)
    ra2 = a2 * (1 + e2); # radius of apocenter target orbit(km)
    cosiREL = cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(W1)*cos(W2)+ sin(i1)*sin(i2)*sin(W1)*sin(W2);
    if cosiREL > 1 or cosiREL < -1:
        cosiREL = 1
    rp2 = a2 * (1 - e2); # radius of apocenter target orbit(km)
    rp1 = a1 * (1 - e1);

    if ra1>ra2: #Strategy is Apocenter-Pericenter
        Vi = sqrt(MU*(2/ra1-1/a1));
        Vf = sqrt(MU*(2/ra1-2/(rp2+ra1)));
        DV1 = sqrt(Vi**2+Vf**2-2*Vi*Vf*cosiREL);   #Change Inclination + pericenter change
        DV2 = sqrt(MU) * abs(sqrt(2/rp2-2/(rp2+ra1))-sqrt(2/rp2-2/(rp2+ra2))); #Apocenter Change
    else:  #(ra1<ra2) Strategy is Pericenter-Apocenter   
        DV1 = sqrt(MU) * abs(sqrt(2/rp1-2/(rp1+ra1))-sqrt(2/rp1-2/(rp1+ra2))); #Apocenter Raise
        Vi = sqrt(MU*(2/ra2-2/(rp1+ra2)));
        Vf = sqrt(MU*(2/ra2-1/a2));
        DV2 = sqrt(abs(Vi*Vi+Vf*Vf-2*Vi*Vf*cosiREL));   #Change Inclination + apocenter change
    return DV1+DV2

def tle2tsp(tlefilename, verbose=False):
    from PyKEP import planet_tle
    planet_list = []

    with open(tlefilename,'r') as f:
        for line in f:
            line1 = f.next()[:69]
            line2 = f.next()[:69]
            planet_list.append( planet_tle(line1,line2) )

    weights = []
    for source_idx,source in enumerate(planet_list):
        row = []
        if verbose:
            print("\rTrying source debris N.: {0}".format(source_idx))
        count = 0
        for target_idx,target in enumerate(planet_list):
            if source == target:
                row.append(0)
            else:
                dist = _three_impulse(source,target)
                row.append(dist)
        weights.append(row)
    return weights