#!/usr/bin/python
"""
The tsplib.py module imports TSPLIB XML Files and prints the resulting Adjacency Matrix.
(http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/)

Usage: tsp.py burma14.xml

Author: Florin Schimbinschi (florinsch@gmail.com)
"""


def read_tsplib(file_name):
    """
    This function parses a TSPLIB XML file and returns an [[Adjacency Matrix]]
    containing costs from one vertex to another stored as double. M[i][j] = 3.141592

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
    except ET.ParseError, e:
            print 'There was a problem parsing', fileName, ':\n', e
            return
    except IOError, e:
            print 'There was a problem opening the file:\n', e
            return
    
    # get root
    root = tree.getroot()

    # problem properties dictionary
#        graph_props = dict()
#        for child in root:
#                graph_props[child.tag] = child.text
#        del graph_props['graph'] # remove graph property
#        graph_props['type'] = root.tag # add problem type
    
#        print 'Problem properties:'
#        for prop, value in graph_props.items():
#                print prop.capitalize(), '->', value

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
            print 'One of the values of the graph attributes is not valid.'
            print 'Hint:', idx_from, '->', idx_to, '=', edge.get('cost')
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


def _print_matrix(mat):
    """
    Prints an Adjacency Matrix to console.

    Args:
        mat ([[]]): A square matrix.
    """
    print ' ',
    for i in range(len(mat[1])):
          print i,
    print
    for i, element in enumerate(mat):
          print i, ''.join(str(element))

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
            _print_matrix( _symmetric_tril( read_tsplib(sys.argv[1]) ))
            #import PyGMO.problem
            #tsp = PyGMO.problem.tsp( read_tsplib(sys.argv[1]) )
            #print tsp
    else:
            print 'No file names given as argument.'
            print __doc__
