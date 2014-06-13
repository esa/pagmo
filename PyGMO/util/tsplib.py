#!/usr/bin/python
"""
The tsplib.py module imports TSPLIB XML Files and prints the resulting Adjacency Matrix.
(http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/)

Usage: tsplib.py burma14.xml

Author: Florin Schimbinschi (florinsch@gmail.com)
"""

import xml.etree.ElementTree as ET

def import_XML_TSPLIB(file_name):
        """
import_XML_TSPLIB parses a TSPLIB XML file and returns an [[Adjacency Matrix]]
containing costs from one vertex to another stored as double. M[i][j] = 3.141592

Args:
        file_name (string): The XML file to be opened for parsing.
Returns:
        adj_mat (double): Adjacency Matrix, 0 per diagonal.
Raises:
        IOError: The input file was not found.
        TypeError: At least one of the attributes in an edge 
                of the XML file is missing or of the wrong type.
        """
        try:
                tree = ET.parse(file_name)
        except ET.ParseError, e:
                print "There was a problem parsing", fileName, ":\n", e
                return
        except IOError, e:
                print "There was a problem opening the file:\n", e
                return
        
        # get root
        root = tree.getroot()

        # problem properties dictionary
        graph_props = dict()
        for child in root:
                graph_props[child.tag] = child.text
        del graph_props['graph'] # remove graph property
        graph_props['type'] = root.tag # add problem type

        # graph data (stored in dict with tuple as keys)
        adj_mat = []
        try: # in case vertex.cost attribute is not set or incorrect type
                for idx_from, vertice in enumerate(root.iter('vertex')):
                        tmp = []
                        for idx_to, edge in enumerate(vertice):
                                tmp.append(float(edge.get('cost')))
                                if idx_from == idx_to: # insert diagonal 0's
                                       tmp.insert(idx_to, 0)
                        adj_mat.append(tmp)
                adj_mat[idx_to+1].append(0) # last diag element 0
        except TypeError:
                print 'One of the values in the graph distances is not valid.'
                print 'Hint:', idx_from, '->', idx_to, '=', edge.get('cost')
                return
        
#        print 'Problem properties:'
#        for prop, value in graph_props.items():
#                print prop.capitalize(), '->', value
        
        return adj_mat
        
def print_matrix(mat):
        """
Prints an Adjacency Matrix to console.
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
                print_matrix( import_XML_TSPLIB(sys.argv[1]) )
                import PyGMO.problem
                tsp = PyGMO.problem.tsp( import_XML_TSPLIB(sys.argv[1]) )
                print tsp
        else:
                print 'No file names given as argument.'
                print __doc__
