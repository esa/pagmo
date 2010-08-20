#!/usr/bin/python


#usage
# double left click --> add node
# drag nodes --> move nodes
# right click on 2 nodes --> add edge
# left click on an edge --> edit edge's weight

import sys
import math
from PyQt4 import Qt,QtGui,QtCore


class Node(QtGui.QGraphicsItem):
	def __init__(self,parent = None):
		QtGui.QGraphicsItem.__init__(self,parent)
		self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
		#self.setFlag(QtGui.QGraphicsItem.ItemSendsGeometryChanges)
    		self.setCacheMode(QtGui.QGraphicsItem.DeviceCoordinateCache)
     		self.setZValue(-1)
		self.edgeList = []
		global counter
		self.id = counter
		counter = counter + 1
		global nodeList
		nodeList.append(self)
		global scene
		scene.addItem(self)
	def boundingRect(self):
		adjust = 0
		return QtCore.QRectF(-7 - adjust, -7 - adjust, 20 + adjust, 20 + adjust)
	def paint(self, painter, option, widget):
		painter.drawEllipse(-7,-7,20,20)
		painter.drawText(self.boundingRect(), " " + str(self.id))
	def addEdge(self,edge):
		self.edgeList.append(edge)
		edge.adjust()
	def itemChange(self, change, value):
		if(change ==  QtGui.QGraphicsItem.ItemPositionHasChanged ) :
			for edge in self.edgeList :
            			edge.adjust()
		return QtGui.QGraphicsItem.itemChange(self,change,value)
	def mousePressEvent(self,event):
		if event.button() == QtCore.Qt.RightButton: #right click on 2 nodes --> create an edge between them
			global rightPressed
			global firstNode
			if rightPressed:
				if not(firstNode.id == self.id): #avoid to create edge from a node to itself
					global scene
					secondNode = self
					newEdge = Edge(firstNode,secondNode)
					self.addEdge(newEdge)
					rightPressed = 0
			else:
				firstNode = self
				rightPressed = 1

class Edge(QtGui.QGraphicsItem):
	def __init__(self, sourceNode, destNode):
		QtGui.QGraphicsItem.__init__(self)
		self.source = sourceNode
		self.dest = destNode
		sourceNode.addEdge(self)
		destNode.addEdge(self)
		#self.adjust()
		self.arrowsize = 10
		self.weight = 0
		self.isThick = 0
		global scene
		scene.addItem(self)
	def adjust(self):
		self.line = QtCore.QLineF(self.mapFromItem(self.source, 0, 0), self.mapFromItem(self.dest, 0, 0));
		self.prepareGeometryChange()

		if (self.line.length() > 20):
			edgeOffset = QtCore.QPointF((self.line.dx() * 10) / self.line.length(), (self.line.dy() * 10) / self.line.length())
			self.sourcePoint = self.line.p1() + edgeOffset
         		self.destPoint = self.line.p2() - edgeOffset
		else:
         		self.sourcePoint = self.line.p1()
			self.destPoint = self.line.p2()
	def paint(self, painter, option, widget):
		if self.isThick:
			pen = QtGui.QPen()
			pen.setWidth(2)
			painter.setPen(pen)
		self.line = QtCore.QLineF(self.sourcePoint, self.destPoint)
		angle = math.acos(self.line.dx() / self.line.length())
		if self.line.dy() >= 0:
         		angle = 2*math.pi - angle
		destArrowP1 = self.destPoint + QtCore.QPointF(math.sin(angle - math.pi / 3) * self.arrowsize, math.cos(angle - math.pi / 3) * self.arrowsize)
     		destArrowP2 = self.destPoint + QtCore.QPointF(math.sin(angle - math.pi + math.pi / 3) * self.arrowsize, math.cos(angle - math.pi + math.pi / 3) * self.arrowsize)
		painter.drawLine(self.line)
		painter.drawLine(self.line.p2(), destArrowP1)
		painter.drawLine(self.line.p2(), destArrowP2)
		painter.drawText((self.dest.pos()+self.source.pos())/2, str(self.weight))
		self.adjust()
	def shape(self): #create an invisible rectangular shape around the arrow to catch clicks on the edge
		extra = 5
		if (self.dest.x() - self.source.x() == 0):
			angle = math.pi / 2
		else:
			angle = math.atan((self.source.y() - self.dest.y())/ (self.source.x() - self.dest.x()))
		points_list = []
		points_list.append(QtCore.QPointF(self.sourcePoint.x() - extra*math.sin(angle), self.sourcePoint.y() + extra * math.cos(angle)))
		points_list.append(QtCore.QPointF(self.sourcePoint.x() + extra*math.sin(angle), self.sourcePoint.y() - extra * math.cos(angle)))
		points_list.append(QtCore.QPointF(self.destPoint.x()   + extra*math.sin(angle), self.destPoint.y()   - extra * math.cos(angle)))
		points_list.append(QtCore.QPointF(self.destPoint.x()   - extra*math.sin(angle), self.destPoint.y()   + extra * math.cos(angle)))
		points_list.append(QtCore.QPointF(self.sourcePoint.x() - extra*math.sin(angle), self.sourcePoint.y() + extra * math.cos(angle)))
		polygon = QtGui.QPolygonF(points_list)
		path = QtGui.QPainterPath()
		path.addPolygon(polygon)
		return path

	def boundingRect(self):
		penwidth = 1
		extra = (self.arrowsize + penwidth)/2.0
		return QtCore.QRectF(self.sourcePoint, QtCore.QSizeF(self.dest.x() - self.source.x(), self.dest.y() - self.source.y())).normalized().adjusted(-extra, -extra, extra, extra)
	def mousePressEvent(self,event): #right click to an edge -> edit edge's weight
		widget = QtGui.QWidget()
		self.weight = float(QtGui.QInputDialog.getText(widget,"Weight", "Edge's weight")[0])


class GraphScene(QtGui.QGraphicsScene):
	def __init__(self):
		QtGui.QGraphicsScene.__init__(self)
	def mouseDoubleClickEvent(self,mouseEvent):
		if mouseEvent.button() == QtCore.Qt.LeftButton:
			global nodeList
			newNode = Node()
			newNode.setPos(mouseEvent.scenePos())

class MainWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)

        self.resize(640, 480)
        self.setWindowTitle('Traveling Salesman Problem')

        self.save = QtGui.QAction('Save graph', self)
        self.load = QtGui.QAction('Load graph', self)
        self.solve = QtGui.QAction('Solve problem', self)
        self.connect(self.save, QtCore.SIGNAL('triggered()'), save)
        self.connect(self.load, QtCore.SIGNAL('triggered()'), load)
        self.connect(self.solve, QtCore.SIGNAL('triggered()'), solve)

        self.saveButton = self.addToolBar('Save graph')
        self.loadButton = self.addToolBar('Load graph')
        self.solveButton = self.addToolBar('Solve problem')

        self.saveButton.addAction(self.save)
        self.loadButton.addAction(self.load)
        self.solveButton.addAction(self.solve)

def save():
	fileName = QtGui.QFileDialog.getSaveFileName()
	ofile = open(fileName,'w')
	writeMatrixToFile(ofile, generateAdjMatrix(nodeList))

def load():
	fileName = QtGui.QFileDialog.getOpenFileName()
	ifile = open(fileName,'r')
	fromFileToMatrix(ifile)
	
#This could be use to delete a node. Actually the interface still doesn't allow to do that. 
#It is also necessary to edit the generateAdjMatrix function to manage the case in which some node has been deleted
def deleteNode(node):
	global nodeList
	for edge in node.edgeList:
		list = []
		if (edge.source.id == node.id):
			list = edge.dest.edgeList
		else:
			list = edge.source.edgeList
		for e in list:
			if ((e.source.id == node.id) or (e.dest.id == node.id)):
					edge.dest.edgeList.remove(e)
	nodeList.remove(node)

def solve():
	solve_problem(generateAdjMatrix(nodeList))

def generateAdjMatrix(nodeList):
	length = len(nodeList)
	matrix = [[noEdge for col in range(length)] for row in range(length)]
	for node in nodeList:
		matrix[node.id][node.id] = 0
		for edge in node.edgeList:
			matrix[edge.source.id][edge.dest.id] = edge.weight
			matrix[edge.dest.id][edge.source.id] = edge.weight
	return matrix

def writeMatrixToFile(ofile,matrix):
	n = len(matrix[0])
	ofile.write(str(n) + '\n') #in the first line write the number of nodes
	for i in range(n):
		for j in range(n):
			ofile.write(str(matrix[i][j]) + '\n')

def fromFileToMatrix(ifile):
	global nodeList
	global edgeList
	n = int(ifile.readline())
	for i in range(n):
		node = Node()
		node.setPos(100*math.cos((i*2*math.pi)/n), 100*math.sin((i*2*math.pi)/n)) #put all the nodes in a circle
	for i in range(n):
		for j in range(n):
			weight = float(ifile.readline())
			if not(weight == noEdge):
				edge = Edge(nodeList[i],nodeList[j])
				edge.weight = weight

def solve_problem(matrix):
	from PyGMO import problem, algorithm, island
	prob = problem.tsp(matrix)
	algo = algorithm.aco(20)
	isl = island(prob,algo,40)
	isl.evolve(1)
	isl.join()
	global nodeList
	path = isl.population.champion.x
	print "Optimal path: " + str(path)
	print "Path weight: "  + str(isl.population.champion.f[0])
	for i in range(len(nodeList)):
		for edge in nodeList[int(path[i])].edgeList:
			if edge.source.id == nodeList[int(path[i])].id and edge.dest.id == nodeList[int(path[ (i+1)%len(path) ])].id:
				edge.isThick = 1
				break
			if edge.dest.id == nodeList[int(path[i])].id and edge.source.id == nodeList[int(path[ (i+1)%len(path) ])].id:
				edge.isThick = 1
				break


def highlightPath(ifile): #highlight the optimal path. The file conteins an ordered list of nodes to visit
	global nodeList
	path = []
	for i in range(len(nodeList)+1):
		path.append(int(ifile.readline()))
	for i in range(len(nodeList)):
		for edge in nodeList[path[i]].edgeList:
			if edge.dest.id == nodeList[path[i+1]].id:
				edge.isThick = 1
				break


app = QtGui.QApplication(sys.argv)


global counter
counter = 0

global noEdge #set high weight when there is no edge in the graph
noEdge =  32767 

global rightPressed
global firstNode
global secondNode

rightPressed = 0


global scene
scene = GraphScene()

global nodeList
nodeList = []


view = QtGui.QGraphicsView(scene)
main_window = MainWindow()
main_window.setCentralWidget(view)
main_window.show()


sys.exit(app.exec_())


