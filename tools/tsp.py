#!/usr/bin/python

import sys
import math
from PyQt4 import Qt,QtGui,QtCore


class Node(QtGui.QGraphicsItem):
	def __init__(self, scn, parent = None):
		QtGui.QGraphicsItem.__init__(self,parent)
		self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
		self.setFlag(QtGui.QGraphicsItem.ItemSendsGeometryChanges)
    		self.setCacheMode(QtGui.QGraphicsItem.DeviceCoordinateCache)
     		self.setZValue(-1)
		self.edgeList = []
		global counter
		self.id = counter
		counter = counter + 1
		self.scene = scn
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
		if event.button() == QtCore.Qt.RightButton:
			global rightPressed
			if rightPressed:
				if not(firstNode.id == self.id): #avoid to create edge from a node to itself
					global edgeList
					secondNode = self
					newEdge = Edge(firstNode,secondNode)
					self.addEdge(newEdge)
					self.scene.addItem(newEdge)
					edgeList.append(newEdge)
					rightPressed = 0
			else:
				global firstNode
				firstNode = self
				rightPressed = 1

class Edge(QtGui.QGraphicsItem):
	def __init__(self, sourceNode, destNode):
		QtGui.QGraphicsItem.__init__(self)
		self.source = sourceNode
		self.dest = destNode
		sourceNode.addEdge(self)
		destNode.addEdge(self)
		self.adjust()
		self.arrowsize = 10
		self.weight = 0
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
	def shape(self):
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
		dialog = QtGui.QInputDialog()
		self.weight = float(dialog.getText(widget,"Weight", "Edge's weight")[0])


class GraphScene(QtGui.QGraphicsScene):
	def __init__(self):
		QtGui.QGraphicsScene.__init__(self)
	def mouseDoubleClickEvent(self,mouseEvent):
		if mouseEvent.button() == QtCore.Qt.LeftButton:
			global nodeList
			newNode = Node(self)
			self.addItem(newNode)
			newNode.setPos(mouseEvent.scenePos())
			nodeList.append(newNode)
	def keyPressEvent(self,keyEvent):
		if (keyEvent.key() == QtCore.Qt.Key_S):
			print generateAdjMatrix(nodeList,edgeList)

def generateAdjMatrix(nodeList, edgeList):
	length = len(nodeList)
	matrix = [[0 for col in range(length)] for row in range(length)]
	for edge in edgeList:
		matrix[edge.source.id][edge.dest.id] = edge.weight
	return matrix


app = QtGui.QApplication(sys.argv)


global counter
counter = 0

global rightPressed
global firstNode
global secondNode

rightPressed = 0


global scene
scene = GraphScene()

global nodeList
nodeList = []

global edgeList
edgeList = []

generateAdjMatrix(nodeList,edgeList)

view = QtGui.QGraphicsView(scene)
view.show()

sys.exit(app.exec_())


