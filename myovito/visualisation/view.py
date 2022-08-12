from ovito import dataset
from ovito.io import import_file
import ovito
from ovito.vis import *
# from PyQt5.QtGui import QPainter
# import math
from ovito.vis import Viewport,TachyonRenderer, OpenGLRenderer 
import sys


class Viewer:
	def __init__(self, filename):
		self.filename = filename

		self.pipeline = import_file(filename, multiple_frames=True)


		self.particles = self.pipeline.source.data.particles.vis


 
	def view(self, vtype=Viewport.Type.Ortho, fov=None, camera_dir= None, camera_pos=None, zoom_all=False,frame = 0):
		data = self.pipeline.compute(frame)
		self.pipeline.add_to_scene()

		vp = Viewport(type=vtype)


		if zoom_all:
			vp.zoom_all()

		if fov!=None:
			vp.fov = fov # vertical size of the visible area in units of length
		if camera_dir!=None:
			vp.camera_dir= camera_dir
		if camera_pos!=None:
			vp.camera_pos= camera_pos
		self.vp = vp

	def render(self,size=(800,400), renderer='opengl'):
		renderer = OpenGLRenderer ()

		if renderer =='opengl':
			renderer = OpenGLRenderer()
		self.vp.render_image(size=size,filename=self.filename.split('/')[-1]+".png", background=(1,1,1), renderer=renderer)
		


