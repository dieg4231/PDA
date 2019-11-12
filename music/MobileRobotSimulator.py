from Tkinter import *
from tkFont import Font
import threading
import ttk
import time
import math
from PIL import Image
from PIL import ImageDraw
import tkMessageBox
import os
import signal
import math
import numpy as np


 
class MobileRobotSimulator(threading.Thread):

	def __init__(self):
		threading.Thread.__init__(self)

		self.spectrum =  np.zeros(360)
		self.angles =360;
		self.stopped = False 
		# map size in meters
		self.mapX = 1 
		self.mapY = 1
		# canvas size in pixels
		self.canvasX= (self.angles*3)+20
		self.canvasY= 600
		# robot position and angle
		self.robotAngle=0
		self.robotX=-100
		self.robotY=-100

		self.p_giro=0
		self.p_distance=0
		
		self.polygonMap = []
		self.nodes_image = None
		self.light=-1
		self.robot=-1

		self.flagOnce=False

		self.light_x = 0
		self.light_y = 0
		self.startFlag = False

		self.lasers = []
		self.sensors_value = [512];
		self.sensors_values = [512];
		self.sensors_values_aux = [512];
		self.sensors_values_aux_old = [512];
		for i in range(512):
			self.sensors_value.append(0)
			self.sensors_values.append(0)
			self.sensors_values_aux.append(0)


		self.graph_list = [200]
		for i in range(200):
			self.graph_list.append(0)

		self.rewind=[]
		self.trace_route= []
		self.varShowNodes   = False
		self.grid =[]
		self.contador = 0;
		self.contador_ = 0;
		self.bandera = True

		self.X_pose = 0
		self.Y_pose = 0

		self.num_polygons = 0  #How many polygons exist in the field.
		self.polygons = []	   #Stors polygons vertexes
		self.polygons_mm = []  #Stors 2 vertexses  for each polygon the maximum and minimum  x and y  points.

		self.objects_data = []
		self.grasp_id = False
		self.current_object = -1;
		self.current_object_name = -1;

		self.initX = 0
		self.initY = 0
		self.initR = 0

		self.ang_arr = 0

		self.b = 999
		self.c = 999
		self.d = 999

		self.bars = []

		self.start()

	def kill(self):  # When press (x) window
		print("Bye")
		self.root.quit()

	def pause(self,*args):
		self.w.itemconfig(self.pause_line, fill='#EE0000')
		for i in self.bars:
			self.w.delete(i)
		self.bars = []
		self.w.update() 

	def new_angle(self,*args):
		maxx = 0;

		for i in range(1,self.angles):
			if self.spectrum[maxx] < self.spectrum[i]:
				maxx=i

		#print(maxx)
		for i in self.bars:
			self.w.delete(i)
		self.bars = []
#		300 = self.spectrum[maxx]	
		cta=0
		for i in xrange(0,(self.angles*3),3):
			
			try:
				val =(self.spectrum[cta]*300)/self.spectrum[maxx]
			except:
				val = 0;
				print("no mms")


			self.bars.append(self.w.create_line(i+10,self.canvasY/2-2 ,i +10,self.canvasY/2-val  ,fill = "#0000FF")) 
			cta = cta +1

		self.w.itemconfig(self.pause_line, fill='#00EE00')
		self.w.update() 

  			

	def handle_service(self,spectrum):
		self.spectrum = ( np.array(spectrum)*.8+(np.array(self.spectrum)) )
		self.a.set(1)

	def handle_pausa(self):
		self.b.set(1)

	def gui_init(self):

		self.backgroundColor = '#EDEDED';#"#FCFCFC";
		self.entrybackgroudColor = "#FBFBFB";##1A3A6D";
		self.entryforegroundColor = '#37363A';
		self.titlesColor = "#303133"
		self.menuColor = "#ECECEC"
		self.menuButonColor = "#375ACC"
		self.menuButonFontColor = "#FFFFFF"
		self.obstacleInnerColor = '#447CFF'
		self.obstaclesOutlineColor="#216E7D"#'#002B7A'
		self.buttonColor = "#1373E6"
		self.buttonFontColor = "#FFFFFF"
		self.canvasColor = "#FFFFFF"
		self.gridColor = "#D1D2D4"
		self.wheelColor  = '#404000'  
		self.robotColor  = '#F7CE3F'  
		self.hokuyoColor = '#4F58DB' 
		self.arrowColor  = '#1AAB4A' 
		self.laserColor  = "#00DD41" 

		self.root = Tk()
		self.root.protocol("WM_DELETE_WINDOW", self.kill)
		self.root.title("Mobile Robot Simulator")

		self.barMenu = Menu(self.root)
		self.settingsMenu = Menu(self.barMenu, tearoff=0)
		self.submenuTheme = Menu(self.settingsMenu, tearoff=0)
		self.submenuCanvas = Menu(self.settingsMenu, tearoff=0)
		self.root.config(menu=self.barMenu )

		self.content   = Frame(self.root)
		self.frame     = Frame(self.content,borderwidth = 5, relief = "flat", width = 600, height = 900 ,background = self.backgroundColor)
		self.rightMenu = Frame(self.content,borderwidth = 5, relief = "flat", width = 300, height = 900 ,background = self.backgroundColor)
		self.w = Canvas(self.frame, width = self.canvasX, height = self.canvasY, bg=self.canvasColor)
		self.w.pack()
		
		self.content  .grid(column = 0 ,row = 0 ,sticky = (N, S, E, W))
		self.frame    .grid(column = 0 ,row = 0 ,columnspan = 3 ,rowspan = 2 ,sticky = (N, S, E, W))
		
		self.root.columnconfigure(0, weight=1)
		self.root.rowconfigure(0, weight=1)
		self.content.columnconfigure(0, weight = 3)
		self.content.columnconfigure(1, weight = 3)
		self.content.columnconfigure(2, weight = 3)
		self.content.columnconfigure(3, weight = 1)
		self.content.columnconfigure(4, weight = 1)
		self.content.rowconfigure(1, weight = 1)

		self.scale = self.w.create_line(10,self.canvasY/2 ,self.canvasX-10 ,self.canvasY/2  ,fill = "#0000FF") 
		
		for i in xrange(0,(self.angles*3),18):
			self.w.create_line(i+10,self.canvasY/2-2 ,i +10,self.canvasY/2+2  ,fill = "#0000FF") 
			self.w.create_text(i+10,self.canvasY/2 +30,fill="darkblue",font="Times 12 italic ",angle=90,text=str( (i/3)-(self.angles/2) ) )
		radio=20
		self.pause_line = self.w.create_oval(100-radio,self.canvasY-100-radio, 100+radio,self.canvasY-100+radio   , outline="#FFFFFF", fill="#111111", width=1)

		self.a = IntVar(value=3)
		self.a.trace("w", self.new_angle)

		self.b = IntVar(value=3)
		self.b.trace("w", self.pause)


	def run(self):	
		self.gui_init()
		self.root.mainloop()



def signal_handler(sig, frame):
        print('You pressed Ctrl+C!')
        sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)


#m = MobileRobotSimulator()
#m.run()

