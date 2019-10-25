#!/usr/bin/env python

from MobileRobotSimulator import *
import time


gui=MobileRobotSimulator()


if __name__ == "__main__":
	print("sss")
	time.sleep(1)
	#print("ddd")
	i =0.0
	while True:
		i = float(sys.stdin.readline())
		#print(i)
		print("xxx")
		gui.handle_service(i)
		
		
