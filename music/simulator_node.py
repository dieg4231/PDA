#!/usr/bin/env python

from MobileRobotSimulator import *
import time


gui=MobileRobotSimulator()


if __name__ == "__main__":
	print("sss")
	val = []
	time.sleep(1)
	#print("ddd")
	i =0.0
	while True:
		i = sys.stdin.readline()

		if str(i) == "-----\n":
			print("xxx")
			for i in val:
				if i > -100:
					gui.handle_service(val)
					break
			val = []
		elif str(i) == "*\n":
			gui.handle_pausa()
		else:
			val.append(float(i))

		#print(i)
		
		
