#!/usr/bin/env python
import signal
import sys
def signal_handler(sig, frame):
        print('You pressed Ctrl+C!')
        sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)
value = 0

while(True):
    value = sys.stdin.readline()
    print(value)
