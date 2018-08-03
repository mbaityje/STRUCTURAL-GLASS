#!/usr/bin/env python
#
# Module for time lists
# 
from __future__ import print_function
import numpy as np

def trace(frame, event, arg):
	'''
	Function to trace at what line things (e.g. segmentation fault) happen.
	Before launching the routine I want to test, execute the command
	
	sys.settrace(trace)
	'''
	print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
	return trace
