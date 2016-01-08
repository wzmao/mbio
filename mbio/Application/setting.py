# -*- coding: utf-8 -*-
"""This module contains some setting functions.
"""

__author__ = 'Wenzhi Mao'

__all__ = ['getMatplotlibDisplay']

def getMatplotlibDisplay(**kwargs):
	"""Check if the current job is under no display status."""
	import os

	return (os.getenv('DISPLAY') is not None)
	# pid = str(os.getpid())
	# findssh=False
	# while pid!="0":
	# 	f = open('/proc/' + str(pid) + '/cmdline', 'r')
	# 	temp = f.read().replace('\0', ' ')
	# 	f.close()
	# 	if temp.find("ssh")!=-1:
	# 		findssh=True
	# 		break
	# 	pid=os.popen("ps -p {0} -oppid=".format(pid)).read().strip()
	# return findssh