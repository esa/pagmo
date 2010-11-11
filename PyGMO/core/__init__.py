# -*- coding: utf-8 -*-
from _core import *
import threading as _threading
import signal as _signal
import os as _os

_orig_signal = _signal.getsignal(_signal.SIGINT)
_main_pid = _os.getpid()

# Alternative signal handler which ignores sigint if called from a child process.
def _sigint_handler(signum,frame):
	import os
	if os.getpid() == _main_pid:
		_orig_signal(signum,frame)

_signal.signal(_signal.SIGINT,_sigint_handler)

# Global lock used when starting processes.
_rlock = _threading.RLock()

class base_island(_core._base_island):
	def __init__(self,*args):
		if len(args) == 0:
			raise ValueError("Cannot initialise base island without parameters for the constructor.")
		_core._base_island.__init__(self,*args)
	def get_name(self):
		return str(type(self))

# This is the function that will be called by the separate process
# spawned from py_island.
def _process_target(q,a,p):
	try:
		tmp = a.evolve(p)
		q.put(tmp)
	except BaseException as e:
		print('Exception caught during evolution:')
		print(e)
		q.put(0)

class py_island(base_island):
	from PyGMO import migration as _migr
	def __init__(self,prob, algo, pop = None, n = 0, migr_prob = 1., s_policy = _migr.best_s_policy(), r_policy = _migr.fair_r_policy()):
		if pop is None:
			super(py_island,self).__init__(prob,algo,n,migr_prob,s_policy,r_policy)
		else:
			super(py_island,self).__init__(pop,algo,migr_prob,s_policy,r_policy)
	def __copy__(self):
		retval = py_island(None,self.algorithm,self.population,None,self.migration_probability,self.s_policy,self.r_policy)
		return retval
	def _perform_evolution(self,algo,pop):
		try:
			import multiprocessing as mp
			q = mp.Queue()
			# Apparently creating/starting processes is _not_ thread safe:
			# http://stackoverflow.com/questions/1359795/error-while-using-multiprocessing-module-in-a-python-daemon
			# Protect with a global lock.
			with _rlock:
				process = mp.Process(target = _process_target, args = (q,algo,pop))
				process.start()
			retval = q.get()
			if isinstance(retval,int):
				raise RuntimeError()
			return retval
		except BaseException as e:
			print(e)
			raise RuntimeError()
