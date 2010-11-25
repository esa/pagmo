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
_process_lock = _threading.Lock()

class base_island(_core._base_island):
	def __init__(self,*args):
		if len(args) == 0:
			raise ValueError("Cannot initialise base island without parameters for the constructor.")
		_core._base_island.__init__(self,*args)
	def get_name(self):
		return str(type(self))
	def __get_deepcopy__(self):
		from copy import deepcopy
		return deepcopy(self)

# This is the function that will be called by the separate process
# spawned from py_island.
def _process_target(q,a,p):
	try:
		tmp = a.evolve(p)
		q.put(tmp)
	except BaseException as e:
		q.put(e)

class py_island(base_island):
	from PyGMO import migration as _migr
	def __init__(self,prob, algo, pop = None, n = 0, migr_prob = 1., s_policy = _migr.best_s_policy(), r_policy = _migr.fair_r_policy()):
		if pop is None:
			super(py_island,self).__init__(prob,algo,n,migr_prob,s_policy,r_policy)
		else:
			super(py_island,self).__init__(pop,algo,migr_prob,s_policy,r_policy)
	def _perform_evolution(self,algo,pop):
		try:
			import multiprocessing as mp
			q = mp.Queue()
			# Apparently creating/starting processes is _not_ thread safe:
			# http://bugs.python.org/issue1731717
			# http://stackoverflow.com/questions/1359795/error-while-using-multiprocessing-module-in-a-python-daemon
			# Protect with a global lock.
			with _process_lock:
				process = mp.Process(target = _process_target, args = (q,algo,pop))
				process.start()
			retval = q.get()
			with _process_lock:
				process.join()
			if isinstance(retval,BaseException):
				raise retval
			return retval
		except BaseException as e:
			print('Exception caught during evolution:')
			print(e)
			raise RuntimeError()
	def get_name(self):
		return "Python island"

# This is the function that will be called by the task client
# in ipy_island.
def _maptask_target(a,p):
	try:
		return a.evolve(p)
	except BaseException as e:
		return e

class ipy_island(base_island):
	# NOTE: when using ipython island, on quitting IPython there might be a warning message
	# reporting an exception being ignored. This seems to be a problem in the foolscap library:
	# http://foolscap.lothar.com/trac/ticket/147
	# Hopefully it will be fixed in the next versions of the library.
	from PyGMO import migration as _migr
	def __init__(self,prob, algo, pop = None, n = 0, migr_prob = 1., s_policy = _migr.best_s_policy(), r_policy = _migr.fair_r_policy()):
		if pop is None:
			super(ipy_island,self).__init__(prob,algo,n,migr_prob,s_policy,r_policy)
		else:
			super(ipy_island,self).__init__(pop,algo,migr_prob,s_policy,r_policy)
	def _perform_evolution(self,algo,pop):
		try:
			from IPython.kernel.client import TaskClient, MapTask
			# Create task client.
			tc = TaskClient()
			# Create the task.
			mt = MapTask(_maptask_target,args = (algo,pop))
			# Run the task.
			task_id = tc.run(mt)
			# Get retval.
			retval = tc.get_task_result(task_id,block = True)
			if isinstance(retval,BaseException):
				raise retval
			return retval
		except BaseException as e:
			print('Exception caught during evolution:')
			print(e)
			raise RuntimeError()
	def get_name(self):
		return "Parallel IPython island"
