# -*- coding: utf-8 -*-
from _core import *
import threading as _threading
import signal as _signal
import os as _os

__doc__ = 'PyGMO core module.'
__all__ = ['archipelago','base_island','champion','distribution_type','individual','ipy_island','island','local_island','migration_direction','population','py_island']

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

# Raw C++ base island class.
_base_island = _core._base_island

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

def _generic_island_ctor(self,*args,**kwargs):
	"""Unnamed arguments:
	
		#. algorithm
		#. problem or population
		#. number of individuals (optional and valid only if the second argument is a problem, defaults to 0 if not specified)
	
	Keyword arguments:
	
		* *migr_prob* -- migration probability (defaults to 1)
		* *s_policy* -- migration selection policy (defaults to 'best selection' policy)
		* *r_policy* -- migration replacement policy (defaults to 'fair replacement' policy)
	
	"""
	from PyGMO.algorithm import _base as _base_algorithm, base as base_algorithm
	from PyGMO.problem import _base as _base_problem, base as base_problem
	from PyGMO.migration import best_s_policy, fair_r_policy, _base_s_policy, _base_r_policy
	if len(args) < 2 or len(args) > 3:
		raise ValueError("Unnamed arguments list must have either 2 or three elements, but %d elements were found instead." % (len(args),))
	if not isinstance(args[0],_base_algorithm):
		raise TypeError("The first unnamed argument must be an algorithm.")
	ctor_args = [args[0]]
	if isinstance(args[1],_base_problem):
		ctor_args.append(args[1])
		if len(args) == 3:
			if not isinstance(args[2],int):
				raise TypeError("Please provide an integer for the number of individuals in the island.")
			ctor_args.append(args[2])
		else:
			ctor_args.append(0)
	elif isinstance(args[1],population):
		if len(args) == 3:
			raise ValueError("When the second unnamed argument is a population, there cannot be a third unnamed argument.")
		ctor_args.append(args[1])
	else:
		raise TypeError("The second unnamed argument must be either a problem or a population.")
	ctor_args.append(kwargs.pop('migr_prob',1.))
	if not isinstance(ctor_args[-1],float):
		raise TypeError("Migration probability must be a float.")
	ctor_args.append(kwargs.pop('s_policy',best_s_policy()))
	if not isinstance(ctor_args[-1],_base_s_policy):
		raise TypeError("s_policy must be a migration selection policy.")
	ctor_args.append(kwargs.pop('r_policy',fair_r_policy()))
	if not isinstance(ctor_args[-1],_base_r_policy):
		raise TypeError("r_policy must be a migration replacement policy.")
	if isinstance(self,base_island):
		super(type(self),self).__init__(*ctor_args)
	elif isinstance(self,_base_island):
		self.__original_init__(*ctor_args)
	else:
		assert(self is None)
		n_pythonic_items = 0
		if isinstance(args[0],base_algorithm):
			n_pythonic_items += 1
		if isinstance(args[1],base_problem):
			n_pythonic_items += 1
		elif isinstance(args[1],population) and isinstance(args[1].problem,base_problem):
			n_pythonic_items += 1
		if n_pythonic_items > 0:
			return py_island(*args,**kwargs)
		else:
			return local_island(*args,**kwargs)

local_island.__original_init__ = local_island.__init__
local_island.__init__ = _generic_island_ctor

# This is the function that will be called by the separate process
# spawned from py_island.
def _process_target(q,a,p):
	try:
		tmp = a.evolve(p)
		q.put(tmp)
	except BaseException as e:
		q.put(e)

class py_island(base_island):
	"""Python island.
	
	This island will launch evolutions using the multiprocessing module, available since Python 2.6.
	Each evolution is transparently dispatched to a Python interpreter in a separate process.

	"""
	__init__ = _generic_island_ctor
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
		return "Python multiprocessing island"

# This is the function that will be called by the task client
# in ipy_island.
def _maptask_target(a,p):
	try:
		return a.evolve(p)
	except BaseException as e:
		return e

class ipy_island(base_island):
	"""Parallel IPython island.
	
	This island will launch evolutions using IPython's MapTask interface. The evolution will be dispatched
	to IPython engines that, depending on the configuration of IPython/ipcluster, can reside either on the
	local machine or on other remote machines.
	
	See: http://ipython.scipy.org/doc/stable/html/parallel/index.html
	
	"""
	# NOTE: when using an IPython island, on quitting IPython there might be a warning message
	# reporting an exception being ignored. This seems to be a problem in the foolscap library:
	# http://foolscap.lothar.com/trac/ticket/147
	# Hopefully it will be fixed in the next versions of the library.
	__init__ = _generic_island_ctor
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

def island(*args,**kwargs):
	return _generic_island_ctor(None,*args,**kwargs)

island.__doc__ = '\n'.join(['Island factory function.\n\nThis function will return an instance of an island object\nbuilt according to the following rule: '+
	'if the arguments include\neither a pythonic problem or a pythonic algorithm, then an instance\nof :class:`py_island` will be returned; '+
	'otherwise, an instance of\n:class:`local_island` will be returned.'] + [s.replace('\t','') for s in _generic_island_ctor.__doc__.split('\n')])

del s

def _get_island_list():
	from PyGMO import core
	names = filter(lambda n: not n.startswith('_') and not n.startswith('base') and n.endswith('_island'),dir(core))
	try:
		from IPython.kernel.client import TaskClient, MapTask
	except ImportError:
		names = filter(lambda n: n != 'ipy_island',names)
	return [core.__dict__[n] for n in names]
