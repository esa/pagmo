# -*- coding: utf-8 -*-
from _core import *
import multiprocessing as _mp

class base_island(_core._base_island):
	def __init__(self,*args):
		if len(args) == 0:
			raise ValueError("Cannot initialise base island without parameters for the constructor.")
		_core._base_island.__init__(self,*args)
	def get_name(self):
		return str(type(self))

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
	@staticmethod
	def _process_target(conn,a,p):
		try:
			tmp = a.evolve(p)
		except Exception as e:
			conn.send(e)
			conn.close()
			return
		conn.send(tmp)
		conn.close()
	def _perform_evolution(self,algo,pop):
		print('here')
		return algo.evolve(pop)
		# -----------
		try:
			parent_conn, child_conn = _mp.Pipe()
			process = _mp.Process(target = py_island._process_target, args = (child_conn,algo,pop))
			process.start()
			retval = parent_conn.recv()
			# After receiving the data wait for the process to exit cleanly.
			process.join()
			if isinstance(retval,Exception):
				raise retval
			elif isinstance(retval,type(pop)):
				return retval
			else:
				raise RuntimeError('Slave process returned unknown type: %s' % str(type(retval)))
		except Exception as e:
			print('Island\'s evolution method raised an exception:')
			print(e)
			raise RuntimeError('Evolution interrupted')
