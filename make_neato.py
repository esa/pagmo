def make_neato(arch):
	t_str = arch.topology.__repr__()
	t_dict = {}
	for i in t_str.split():
		tmp = i.split('->')
		tmp_int = int(tmp.pop(0))
		t_dict[tmp_int] = [int(j) for j in tmp]
	max_n_edges = max([len(t_dict[i]) for i in t_dict])
	min_n_edges = min([len(t_dict[i]) for i in t_dict])
	print max_n_edges
	print min_n_edges
	retval = 'strict digraph foo {\n'
	retval += 'edge [len=' + str(max_n_edges) + '];\n'
	for i in t_dict:
		retval += str(i) + ' [height=' + str(len(t_dict[i])/5.) + ',width=' + str(len(t_dict[i])/5.) + ',label=' + str(len(t_dict[i])) + ',fontsize=' + \
			str(len(t_dict[i])*10) + '];\n'
		#retval += 'edge [w=' + str(len(t_dict[i])*10./max_n_edges) + '];\n'
		for j in t_dict[i]:
			retval += str(i) + '->' + str(j) + ';\n'
	retval += '}'
	return retval
