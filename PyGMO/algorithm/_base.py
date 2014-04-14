# -*- coding: utf-8 -*-
from PyGMO.algorithm._algorithm import _base


class base(_base):

    """
    All Algorithms written in Python derive from this class
    """

    def __init__(self):
        _base.__init__(self)

    def get_name(self):
        return str(type(self))

    def __get_deepcopy__(self):
        from copy import deepcopy
        return deepcopy(self)
