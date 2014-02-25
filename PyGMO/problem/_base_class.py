from PyGMO.problem._problem import _base


class base(_base):

    def __init__(self, *args):
        if len(args) == 0:
            raise ValueError(
                "Cannot initialise base problem without parameters for the constructor.")
        _base.__init__(self, *args)

    def _get_typename(self):
        return str(type(self))

    def __get_deepcopy__(self):
        from copy import deepcopy
        return deepcopy(self)

    def get_name(self):
        return self._get_typename()
