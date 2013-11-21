# Defines common behavior of object in nampy.core
# This script is developed from Object.py
# in COBRApy, which was distributed
# under GNU GENERAL PUBLIC LICENSE Version 3

# Ebrahim, A., Lerman, J. a, Palsson, B. O., & Hyduke, D. R. (2013). 
# COBRApy: COnstraints-Based Reconstruction and Analysis for Python. 
# BMC systems biology, 7(1), 74. doi:10.1186/1752-0509-7-74

# Could replace this with collections.OrderedDict
# as ebrahim et al. suggest 
# but there are some nice default methods defined here
# like _check

# Could rewrite this but there are some
# nice default methods defined here
# like _check

class Object(object):
    #__slots__ = ['id']
    def __init__(self, id=None):
        """
        id: None or a string
        
        """
        self.id = id
        #The following two fields will eventually
        #be objects that enforce basic rules about
        #formatting notes and annotation
        self.notes = {}
        self.annotation = {}

    def __getstate__(self):
        """To prevent excessive replication during deepcopy.
        """
        state = self.__dict__.copy()
        if '_model' in state:
            state['_model'] = None
        return state
    
    def guided_copy(self):
        """Trying to make a faster copy procedure for cases where large
        numbers of metabolites might be copied.  Such as when copying reactions.

        This function allows us to manipulate how specific attributes are copied.

        """
        the_copy = self.__class__(self.id)
        [setattr(the_copy, k, v)
         for k, v in self.__dict__.iteritems()]
        return(the_copy)

    ## def __setstate__(self, state):
    ##     self.__dict__.update(state)
    ## def __getstate__(self):
    ##     the_dict = dict([(x, eval('self.%s'%x))
    ##                      for x in self.__slots__])
    ##     return(the_dict)
    ## def __setstate__(self, the_dict):
    ##     self.id = the_dict.pop('id')
    #Allows comparison of Objects based on ids and with ids
    #
    #Not the best idea.  This will be removed in the next major
    #release
    #

    def __lt__(self, other):
        if hasattr(other, 'id'):
            x = self.id < other.id
        elif type(other) == type(self.id):
            x = self.id < other
        return x
    
    def __le__(self, other):
        if hasattr(other, 'id'):
            x = self.id <= other.id
        elif type(other) == type(self.id):
            x = self.id <= other
        return x
    
    def __gt__(self, other):
        if hasattr(other, 'id'):
            x = self.id > other.id
        elif type(other) == type(self.id):
            x = self.id > other
        return x
    
    def __ge__(self, other):
        if hasattr(other, 'id'):
            x = self.id >= other.id
        elif type(other) == type(self.id):
            x = self.id >= other
        return x
    
    def __ne__(self, other):
        x = True
        if hasattr(other, 'id'):
            x = self.id != other.id
        elif type(other) == type(self.id):
            x = self.id != other
        return x
    
    def __eq__(self, other):
        x = False
        if hasattr(other, 'id'):
            x = self.id == other.id
        elif type(other) == type(self.id):
            x = self.id == other
        return x

    def startswith(self, x):
        return self.id.startswith(x)

    def endswith(self, x):
        return self.id.endswith(x)

    def __contains__(self, x):
        return self.id.__contains__(x)

    def __iter__(self):
       return list(self.id).__iter__()

    def __getitem__(self, index):
        return self.id[index]
    
    def __getslice__(self, i, j):
        return self.id[i:j]
    
    def __repr__(self):
        return repr(self.id)

    def __str__(self):
        return str(self.id)
