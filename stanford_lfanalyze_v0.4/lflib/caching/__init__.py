import os, sys
import cPickle as pickle

class Cache(object):
    """
    Abstract base class for cache levels.
    """
    serialize_values = True # flip this in subclasses to override serialization

    def write(self, key, value):
        "Override to write the already-serialized valuye to your cache"
        raise NotImplementedError

    def read(self, key):
        "Override to read the serialized value from your cache and return it"
        raise NotImplementedError

    def __delitem__(self, key):
        "Override this to delete an item from the cache"
        raise NotImplementedError

    def __contains__(self, key):
        "You can override this if there's a cheap way to check for the existence of a key in your cache"
        raise NotImplementedError

    def serialize(self, value):
        return pickle.dumps(value)

    def deserialize(self, value):
        return pickle.loads(value)

    def __setitem__(self, key, value):
        if self.serialize_values:
            value = self.serialize(value)
        return self.write(key, value)

    def __getitem__(self, key):
        value = self.read(key)
        if self.serialize_values:
            value = self.deserialize(value)
        return value

class MultiCache(object):
    """
    Set and get from multiple cache levels in sequence.
    Initialize this with a list of Cache instances, in the order you want them checked.
    Reads will return the first cache result found.
    Writes will write to all caches in order.  Some writes may be made asyncronously.
    """
    def __init__(self, *cache_levels):
        self.cache_levels = cache_levels

    def __setitem__(self, key, value):
        for cl in self.cache_levels:
            cl[key] = value

    def __getitem__(self, key):
        for cl in self.cache_levels:
            try:
                return cl[key]
            except KeyError:
                continue
        else:
            raise KeyError

    def __delitem__(self, key):
        """
        Delete the key from all caches, ignoring KeyErrors.
        """
        for cache in self.cache_levels:
            try:
                del cache[key]
            except KeyError:
                pass

class FileSystemCache(Cache):
    """
    Dead-simple filesystem cache
    
    >>>
    """
    def __init__(self, root_dir):
        self.root = root_dir
        if not os.path.exists(self.root):
            os.makedirs(self.root)

    def __repr__(self):
        return "<FileSystemCache at %s>" % self.root

    def __contains__(self, key):
        return os.path.exists( os.path.join( self.root, key) )

    def __delitem__(self, key):
        if key not in self:
            raise KeyError
        path = os.path.join(self.root, key)
        os.remove(path)
        try:
            os.removedirs(os.path.split(path)[0])
        except OSError:
            # leaf directory is non-empty.  ignore.
            pass

    def write(self, key, value):
        head, tail = os.path.split(key)
        head = os.path.join(self.root, head)
        path = os.path.join(head, tail)
        if os.path.isdir(path):
            raise Exception("Cannot write to file because it is a directory: %s" % path)
        if not os.path.exists(head):
            os.makedirs(head)
        if not os.path.isdir(head):
            raise IOError("Key Rejected by %s because %s is not a directory" % (self, head) )
        with open(path, 'w') as outfile:
            outfile.write(value)

    def read(self, key):
        path = os.path.join(self.root, key)
        if not os.path.exists(path):
            raise KeyError
        if os.path.isdir(path):
            raise Exception("Cannot read from file because it is a directory: %s" % path)
        with open(path) as infile:
            value = infile.read()
        return value

class DictCache(dict):
    "Dead-simple in-memory python cache"
    pass

