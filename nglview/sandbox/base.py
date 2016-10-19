import abc, six

@six.add_metaclass(abc.ABCMeta)
class BaseMD(object):
    '''Unstable API
    '''
    @abc.abstractmethod
    def initialize(self):
        pass

    @abc.abstractmethod
    def update(self):
        pass

    @abc.abstractmethod
    def stop(self):
        pass

