
from time import time
from functools import wraps
import os
import shutil
import datetime
import getpass
import colored

class LazyLogger(object):
    def __init__(self, name=''):
        self._name = name
        pass

    def debug(self, msg):
        print(("{0} - " + colored.fg("blue") +
               "[Debug]: " + colored.fg("black")).format(self._name) , msg)

    def notice(self, msg):
        print(("{0} - " + colored.fg("green") +
               "[Info]: " + colored.fg("black")).format(self._name) , msg)

    def info(self, msg):
        self.notice(msg)

    def warning(self, msg):
        print(("{0} - " + colored.fg("orange_1") +
               "[Warning]: " + colored.fg("black")).format(self._name) , msg)

    def warn(self, msg):
        print(("{0} - " + colored.fg("orange_1") +
               "[Warning]: " + colored.fg("black")).format(self._name) , msg)

    def error(self, msg):
        print(("{0} - " + colored.fg("red") +
               "[Error]: " + colored.fg("black")).format(self._name) , msg)


def timeit(func):
    @wraps(func)
    def _time_it(*args, **kwargs):
        start = int(round(time() * 1000))
        try:
            return func(*args, **kwargs)
        finally:
            end_ = int(round(time() * 1000)) - start
            print(f"Total execution time: {end_ if end_ > 0 else 0} ms")
    return _time_it


def build_directory_name(config_file, basedir):
    rawfilename = os.path.splitext(os.path.basename(config_file))[0]
    bdir = basedir
    now = datetime.datetime.now()
    directory_name =  getpass.getuser()
    directory_name += '_' + rawfilename
    directory_name += '_' + now.strftime("%Y-%m-%dT%H:%M:%S")
    directory_name += '/'
    directory_name = os.path.join(bdir,directory_name)
    return directory_name

def copy_cfgFileToDir(directory, config_file):
    destination_file = directory +  'config/' + os.path.basename(config_file)
    os.mkdir(directory + 'config/')
    shutil.copyfile(config_file, destination_file)
