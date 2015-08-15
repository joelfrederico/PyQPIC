import pkg_resources as _pkg
import shutil as _shutil
import os as _os


def get_bin(path=None):
    with _pkg.resource_stream('PyQPIC', 'bin/qpic.e') as fid_read:
        if path is None:
            path = _os.getcwd()
        with open(path) as fid_write:
            _shutil.copyfileobj(fid_read, fid_write)
