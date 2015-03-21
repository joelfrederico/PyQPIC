import unittest
import PyQPIC as _pq
import numpy as _np
import hashlib
import pkg_resources as pkg


class TestTemplateGeneration(unittest.TestCase):
    def setUp(self):
        pass

    def test_result(self):
        qpic    = _pq.QuickPICSettings()
        magic   = _pq.MagicSettings()
        bunches = _np.array(
            [
                _pq.BunchSettings(beamtype='drive'), _pq.BunchSettings(beamtype='witness')
            ]).flatten()
        plasma  = _pq.PlasmaSettings(qpic=qpic, magic=magic, bunches=bunches)
        box     = _pq.BoxSettings()
        
        _pq.deckgen(
            plasma  = plasma,
            bunches = bunches,
            box     = box,
            magic   = magic,
            qpic    = qpic
            )
        
        md5_new = hashlib.md5()
        with open('rpinput', 'rb') as fid:
            md5_new.update(fid.read())
        
        hash_new = md5_new.hexdigest()
        
        md5_orig = hashlib.md5()
        with pkg.resource_stream('PyQPIC', 'tests/files/rpinput') as fid:
            md5_orig.update(fid.read())
        
        hash_orig = md5_orig.hexdigest()
        self.assertEqual(hash_orig, hash_new)
