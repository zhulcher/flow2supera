import numpy as np
import os

with np.load(os.path.join(os.path.dirname(__file__),'pdg_data/pdg.npz'),'r') as f:
    _PDG_DATA = dict(f)

def pdg2mass(pdg_code):
    '''
    Given a PDG code, return the mass in MeV.
    If the PDG code is 10 digits = a nucleus, it returns atomic number x 1000.
    '''
    if pdg_code > 1000000000:
        return int(str(pdg_code)[-4:-1])*1000.

    where = np.where(_PDG_DATA['pdg_code']==pdg_code)[0]
    if len(where)<1:
        return -1
    assert len(where) == 1
    return _PDG_DATA['mass'][where[0]]*1000.


'''
import ROOT
import numpy as np
pdg = ROOT.TDatabasePDG()
pdg.ReadPDGTable()

data=dict(pdg_code=[],mass=[])
for p in pdg.ParticleList():
    data['pdg_code'].append(p.PdgCode())
    data['mass'].append(p.Mass())
np.savez('pdg.npz',**data)
'''