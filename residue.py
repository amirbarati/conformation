import pandas as pd

class Residue(object):
  def __init__(self, resSeq, chain_id=None, chain_name=None, ballosteros_weinstein=None):
    self.resSeq = resSeq
    self.chain_id = chain_id
    self.chain_name = chain_name
    self.ballosteros_weinstein = ballosteros_weinstein

  def conformation_res_mdtraj_res_equivalent(self, mdtraj_res):
    if self.resSeq != mdtraj_res.resSeq:
      return False

    if self.chain_id is not None and (self.chain_id != mdtraj_res.chain_id):
      return False

    if self.chain_name is not None and (self.chain_name != mdtraj_res.chain_name):
      return False 

    return True

class ContactFeature(object):
  def __init__(self, residue_i, residue_j):
    self.residue_i = residue_i
    self.residue_j = residue_j