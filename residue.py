import pandas as pd

class Residue(object):
  def __init__(self, resSeq, chain_id=None, chain_name=None, ballosteros_weinstein=None):
    self.resSeq = resSeq
    self.chain_id = chain_id
    self.chain_name = chain_name
    self.ballosteros_weinstein = ballosteros_weinstein

class ContactFeature(object):
  def __init__(self, residue_i, residue_j):
    self.residue_i = residue_i
    self.residue_j = residue_j