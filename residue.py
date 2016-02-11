import pandas as pd

class Residue(object):
  def __init__(self, resSeq, chain_id=None, chain_name=None, res_name=None, ballosteros_weinstein=None):
    self.resSeq = resSeq
    self.chain_id = chain_id
    self.chain_name = chain_name
    self.ballosteros_weinstein = ballosteros_weinstein
    self.res_name = res_name

  def __repr__(self):
    if self.chain_id is not None:
      return "%s%d" % (self.chain_id, self.resSeq)
    else:
      return str(self.resSeq)

  def __eq__(self, other): 
    return self.__dict__ == other.__dict__

  def __hash__(self):
    return hash(self.__repr__())

  def __ne__(self, other):
    return (not self.__eq__(other))

  def is_mdtraj_res_equivalent(self, mdtraj_res):
    if self.resSeq != mdtraj_res.resSeq:
      return False

    if self.chain_id is not None and (self.chain_id != mdtraj_res.chain.id):
      return False

    if self.chain_name is not None and (self.chain_name != mdtraj_res.chain.id):
      return False 

    return True

class Atom(object):
  def __init__(self, atom_id=None, resSeq=None, chain_id=None, atom_name=None, res_name=None):
    self.resSeq = resSeq
    self.chain_id = chain_id
    self.atom_name = atom_name
    self.atom_id = atom_id
    self.res_name = res_name

  def __repr__(self):
    return "%s-%s" %(self.res_name, self.atom_name)

  def __eq__(self, other):
    return self.__dict__ == other.__dict__

  def __hash__(self):
    return hash(self.__repr__())

  def __ne__(self, other):
    return (not self.__eq__(other))

  def is_mdtraj_atom_equivalent(self, mdtraj_atom):
    if self.resSeq != mdtraj_atom.residue.resSeq:
      return False
    if self.chain_id is not None and (self.chain_id != mdtraj_atom.residue.chain.id):
      return False
    if self.atom_name is not None and (self.atom_name != mdtraj_atom.name):
      return False
    if self.atom_id is not None and (self.atom_id != mdtraj_atom.index):
      return False
    if self.res_name is not None and (self.res_name != mdtraj_atom.residue.name):
      return False

    return True


class ContactFeature(object):
  def __init__(self, residue_i, residue_j):
    self.residue_i = residue_i
    self.residue_j = residue_j
