from residue import Residue

def convert_list_to_resobject_list(contact_residues):
  resobject_list = []
  for residue in contact_residues:
    new_residue = Residue(resSeq = residue[1], chain_id = residue[0])
    resobject_list.append(new_residue)
  return(resobject_list)

#CHOOSE RESIDUES:
helix_residues = {}
helix_residues["tm1"] = range(29,61)
helix_residues["icl1"] = range(61,67)
helix_residues["tm2"] = range(67,97)
helix_residues["ecl1"] = range(97,102)
helix_residues["tm3"] = range(102,137)
helix_residues["icl2"] = range(137,147)
helix_residues["tm4"] = range(147,171)
helix_residues["ecl2"] = range(171,197)
helix_residues["tm5"] = range(197,230)
helix_residues["icl3"] = range(230,267)
helix_residues["tm6"] = range(267,299)
helix_residues["ecl3"] = range(299,305)
helix_residues["tm7"] = range(305,330)

cutoff = 1.0


helix_residues["tm8"] = range(330,340)

switch_residues = [130, 131, 208, 211, 219, 268, 272, 286, 288, 316, 322, 323, 326]
switch_npxx = [130, 131, 208, 211, 219, 268, 272, 286, 288, 316] + range(322,328)
switch_pp_npxx = set(switch_npxx + [51, 79, 106, 113, 118, 121, 130, 131, 132, 141, 158, 208, 211, 219, 268, 272, 282, 285, 286, 288, 316, 318, 319, 320, 323, 326, 282])
tm6_residues = range(270,299)
bp_residues = [82, 86, 93, 106, 110, 113, 114, 117, 118, 164, 191, 192, 193, 195, 199, 200, 203, 206, 208, 286, 289, 290, 293, 305, 308, 309, 312, 316]
dihedral_residues = list(set(switch_npxx + tm6_residues))
skip_5_residues = range(30,340,5)
skip_3_residues = range(30,340,3)
skip5_switches_pp_npxx = list(set(skip_5_residues + list(switch_pp_npxx)))
skip5_switches_pp_npxx_ser = list(set(skip_5_residues + list(switch_pp_npxx) + [207]))
skip3_switches_pp_npxx = list(set(skip_3_residues + list(switch_pp_npxx)))
#print(len(skip5_switches_pp_npxx))
all_residues = range(29,340)
tm_residues = helix_residues["tm1"] + helix_residues["tm2"] + helix_residues["tm3"] + helix_residues["tm4"] + helix_residues["tm5"] + helix_residues["tm6"] + helix_residues["tm7"] + helix_residues["tm8"]
sampling_method = "random"
precision = "SP"

#feature_types = "_switches_tm6"
#feature_types = "_switches_npxx_tm6_bp"
#feature_types = "_switches_npxx_tm6_dihedrals_switches_npxx_contact"
#feature_types = "_switches_npxx_tm6_dihedrals_switches_pp_npxx_contact"\
#feature_types = "_skip5_switches_pp_npxx_contact"
#feature_types = "_skip5_switches_pp_npxx_contact_cutoff%dnm" %(int(cutoff))
#feature_types = "_skip3_switches_pp_npxx_contact_cutoff20"
#feature_types = "switches_pp_npxx_contact_cutoff10000"
#feature_types = "skip5_switches_pp_npxx_ser_cutoff%dnm" %(int(cutoff))
#feature_types = "all_residues_under_cutoff%dnm" %(int(cutoff))
#feature_types = "all_residues_under_cutoff%dnm_allframes" %(int(cutoff))
#feature_types = "all_tm_residues_under_cutoff%dnm" %(int(cutoff))
#feature_types = "reimaged_notrajfix_tm_residues_under_cutoff%dnm" %(int(cutoff))
#feature_types = "reimaged_notrajfix_tm_residues_2rh1_3sn6_under_cutoff%dnm" %(int(cutoff))
#feature_types = "reimaged_notrajfix_all_residues_under_cutoff%dnm" %(int(cutoff))
#feature_types = "all_residues_2rh1_3sn6_under_cutoff%dnm" %(int(cutoff))
feature_name = "tm_residues_2rh1_3sn6_under_cutoff%dnm" %(int(cutoff))
#feature_types = "reference_receptors"
contact_residues = [(0, res) for res in tm_residues]
#contact_residues = [(0, res) for res in all_residues]

contact_residues = convert_list_to_resobject_list(contact_residues)