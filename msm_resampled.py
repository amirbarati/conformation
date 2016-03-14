import numpy as np
from msmbuilder.utils import verbosedump, verboseload


def resample_by_msm(total_samples, msm_object, clusters_map, num_trajs, save_file):
  equilibrium_populations = msm_object.populations_
  num_to_sample_per_cluster = {}
  for cluster_id in msm_object.mapping_.keys():
    state_id = msm_object.mapping_[cluster_id]
    num_to_sample_per_cluster[cluster_id] = np.rint(equilibrium_populations[state_id] * total_samples)

  print("Found number to sample per cluster based on equilibrium proporrtions.")
  sample_pairs = []
  for cluster_id in msm_object.mapping_.keys():
    traj_index_pairs = list(clusters_map[cluster_id])
    if len(traj_index_pairs) == 0:
      continue
    num_to_sample = num_to_sample_per_cluster[cluster_id]
    random_indices = np.random.choice(range(0, len(traj_index_pairs)), size=num_to_sample, replace=True)
    clusters_sample_pairs = [traj_index_pairs[i] for i in random_indices]
    sample_pairs += clusters_sample_pairs

  print("Obtained random (trajectory, frame) pairs based on equilibrium populations")
  #if there exists some fancy numpy way to index a 3d array by 2d tuples, then great, else:
  traj_to_frames = {}
  for i in range(0, num_trajs):
    traj_to_frames[i] = []

  for sample_pair in sample_pairs:
    traj_to_frames[sample_pair[0]].append(sample_pair[1])

  print("Rearranged equilibrium sampled frames based on trajectories")

  verbosedump(traj_to_frames, save_file)
  return traj_to_frames

def resample_features_by_msm_equilibirum_pop(features, traj_to_frames, save_file):
  resampled_features = []

  for traj_index, frames in traj_to_frames.iteritems():
    resampled_features.append(features[traj_index][frames, :])

  resampled_features = np.concatenate(resampled_features)

  verbosedump(resampled_features, save_file)

def generate_msm_traj_index_series(msm_object, start_cluster, n_steps, clusters_map, save_file):
  msm_trajectory = msm.sample_discrete(state=start_cluster, n_steps=n_steps)

  traj_index_pairs = []
  for cluster in msm_trajectory:
    traj_index_pair = np.random.choice(clusters_map[cluster], size=1)[0]
    traj_index_pairs.append(traj_index_pair)

  verbosedump(traj_index_pairs, save_file)

def generate_msm_traj_index_series(msm_object, start_cluster, end_cluster, clusters_map, save_file):
  msm_trajectory = msm.sample_discrete(state=start_cluster, n_steps=n_steps)

  traj_index_pairs = []
  for cluster in msm_trajectory:
    traj_index_pair = np.random.choice(clusters_map[cluster], size=1)[0]
    traj_index_pairs.append(traj_index_pair)

  verbosedump(traj_index_pairs, save_file)

def resample_features_by_msm_trajectory(features, traj_index_pairs, save_file):
  return
  