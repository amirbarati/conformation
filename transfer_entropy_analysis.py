from mdentropy.core.information import ncmutinf

"""
Compute Transfer Entropy between all pairs of features
from a list of timeseries data.
----------
Parameters:

featurized_timeseries: list
  List of Numpy arrays, each with shape (n_frames, n_features).
lag_time: int
  Lag time to be used in computation of Transfer Entropy
feature_names: list of str
  If your features have names (e.g., "Arg325" or "Temp.in Kansas"),
  you can optionally include them here.

Returns:
tentropy_pairs: Numpy array of shape (1, n_features^2-n_features), containing
  the transfer entropy between each possible pair of features.
tentropy_pairs_id_tuples: List of tuples of ints. 
 In same order as columns of tentropy_pairs.
 Each tuple contains two ints describing the ids of the features for which
 tentropy was calculated. 
tentropy_pairs_names: List of tuples of strings.

Caveats:
Implemented weighted mean is not optimal and may not be numerically stable.
This is an area for improvement.
"""
def compute_tentropy_pairs(featurized_timeseries, lag_time,
                           feature_names=None):

  total_frames = 0.
  n_vars = featurized_timeseries[0].shape[1]
  n_tent_pairs = n_vars ** 2 - n_vars
  tentropy_pairs = np.zeros((1, n_tent_pairs))
  tentropy_pairs_names = []

  for t_id, t in enumerate(featurized_timeseries):
    k = 0
    for i in range(0, t.shape[1]):
      for j in range(0, t.shape[1]):
        if i == j: 
          continue
        x = t[lag_time::lag_time, i]
        y = t[::lag_time, j][:-1]
        z = t[::lag_time, i][:-1]
        n_frames = x.shape[0]
        n_bins = np.floor(np.sqrt(n_frames / 5.))
        tent = ncmutinf(n_bins, x, y, z, n_bins)

        tentropy_pairs[k] += tent * n_frames
        if k==0:
          total_frames += n_frames

        k += 1

        if t_id == 0 and feature_names is not None:
          tentropy_pairs_names.append((feature_names[i], 
                                       feature_names[j]))

  tentropy_pairs /= total_frames

  return tentropy_pairs, tentropy_pairs_id_tuples, tentropy_pairs_names



