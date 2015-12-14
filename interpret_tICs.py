from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
import numpy as np 
import matplotlib.pyplot as plt
import pickle
import mdtraj as md
from copy import deepcopy
from io_functions import *
import pandas as pd

from matplotlib.backends.backend_pdf import PdfPages
from sklearn import preprocessing
import matplotlib
matplotlib.style.use('ggplot')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

'''Pseudocode:

compute_random_forests_fxn:
  Load and concat all features
  Load and concat tIC projections

  for each column in tIC array:
    do random forest regression from all feature space onto that tIC
    save resulting feature importances. 

interpret_tIC_rf:
  for each tIC:
    load feature importances 
    get feature names, combine and sort (use old code)
    plot a bar graph. 
    print out a matlines (or equivalent) plot for top 10 residues that matter.
    figure out a way to superimpose these on top of tIC heatmaps

'''

def compute_random_forests(features_dir, projected_tica_coords, rf_dir, R=True, n_trees=500, 
  n_tica_components=25):

  features = np.concatenate(load_file_list(None, directory = features_dir, ext = ".dataset"))
  tics = np.concatenate(load_file(projected_tica_coords))

  for j in range(0, np.shape(tics)[1]):
    rfr = RandomForestRegressor(n_estimators=n_trees, max_features="sqrt", oob_score=True,
                                n_jobs=16, verbose=1)
    print("Fitting tree for tIC %d" %(j+1))
    rfr.fit(features, tics[:,j])
    print("Done fitting RF model. Saving now...")
    with open("%s/tIC.%d_rfr.pkl" %(rf_dir, j), "wb") as f:
      pickle.dump(rfr.feature_importances_, f)
    print("Saved random forest feature importances")

  return

def get_feature_list(feature_residues_csv, structure_file):
  feature_names = generate_features(feature_residues_csv)
  structure = md.load_frame(structure_file, index=0)
  all_resSeq = [r.resSeq for r in structure.topology.residues if r.is_protein]
  all_res = [str(r).title() for r in structure.topology.residues if r.is_protein]

  feature_list = []
  for i, feature_name in enumerate(feature_names):
    try:
      res_i = int(feature_name[0])
      res_j = int(feature_name[1])
    except:
      res_i = int(feature_name[0][1])
      res_j = int(feature_name[1][1])

    res_i_idx = all_resSeq.index(res_i)
    res_i_name = all_res[res_i_idx]

    res_j_idx = all_resSeq.index(res_j)
    res_j_name = all_res[res_j_idx]

    feature_list.append((res_i_name, res_j_name))

  return(feature_list)

def plot_importance_df(df, save_string, j, save_dir):
  bar_width = 0.4
  opacity = 0.8
  index = np.arange(50)

  df.sort(columns='importance', inplace=True, ascending=0)
  df = df[:50]
  plt.barh(np.arange(50), pd.np.array(df["importance"]), bar_width, 
          alpha=opacity, color='b',label='Feature importance')
  plt.ylabel('Feature')
  plt.xlabel('Overall Importance')
  plt.title("Random-Forest + GMM Model of tIC %d Feature Importance" %(j+1))
  plt.yticks(index + bar_width, df["feature_name"].tolist(), fontsize=8)
  pp = PdfPages("%s/tIC%d_%s.pdf" %(save_dir, (j+1), save_string))
  pp.savefig()
  pp.close()
  plt.clf()

def compute_residue_importances(df, percentile=95):
  residue_importances = {}
  for row in df.iterrows():
    row = row[1]
    res_i = row['res_i']
    if res_i not in residue_importances.keys():
      residue_importances[res_i] = []
    residue_importances[res_i].append(row['importance'])

    res_j = row['res_j']
    if res_j not in residue_importances.keys():
      residue_importances[res_j] = []
    residue_importances[res_j].append(row['importance'])

  for residue in residue_importances.iterkeys():
    residue_importances[residue] = np.percentile(residue_importances[residue], percentile)

  return(residue_importances)


def interpret_tIC_rf(rf_dir, feature_residues_csv, structure_file, n_tica_components=25):
  df = pd.DataFrame(columns=('index', 'res_i', 'res_j', 'feature_name', 'importance'))
  feature_list = get_feature_list(feature_residues_csv, structure_file)
  for i, feature in enumerate(feature_list):
    df.loc[i] = [i, feature[0], feature[1], "", 0.0]

  for j in range(0,n_tica_components):
    print("Interpreting tIC %d" %(j+1))
    df_j = deepcopy(df)
    with open("%s/tIC.%d_rfr.pkl" %(rf_dir, j), "rb") as f:
      feature_importances = pickle.load(f)

    df_j['importance'] = feature_importances
    concatenated_feature_names = zip(df['res_i'].tolist(), df['res_j'].tolist())
    concatenated_feature_names = ["%s-%s" %(res_i, res_j) for res_i, res_j in concatenated_feature_names]
    df_j['feature_name'] = concatenated_feature_names
    with open("%s/tIC.%d_rfr_feature_importance_df.pkl" %(rf_dir, j), "wb") as f:
      pickle.dump(df_j, f)
    plot_importance_df(df_j, save_string="per_feature_rf_importance", j=j, save_dir=rf_dir)



    residue_df = pd.DataFrame(columns=('feature_name', 'importance'))
    residue_importances = compute_residue_importances(df_j)
    for i, (key, value) in enumerate(residue_importances.iteritems()):
      residue_df.loc[i] = [key, value]
    plot_importance_df(residue_df, save_string="per_residue_rf_importance", j=j, save_dir=rf_dir)

    with open("%s/tIC.%d_rfr_residue_importance_df.pkl" %(rf_dir, j), "wb") as f:
      pickle.dump(residue_df, f)

def plot_top_features_per_tIC(projected_tica_coords_file, features_dir, features_ext, 
                              rf_dir, n_tica_components, normalize=False, n_features=5):
  features = np.concatenate(load_file_list(get_trajectory_files(features_dir, features_ext)))
  projected_tica_coords = np.concatenate(load_file(projected_tica_coords_file))

  for j in range(0, n_tica_components):
    print("Examining tIC %d" %(j+1))
    with open("%s/tIC.%d_rfr_feature_importance_df.pkl" %(rf_dir, j), "rb") as f:
      df = pickle.load(f)
    df.sort(columns='importance', inplace=True, ascending=0)
    indices = df[:n_features]['index'].tolist()
    top_features_normalized = features[:,indices]
    if normalize:
      top_features_normalized = preprocessing.scale(top_features_normalized)
    top_features_normalized = top_features_normalized[projected_tica_coords[:,j].argsort(),:]
    projected_tica_coords_sorted = projected_tica_coords[projected_tica_coords[:,j].argsort(),:]
    
    df = pd.DataFrame(top_features_normalized, 
                      index=projected_tica_coords_sorted[:,j], 
                      columns=df[:n_features]['feature_name'])

    plt.figure()
    pd.rolling_mean(df, 1000).plot()
    plt.legend(loc='center left', bbox_to_anchor=(1.0,0.5))
    plt.title("Feature Behavior Over Course of tIC %d" %(j+1))
    plt.xlabel("tIC %d Projection" %(j+1))
    plt.ylabel("Distance (nm)")
    pp = PdfPages("%s/tIC%d_feature_behavior.pdf" %(rf_dir, (j+1)))
    pp.savefig(bbox_inches='tight')
    pp.close()
    plt.clf()



def top_residues_per_tIC():
  return




    


