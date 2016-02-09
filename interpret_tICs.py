from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
import numpy as np 
import matplotlib.pyplot as plt
import pickle
import gzip
import mdtraj as md
from copy import deepcopy
from io_functions import *
import pandas as pd
import os

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
  n_tica_components=25, start_tIC=0):

  features = np.concatenate(load_file_list(None, directory = features_dir, ext = ".dataset"))
  tics = np.concatenate(load_file(projected_tica_coords))

  for j in range(start_tIC, np.shape(tics)[1]):
    #os.system("rm %s/tIC.%d_rfr.pkl" %(rf_dir, j))
    rfr = RandomForestRegressor(n_estimators=n_trees, max_features="sqrt",
                                n_jobs=8, verbose=1)
    print(("Fitting tree for tIC %d" %(j+1)))
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

def plot_importance_df(df, column_name, ylabel, title, analysis_type, save_string, j, save_dir):
  bar_width = [0.4 for i in range(0,50)]
  opacity = 0.8
  index = np.arange(50)
  
  df = df.reindex(df.importance.abs().order(ascending=False).index)
  df = df[:50]
  print(df)
  plt.barh(bottom=np.arange(50), width=df['importance'].values, height=bar_width, 
           alpha=opacity, color='b',label='Feature importance')
  plt.ylabel(ylabel)
  plt.xlabel('Overall Importance')
  plt.title(title)
  plt.yticks(index + bar_width, df[column_name].tolist(), fontsize=8)
  pp = PdfPages("%s/%s%d_%s.pdf" %(save_dir, analysis_type, (j+1), save_string))
  pp.savefig()
  pp.close()
  plt.clf()

def compute_residue_importances(df, percentile=95):
  residue_importances = {}
  for row in df.iterrows():
    row = row[1]
    res_i = row['res_i']
    if res_i not in list(residue_importances.keys()):
      residue_importances[res_i] = []
    residue_importances[res_i].append(np.abs(row['importance']))

    res_j = row['res_j']
    if res_j not in list(residue_importances.keys()):
      residue_importances[res_j] = []
    residue_importances[res_j].append(np.abs(row['importance']))

  for residue in residue_importances.keys():
    residue_importances[residue] = np.percentile(residue_importances[residue], percentile)

  return(residue_importances)

def merge_importances_features(feature_importances_file, feature_residues_map, importances=None):
  with open(feature_residues_map, "rb") as f:
    feature_tuples = pickle.load(f)
  feature_names = ["%s-%s" %(f[0].res_name.title(), f[1].res_name.title()) 
                                for f in feature_tuples]

  if importances is None:
    try:
      with open(feature_importances_file) as f:
        importances = pickle.load(f)
    except:
      with gzip.open(feature_importances_file) as f:
        importances = pickle.load(f)

  feature_importances = list(zip(feature_names, importances.tolist()))
  rows = []
  for i, feature_importance in enumerate(feature_importances):
    row = [feature_importance[0], 
                 feature_tuples[i][0].resSeq, feature_tuples[i][1].resSeq, 
                 feature_tuples[i][0].chain_id, feature_tuples[i][1].chain_id, 
                 feature_importance[1]]
    rows.append(row)
    
  df = pd.DataFrame(rows, columns=('feature_name', 'res_i', 'res_j', 'chain_i', 'chain_j', 'importance'))

  return df

def compute_per_residue_importance(df, percentile):
  net_df = pd.DataFrame(columns=['residue', 'resid', 'importance'])
  for _, row in df.iterrows():
    res_i = row['feature_name'].split("-")[0]
    res_j = row['feature_name'].split("-")[1]
    resid_i = row['res_i']
    resid_j = row['res_j']

    importance = row["importance"]
    if res_i not in net_df.index:
      net_df.loc[res_i] = [res_i, resid_i, [importance]]
    else:
      net_df.loc[res_i]["importance"].append(importance)
    if res_j not in net_df.index:
      net_df.loc[res_j] = [res_j, resid_j, [importance]]
    else:
      net_df.loc[res_j]["importance"].append(importance)

  net_df["importance"] = net_df["importance"].apply(np.percentile, q=percentile)

  return net_df  

def interpret_tIC_components(tica_filename, save_dir, feature_residues_pkl, n_tica_components=25, percentile=95):
  tica_object = verboseload(tica_filename)
  tica_components = tica_object.components_

  for j in range(0, np.shape(tica_components)[0]):
    components = tica_components[j,:]
    print(("Interpreting tIC %d" %(j+1)))
    feature_importances_df = merge_importances_features(None, feature_residues_pkl, components)
    residue_importances_df = compute_per_residue_importance(feature_importances_df, percentile)

    print("feature_importances_df.shape")
    print(feature_importances_df.shape)
    print("residue_importances_df.shape")
    print(residue_importances_df.shape)

    plot_importance_df(feature_importances_df, "feature_name", "Contact", 
                       "tIC %d Contact Component" % (j+1),
                       "tIC", 
                       "per_contact_importances", j, save_dir)

    plot_importance_df(residue_importances_df, "residue", "Residue", 
                       "tIC %d Residue Importances" % (j+1), 
                       "tIC",
                       "per_residue_importances", j, save_dir)

    with open("%s/tIC.%d_feature_importance_df.pkl" %(save_dir, j), "wb") as f:
      pickle.dump(feature_importances_df, f)

    with open("%s/tIC.%d_residue_importance_df.pkl" %(save_dir, j), "wb") as f:
      pickle.dump(residue_importances_df, f)

def interpret_tIC_rf(rf_dir, feature_residues_pkl, n_tica_components=25, percentile=95):
  for j in range(0, n_tica_components):
    print(("Interpreting tIC %d" %(j+1)))
    importances_file = "%s/tIC.%d_rfr.pkl" %(rf_dir, j)
    if os.path.exists(importances_file):
      feature_importances_df = merge_importances_features(importances_file, feature_residues_pkl)
      residue_importances_df = compute_per_residue_importance(feature_importances_df, percentile)

      plot_importance_df(feature_importances_df, "feature_name", "Contact", 
                         "tIC %d Contact Importances" % (j+1),
                         "tIC", 
                         "per_contact_importances", j, rf_dir)

      plot_importance_df(residue_importances_df, "residue", "Residue", 
                         "tIC %d Residue Importances" % (j+1), 
                         "tIC",
                         "per_residue_importances", j, rf_dir)

      with open("%s/tIC.%d_rfr_feature_importance_df.pkl" %(rf_dir, j), "wb") as f:
        pickle.dump(feature_importances_df, f)

      with open("%s/tIC.%d_rfr_residue_importance_df.pkl" %(rf_dir, j), "wb") as f:
        pickle.dump(residue_importances_df, f)


'''
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
'''

def plot_top_features_per_tIC(projected_tica_coords_file, features_dir, features_ext, 
                              rf_dir, n_tica_components, normalize=False, n_features=5):
  features = np.concatenate(load_file_list(get_trajectory_files(features_dir, features_ext)))
  projected_tica_coords = np.concatenate(load_file(projected_tica_coords_file))

  for j in range(0, n_tica_components):
    print(("Examining tIC %d" %(j+1)))
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




    


