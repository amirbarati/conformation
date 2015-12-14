from sklearn import mixture
from sklearn.ensemble import RandomForestClassifier
import numpy as np 
import matplotlib.pyplot as plt
from msmbuilder.dataset import dataset, _keynat, NumpyDirDataset
from msmbuilder.utils import verbosedump, verboseload
import time
import random
import os
import multiprocessing as mp
import csv
from io_functions import *
from functools import partial
import gzip, pickle
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import operator

from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com
from rpy2.robjects import r
import rpy2.robjects.numpy2ri as numpy2ri
numpy2ri.activate()
base = get_base()
R_functions = "%s/conformation/analysis.R" %base
R_analysis = "%s/conformation/b2ar_analysis.R" %base
ro.r.source(R_functions)
ro.r.source(R_analysis)


def select_model(X, tic_j, max_components, save_dir):
  num_trials = 0
  possible_components = []
  for trial in range(0,num_trials):
    train_indices = random.sample(range(0,len(X)), len(X)/2)
    test_indices = list(set(range(0,len(X))) - set(train_indices))
    X_train = X[train_indices,:]
    X_test = X[test_indices,:]
    
    bics = []
    aics = []
    test_likelihoods = []
    models = []
    for n_components in range(max_components, max_components):
        print("For tIC %d looking at GMM model with %d components" %(tic_j, n_components))
        g = mixture.DPGMM(n_components=10)
        #g = mixture.GMM(n_components = n_components, n_init = 3, min_covar = 1e-1, params='mc')
        g.fit(X_train)
        bics.append(g.bic(X_train))
        aics.append(g.aic(X_train))
        test_likelihoods.append(sum(g.score(X_test)))
        models.append(g)
    #plt.plot(range(1,max_components),aics)
    #minbic = bics.index(min(bics)) + 1
    #minaic = aics.index(min(aics)) + 1
    max_likelihood = test_likelihoods.index(max(test_likelihoods)) + 1
    #if(minbic) > 1: min_likelihood = test_likelihoods[1:].index(max(test_likelihoods[1:])) + 2
    #selected_components = min(minbic,minaic,min_likelihood)
    possible_components.append(max_likelihood)

  #num_components = min(possible_components)
  #g = mixture.GMM(n_components = num_components, n_init=5, tol =1e-5, min_covar = 1e-1, params='mc') 	 	
  g = mixture.DPGMM(n_components=5, alpha=10.0)
  g.fit(X)



  '''
  pickle.dump(bics, open("%s/tIC%d_gmm_bics.pkl" %(save_dir, tic_j), "wb"))
  plt.scatter(range(1,max_components),bics)
  pp = PdfPages("%s/tIC%d_gmm_bics.pdf" %(save_dir, tic_j))
  pp.savefig()
  pp.close()
  plt.clf()

  pickle.dump(aics, open("%s/tIC%d_gmm_aics.pkl" %(save_dir, tic_j), "wb"))
  plt.scatter(range(1,max_components),aics)
  pp = PdfPages("%s/tIC%d_gmm_aics.pdf" %(save_dir, tic_j))
  pp.savefig()
  pp.close()
  plt.clf()

  pickle.dump(test_likelihoods, open("%s/tIC%d_gmm_test_likelihoods.pkl" %(save_dir, tic_j), "wb"))
  plt.scatter(range(1,max_components),test_likelihoods)
  pp = PdfPages("%s/tIC%d_gmm_test_likelihoods.pdf" %(save_dir, tic_j))
  pp.savefig()
  pp.close()
  plt.clf()
  '''

  return(g)

def compute_gmm(tic_j_x_tuple, max_components, save_dir):
  print("Analyzing tIC %d" %(tic_j_x_tuple[0]))
  model = select_model(tic_j_x_tuple[1].reshape(-1,1), tic_j_x_tuple[0], max_components, save_dir)
  with gzip.open("%s/tIC%d_gmm.pkl.gz" %(save_dir, tic_j_x_tuple[0]), "wb") as f:
    pickle.dump(model, f)
  return

def compute_gmms(projected_tica_coords, save_dir, max_components):
  tics = np.concatenate(load_file(projected_tica_coords))
  compute_gmm_partial = partial(compute_gmm, max_components = max_components, save_dir = save_dir)
  #for pair in [(j, tics[:,j]) for j in range(0,np.shape(tics)[1])]:
   # compute_gmm_partial(pair)
  pool = mp.Pool(mp.cpu_count())
  pool.map(compute_gmm_partial, [(j, tics[:,j]) for j in range(0,np.shape(tics)[1])])
  pool.terminate()
  return

def compute_gmm_R(tIC_j_x_tuple, max_components, save_dir):
  j = tIC_j_x_tuple[0]
  tIC = tIC_j_x_tuple[1]
  gmm = r['compute.tIC.mixture.model'](tIC, j, save_dir, max_components=max_components, num_repeats=5)
  gmm_dict = { key : np.array(gmm.rx2(key)) for key in gmm.names}
  return(gmm_dict)

def compute_gmms_R(projected_tica_coords, max_components, save_dir, max_j=10):
  tics = np.concatenate(load_file(projected_tica_coords))
  #tics = tics[range(0,np.shape(tics)[0],10),:]
  compute_gmm_partial = partial(compute_gmm_R, max_components = max_components, save_dir = save_dir)
  #for pair in [(j, tics[:,j]) for j in range(0,np.shape(tics)[1])]:
   # compute_gmm_partial(pair)
  pool = mp.Pool(mp.cpu_count())
  gmms = pool.map(compute_gmm_partial, [(j, tics[:,j]) for j in range(0,max_j)])
  pool.terminate()

  with gzip.open("%s/all_gmms.pkl.gz" %save_dir, "wb") as f:
    pickle.dump(gmms, f)

  for j, gmm in enumerate(gmms):
    gmm_means = gmm['mu']
    with open("%s/tIC.%d_means.pkl" %(save_dir, j), "wb") as f:
      pickle.dump(gmm_means, f)

    gmm_classes = gmm['classes']
    with gzip.open("%s/tIC.%d_classes.pkl.gz" %(save_dir, j), "wb") as f:
      pickle.dump(gmm_classes, f)

  return(gmms)

def compute_rf_model(features, gmm):
  return

def plot_importances(feature_importances, save_dir, i, percentile=True):
  bar_width = 0.4
  opacity = 0.8
  index = np.arange(50)
  
  feature_importance_dict = {}

  if percentile:
    for feature_importance in feature_importances:
      residue_i = feature_importance[0].split("_")[0]
      residue_j = feature_importance[0].split("_")[1]
      importance = feature_importance[1]
      if residue_i not in feature_importance_dict.keys():
        feature_importance_dict[residue_i] = []
      if residue_j not in feature_importance_dict.keys():
        feature_importance_dict[residue_j] = []

      feature_importance_dict[residue_i].append(importance)
      feature_importance_dict[residue_j].append(importance)
    for residue in feature_importance_dict.keys():
      feature_importance_dict[residue] = np.percentile(feature_importance_dict[residue], 95.0)

    feature_importances = [(key, value) for key, value in feature_importance_dict.iteritems()]
    feature_importances.sort(key=operator.itemgetter(1),reverse=True)

  
  plt.barh(np.arange(50), [f[1] for f in feature_importances[0:50]], bar_width, alpha=opacity, color='b',label='Feature importance')
  plt.ylabel('Feature')
  plt.xlabel('Overall Importance')
  plt.title("Random-Forest + GMM Model of tIC %d Feature Importance" %(i+1))
  plt.yticks(index + bar_width, [f[0] for f in feature_importances[0:50]],fontsize=8)
  pp = PdfPages("%s/tIC%d_overall_importances.pdf" %(save_dir, i))
  pp.savefig()
  pp.close()
  plt.clf()

def plot_overall_rf_importances(rf_dir, feature_residues_map):
  filename = "%s/tIC%d_overall_rf.pkl" %(rf_dir, 0)
  j = 0
  feature_names = generate_features(feature_residues_map)
  feature_names = ["%d_%d" %(f[0], f[1]) for f in feature_names]
  
  while(os.path.exists(filename)):
    print("Examining tIC %d" %(j+1))
    with gzip.open(filename, "rb") as f:
      rf = pickle.load(f)

    feature_importances = zip(feature_names, rf.feature_importances_.astype(float).tolist())
    feature_importances.sort(key=operator.itemgetter(1),reverse=True)
    print(feature_importances[0:10])
    plot_importances(feature_importances, rf_dir, j)
    j += 1
    filename = "%s/tIC%d_overall_rf.pkl" %(rf_dir, j)


def plot_component_importances(df, save_dir, tic_i, component):
  df = df.sort_values(by="feature_importance")
  top_n = 50
  top_df = df.iloc[0:top_n]
  for i in range(0,top_n):
    if top_df.iloc[i]["component_mean"] < top_df.iloc[i]["non_component_mean"]:
        top_df.iloc[i][feature_importance] = top_df.iloc[i][feature_importance] * -1.0


  bar_width = 0.2
  plt.barh(np.arange(50), top_df["feature_importance"].tolist(), bar_width, alpha=opacity, color='b',label='Feature importance')
  bar_width = 0.2
  plt.xlabel('Feature')
  plt.ylabel('Overall Importance')
  plt.title("Random-Forest + GMM Model of tIC %d Feature Importance" %(tic_i+1))
  plt.xticks(index + bar_width, top_df["feature_name"].tolist(), rotation='vertical')
  pp = PdfPages("%s/tIC%d_c%d_vs_all_importances.pdf" %(save_dir, (tic_i+1), component))
  pp.savefig()
  pp.close()
  plt.clf()



def plot_column_pair(i, num_columns, save_dir, titles, data, gmm_means, refcoords):
  for j in range(i+1, num_columns):
    plt.hexbin(data[:,i],  data[:,j], bins = 'log', mincnt=1)
    print(gmm_means)
    print(gmm_means[i])
    print(gmm_means[j])
    for mean in gmm_means[i]:
      plt.axvline(x=mean,color='k',ls='dashed')
    for mean in gmm_means[j]:
      plt.axhline(y=mean,color='k',ls='dashed')
    if refcoords is not None:
      plt.scatter([refcoords[0,i]], [refcoords[0,j]], marker = 's', c='w',s=15)
      plt.scatter([refcoords[1,i]], [refcoords[1,j]], marker = 'v', c='k',s=15)
    if titles is not None: 
      pp = PdfPages("%s/%s_%s.pdf" %(save_dir, titles[i], titles[j]))
      plt.xlabel(titles[i])
      plt.ylabel(titles[j])
      pp.savefig()
      pp.close()
      plt.clf()
    else:
      pp = PdfPages("%s/tIC.%d_tIC.%d.pdf" %(save_dir, i+1, j+1))
      plt.xlabel("tIC.%d" %(i+1))
      plt.ylabel("tIC.%d" %(j+1))
      pp.savefig()
      pp.close()
      plt.clf()

def plot_tics_gmm(save_dir, data_file, gmm_dir, R=True, titles = None, tICA = False, scale = 1.0, refcoords_file = None):
  data = np.concatenate(load_file(data_file))

  if(refcoords_file is not None and os.path.exists(refcoords_file)):
    refcoords = load_file(refcoords_file)
  else:
    refcoords = None
  print(np.shape(refcoords))
  print(refcoords)

  gmm_means = []
  if not R:
    for j in range(0,np.shape(data)[1]):
      with gzip.open("%s/tIC%d_gmm.pkl.gz" %(gmm_dir, j)) as f:
        gmm = pickle.load(f)
      gmm_means.append(gmm.means_)
  else:
    for j in range(0,np.shape(data)[1]):
      with open("%s/tIC.%d_means.pkl" %(gmm_dir, j)) as f:
        means = pickle.load(f)
      gmm_means.append(means)    

  num_columns = np.shape(data)[1]
  plot_column_pair_partial = partial(plot_column_pair, num_columns = num_columns, save_dir = save_dir, titles = titles, 
    data = data, gmm_means = gmm_means, refcoords = refcoords)
  #for i in range(0,num_columns):
  #  plot_column_pair_partial(i)
  pool = mp.Pool(mp.cpu_count())
  pool.map(plot_column_pair_partial, range(0,num_columns))
  pool.terminate()

  print("Done plotting columns")
  return

def plot_tics_gmm_R(save_dir, data_file, gmm_dir, titles = None, tICA = False, scale = 1.0, refcoords_file = None):
  data = verboseload(data_file)
  data = np.concatenate(data)
  data[:,0] *= scale

  if(refcoords_file is not None):
    refcoords = load_file(refcoords_file)
  else:
    refcoords = None
  print(np.shape(refcoords))
  print(refcoords)

  gmm_means = []
  for j in range(0,np.shape(data)[1]):
    with gzip.open("%s/tIC%d_gmm.pkl.gz" %(gmm_dir, j)) as f:
      gmm = pickle.load(f)
    gmm_means.append(gmm.means_)

  num_columns = np.shape(data)[1]
  plot_column_pair_partial = partial(plot_column_pair, num_columns = num_columns, save_dir = save_dir, titles = titles, 
    data = data, gmm_means = gmm_means, refcoords = refcoords)
  #for i in range(0,num_columns):
  #  plot_column_pair_partial(i)
  pool = mp.Pool(mp.cpu_count())
  pool.map(plot_column_pair_partial, range(0,num_columns))
  pool.terminate()

  print("Done plotting columns")
  return


def compute_one_vs_all_rf_models(features_dir, projected_tica_coords, gmm_dir, save_dir, R=True, n_trees = 10, n_tica_components=25, feature_residues_map=""):
  features = np.concatenate(load_file_list(None, directory = features_dir, ext = ".dataset"))
  tics = np.concatenate(load_file(projected_tica_coords))
  feature_names = generate_features(feature_residues_map)
  feature_names = ["%d_%d" %(f[0], f[1]) for f in feature_names]

  for i in range(0, n_tica_components):
    print("Computing random forest model for tIC.%d" %(i+1))
    
    if R:
      with gzip.open("%s/tIC.%d_classes.pkl.gz" %(gmm_dir, i), "rb") as f:
        Y = pickle.load(f)
    else:
      with gzip.open("%s/tIC%d_gmm.pkl.gz" %(gmm_dir, i), "rb") as f:
        gmm = pickle.load(f)
      Y = gmm.predict(tics[:,i].reshape(-1,1))

    n_components = len(np.unique(Y).tolist())
    print(np.unique(Y).tolist())
    print("n_components %d" %n_components)

    for component in np.unique(Y).tolist():
      print("Analyzing component %d" %component)
      all_indices = range(0,np.shape(tics)[0])
      component_indices = [k for k in all_indices if int(Y[k]) == component]
      non_component_indices = list(set(all_indices)-set(component_indices))
      non_component_indices.sort()
      print(len(component_indices))
      print(len(non_component_indices))
      print(len(component_indices)+len(non_component_indices))
      print("Found indices")
      Z = copy.deepcopy(Y)
      Z[component_indices] = 0
      Z[non_component_indices] = 1
      rf = RandomForestClassifier(max_features="sqrt",bootstrap=True,n_estimators=n_trees,n_jobs=-1)
      print("fitting random forest model")
      r = rf.fit(features, Z)
      print("fit random forest model, saving now.")
      with gzip.open("%s/tIC%d_c%d_vs_all_rf.pkl.gz" %(save_dir, i, component), "wb") as f:
        pickle.dump(rf.feature_importances_.astype(float), f)
      print("Finished saving")

      feature_component_means = np.mean(features[component_indices,:], axis=0)
      print(np.shape(feature_component_means))
      print(feature_component_means[0:100])
      feature_non_component_means = np.mean(features[non_component_indices,:], axis=0)

      feature_importances = zip(feature_names, rf.feature_importances_.astype(float).tolist())
      pickle.dump(feature_importances, open("%s/tIC%d_c%d_importances_list.pkl" %(save_dir, i, component), "wb"))
      '''
      df = pd.DataFrame(columns=('feature_name', 'feature_importance', 'component_mean', 'non_component_mean'))

      for i, feature_importance in enumerate(feature_importances):
        df.iloc[k] = [feature_importance[0], feature_importance[1], feature_component_means[i], feature_non_component_means[i]]

      pickle.dump(df, open("%s/tIC%d_c%d_vs_all_df.pkl" %(save_dir, (i+1), component), "wb"))
      try:
        plot_component_importances(df, save_dir, i, component)
      except:
        continue
    '''


  return

def compute_overall_rf_models(features_dir, projected_tica_coords, gmm_dir, save_dir, R=True, n_trees = 10, n_tica_components=25, feature_residues_map = ""):
  '''
  if not os.path.exists("%s/feature_subset.h5" %save_dir):
    features = np.concatenate(load_file_list(None, directory = features_dir, ext = ".dataset"))
    features = features[range(0,np.shape(features)[0],100),:]
    verbosedump(features, "%s/feature_subset.h5" %save_dir)
  else:
    features = verboseload("%s/feature_subset.h5" %save_dir)
  '''

  features = np.concatenate(load_file_list(None, directory = features_dir, ext = ".dataset"))

  tics = np.concatenate(load_file(projected_tica_coords))
  feature_names = generate_features(feature_residues_map)
  feature_names = ["%d_%d" %(f[0], f[1]) for f in feature_names]

  for i in range(0, n_tica_components):
    print("Computing random forest model for tIC.%d" %(i+1))

    if R:
      with gzip.open("%s/tIC.%d_classes.pkl.gz" %(gmm_dir, i), "rb") as f:
        Y = pickle.load(f)
    else:
      with gzip.open("%s/tIC%d_gmm.pkl.gz" %(gmm_dir, i), "rb") as f:
        gmm = pickle.load(f)
      Y = gmm.predict(tics[:,i].reshape(-1,1))
      if len(gmm.means_) == 1: continue

    rf = RandomForestClassifier(max_features="sqrt",bootstrap=True,n_estimators=n_trees,n_jobs=-1)
    r = rf.fit(features, Y)

    with gzip.open("%s/tIC%d_overall_rf.pkl" %(save_dir, i), "wb") as f:
      pickle.dump(rf, f)
    print("Saved RF model")

    feature_importances = zip(feature_names, rf.feature_importances_.astype(float).tolist())
    feature_importances.sort(key=operator.itemgetter(1),reverse=True)

    print(feature_importances[0:20])
    try:
      plot_importances(feature_importances, save_dir, i)
    except:
      continue





