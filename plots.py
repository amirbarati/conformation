import numpy as np
import matplotlib
matplotlib.style.use('ggplot')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from io_functions import *
import corner
import pandas as pd
import seaborn as sns


def plot_corner(data_file, plot_file, title, x_prefix, chosen_columns=None):
  data = np.concatenate(load_file(data_file))
  if chosen_columns is not None:
    data = data[:,chosen_columns]

  ptitles = []
  for j in range(0,np.shape(data)[1]):
    ptitles.append("%s %d" % (x_prefix, j))
  figure = corner.corner(data, labels=ptitles,
                           show_titles=True, title_args={"fontsize": 12}, color='rgb')
  figure.gca().annotate(title, xy=(0.5, 1.0), xycoords="figure fraction",
                        xytext=(0, -5), textcoords="offset points",
                        ha="center", va="top")
  pp = PdfPages(plot_file)
  pp.savefig(figure)
  pp.close()

def plot_seaborn(data_file, plot_file, title, x_prefix, chosen_columns=None):
  data = np.concatenate(load_file(data_file))
  if chosen_columns is not None:
    data = data[:,chosen_columns]

  ptitles = []
  for j in range(0,np.shape(data)[1]):
    ptitles.append("%s %d" % (x_prefix, j))

  df = pd.DataFrame(data, columns=ptitles)
  g = sns.PairGrid(df)
  g.map_diag(sns.kdeplot)
  g.map_offdiag(sns.kdeplot, cmap="Blues_d", n_levels=6);
  g.savefig(plot_file)

def plot_timescales(operator_file, plot_file, title):
  operator = verboseload(operator_file)
  timescales = operator.timescales_
  print("timescales")
  print(timescales)
  df = pd.DataFrame(columns=['timescales'])
  df['timescales'] = np.array(timescales)
  df.index = list(range(1,len(timescales)+1))
  plt.figure()
  df.plot(kind='bar'); plt.axhline(0, color='k')
  plt.title(title)
  pp = PdfPages(plot_file)
  pp.savefig()
  pp.close()  
