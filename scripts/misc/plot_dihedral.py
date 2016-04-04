from pandas import DataFrame
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns

dir_cm = ''
dir_unmod = ''
fname = "chin_alpha_gamma.pk"
df_cm = DataFrame(pd.read_pickle(dir_cm + fname))
df_unmod = DataFrame(pd.read_pickle(dir_unmod + fname))

keys = df_cm.keys()

for key in keys:
    fig = plt.figure()
    ax = fig.add_subplot(120)
    ax.hist(df_cm[key], normed=1)
    ax.set_xlim([-180., 180.])
    ax.set_ylim([0, 0.1])
    ax = fig.add_subplot(121)
    ax.hist(df_unmod[key], normed=1)
    ax.set_ylim([0, 0.1])
    ax.set_xlim([-180., 180.])
    fname = key.replace(":", "_")
    plt.title(fname)
    plt.savefig("./plots/" + fname + ".png", dpi=300)
