
# coding: utf-8

# In[7]:

# config to get better plot
get_ipython().magic('matplotlib inline')
get_ipython().magic("config InlineBackend.figure_format = 'retina'")

import matplotlib
matplotlib.rcParams['savefig.dpi'] = 1.5 * matplotlib.rcParams['savefig.dpi'] # larger image
matplotlib.rcParams['axes.labelcolor'] =  'green' # default green for label

matplotlib.rcParams['axes.linewidth'] =  0.5

from matplotlib import pyplot as plt
import seaborn as sb # add seaborn for pretty plot (vs default one in matplotlib)

import warnings # just to avoid any warning to make this notebook prettier
warnings.filterwarnings('ignore')


# In[3]:

# import pytraj
import pytraj as pt

# load sample data
traj = pt.datafiles.load_tz2()
traj


# In[4]:

# find hbond
hb = pt.hbond(traj)
print('donor - acceptor: ', hb._amber_mask())

print("")
print(hb.data)


# In[5]:

dist = pt.distance(traj, hb._amber_mask())
print('all hbond distances: ', dist)


# In[6]:

sb.color_palette('deep', n_colors=6, desat=0.5)
sb.set_style(style='white')

# scatter plot for distance between ':1@OG :2@H' and ':5@O :3@HG1'
# the point is colored by frame number (total frame = traj.n_frames (101))
fig = plt.scatter(dist[0], dist[1], marker='o', c=range(traj.n_frames), alpha=0.8, cmap='Spectral')
plt.colorbar()
plt.grid()
plt.xlabel(':1@OG :2@H')
plt.ylabel(':5@O :3@HG1')

