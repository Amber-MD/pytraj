
# coding: utf-8

# In[1]:

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


# In[2]:

# import pytraj
import pytraj as pt

# load sample data
traj = pt.iterload('tz2.nc', 'tz2.parm7')
traj


# In[3]:

# find hbond
hb = pt.hbond(traj)

distance_mask = hb.get_amber_mask()[0]
print('hbond distance mask: {} \n '.format(distance_mask))

angle_mask = hb.get_amber_mask()[1]
print('hbond angle mask: {} \n'.format(angle_mask))

print("hbond data")
print(hb.data) # 1: have hbond; 0: does not have hbond


# In[4]:

dist = pt.distance(traj, hb.get_amber_mask()[0])
print('all hbond distances: ', dist)


# In[5]:

angle = pt.angle(traj, hb.get_amber_mask()[1])
angle


# ### plot demo

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


# ### Do some statistics

# In[7]:

def stat(hb, data, distance_or_angle_mask):
    '''
    
    Parameters
    ----------
    hb : get from `pt.hbond(traj, ...)`
    data : either distance or angle
    distance_mask : distance mask
    '''
    import numpy as np
    arr = hb.data[1:].values == 1
    
    std_ = {}
    mean_ = {}
    for idx, mask in enumerate(distance_or_angle_mask):
        std_[mask] = np.std(data[idx][arr[idx]])
        mean_[mask] = np.mean(data[idx][arr[idx]])
    
    return mean_, std_


# In[8]:

stat(hb, dist, distance_mask)


# In[9]:

stat(hb, angle, angle_mask)

