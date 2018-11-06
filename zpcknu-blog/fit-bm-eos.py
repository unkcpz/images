#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# In[2]:


def birch_murnaghan(V, E0, V0, B0, B01):
    r = (V0 / V) ** (2. / 3.)
    return E0 + 9. / 16. * B0 * V0 * (r - 1.) ** 2 *                 (2. + (B01 - 4.) * (r - 1.))


def fit_birch_murnaghan_params(volumes_, energies_):
    from scipy.optimize import curve_fit

    volumes = np.array(volumes_)
    energies = np.array(energies_)
    params, covariance = curve_fit(
        birch_murnaghan, xdata=volumes, ydata=energies,
        p0=(
            energies.min(),  # E0
            volumes.mean(),  # V0
            0.1,  # B0
            3.,  # B01
        ),
        sigma=None
    )
    return params, covariance


# In[3]:


def plot_eos(volumes_, energies_):
    """
    Plots equation of state taking as input the pk of the ProcessCalculation 
    printed at the beginning of the execution of run_eos_wf 
    """
    import matplotlib.pyplot as plt
    
    volumes = np.array(volumes_)
    energies = np.array(energies_)

    params, covariance = fit_birch_murnaghan_params(volumes, energies)
    
    vmin = volumes.min()
    vmax = volumes.max()
    vrange = np.linspace(vmin, vmax, 300)

    plt.plot(volumes,energies,'o')
    plt.plot(vrange, birch_murnaghan(vrange, *params))

    plt.xlabel("Volume (ang^3)")
    # I take the last value in the list of units assuming units do not change
    plt.ylabel("Energy (eV)")
    plt.savefig("fit-bm.png")


# In[4]:


awithe = np.array(
    [[0.95, -10.54342866],
[0.96, -10.65937674],
[0.97, -10.74563203],
[0.98 , -10.80462694],
[0.99, -10.83863631],
[1.00 , -10.85001443],
[1.01 , -10.84077847],
[1.02, -10.81303532],
[1.03, -10.76877249],
[1.04, -10.70952633],
[1.05 ,-10.63683672]])
a = [i[0] for i in awithe]
e = [i[1] for i in awithe]


# In[5]:


volumes = [i**3 for i in a]


# In[6]:


p, c = fit_birch_murnaghan_params(volumes, e)


# In[7]:


print p, c


# In[9]:


plot_eos(volumes, e)


# In[ ]:




