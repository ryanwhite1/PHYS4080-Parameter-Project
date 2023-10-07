# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 11:04:13 2023

@author: ryanw
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import animation
import camb
plt.style.use('dark_background')
# set LaTeX font for our figures
plt.rcParams.update({"text.usetex": True})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'

def run_camb(Omega_bh2=0.02242, Omega_cdmh2=0.11933, Omega_k=0.0, H0=67.66, DE_EoS=-1.0, 
             scalar_amplitude = 2.105e-9, max_pk_redshift=10.0):
    my_cosmology = camb.set_params(ombh2 = Omega_bh2, omch2 = Omega_cdmh2, omk = Omega_k, H0=H0, 
                                   w=DE_EoS, As=scalar_amplitude, WantCls=True, 
                                   WantTransfer=True, WantDerivedParameters=True, lmax=2500,
                                   redshifts=np.concatenate([np.logspace(np.log10(max_pk_redshift), -2.0, 100),[0.0]]))
    run = camb.get_results(my_cosmology)
    TT, EE, BB, TE = np.split(run.get_cmb_power_spectra(spectra=['total'], CMB_unit='muK')['total'],4,axis=1)
    Pk_interpolator = run.get_matter_power_interpolator()
    return [my_cosmology, run, TT, EE, BB, TE, Pk_interpolator]

Omega_b_vals = np.linspace(0.02, 0.07, 25)
Omega_cdm_vals = np.linspace(0.1, 0.3, 10)
results = []
for i in range(len(Omega_b_vals)):
    for j in range(len(Omega_cdm_vals)):
        results.append(run_camb(Omega_bh2=Omega_b_vals[i]*0.6766**2, Omega_cdmh2=Omega_cdm_vals[j]*0.6766**2))
    print(i)
    
TT_new = results[0][2]
lmax = 2500
l = np.linspace(2, lmax, len(TT_new) - 2)

colours = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=min(Omega_cdm_vals), vmax=max(Omega_cdm_vals)), 
                                       cmap='cool')

# now calculate some parameters for the animation frames and timing
length = 5
nt = len(Omega_b_vals) # number of timesteps
frames = np.arange(0, nt, 1)    # iterable for the animation function. Chooses which frames (indices) to animate.
fps = len(frames) // length  # fps for the final animation

### --- TT power spectra --- ###
fig, axes = plt.subplots(figsize=(12, 7), ncols=2, width_ratios=[20, 1], gridspec_kw={"wspace":0})

def animate(i):
    for ax in axes:
        ax.clear()
    index = i * len(Omega_cdm_vals)
    for j in range(len(Omega_cdm_vals)):
        axes[0].plot(l, results[index + j][2][2:], c=colours.to_rgba(Omega_cdm_vals[j]))
    axes[0].set(xlabel="Multipole Moment $l$", ylabel="Anisotropy Power ($\mu K^2$)",
           xscale='log',
           title=f'$\Omega_b = {round(Omega_b_vals[i], 3)}$',
           ylim=(0, 8500))
    fig.colorbar(mappable=colours, cax=axes[1], label="$\Omega_{cdm}$")
    return fig,

axes[0].set_facecolor('k')   # black background
ani = animation.FuncAnimation(fig, animate, frames=frames, interval=1, cache_frame_data=False, blit=True)
ani.save("TTPowerSpectrum.gif", writer='pillow', fps=fps, dpi=300)


### --- TE power spectra --- ###
fig, axes = plt.subplots(figsize=(12, 7), ncols=2, width_ratios=[20, 1], gridspec_kw={"wspace":0})

def animate(i):
    for ax in axes:
        ax.clear()
    index = i * len(Omega_cdm_vals)
    for j in range(len(Omega_cdm_vals)):
        axes[0].plot(l, results[index + j][5][2:], c=colours.to_rgba(Omega_cdm_vals[j]))
    axes[0].set(xlabel="Multipole Moment $l$", ylabel="Anisotropy Power ($\mu K^2$)",
           xscale='log',
           title=f'$\Omega_b = {round(Omega_b_vals[i], 3)}$',
           ylim=(-190, 190))
    fig.colorbar(mappable=colours, cax=axes[1], label="$\Omega_{cdm}$")
    return fig,

axes[0].set_facecolor('k')   # black background
ani = animation.FuncAnimation(fig, animate, frames=frames, interval=1, cache_frame_data=False, blit=True)
ani.save("TEPowerSpectrum.gif", writer='pillow', fps=fps, dpi=300)


### --- EE power spectra --- ###
fig, axes = plt.subplots(figsize=(12, 7), ncols=2, width_ratios=[20, 1], gridspec_kw={"wspace":0})

def animate(i):
    for ax in axes:
        ax.clear()
    index = i * len(Omega_cdm_vals)
    for j in range(len(Omega_cdm_vals)):
        axes[0].plot(l, results[index + j][3][2:], c=colours.to_rgba(Omega_cdm_vals[j]))
    axes[0].set(xlabel="Multipole Moment $l$", ylabel="Anisotropy Power ($\mu K^2$)",
           xscale='log',
           title=f'$\Omega_b = {round(Omega_b_vals[i], 3)}$',
           ylim=(0, 70))
    fig.colorbar(mappable=colours, cax=axes[1], label="$\Omega_{cdm}$")
    return fig,

axes[0].set_facecolor('k')   # black background
ani = animation.FuncAnimation(fig, animate, frames=frames, interval=1, cache_frame_data=False, blit=True)
ani.save("EEPowerSpectrum.gif", writer='pillow', fps=fps, dpi=300)


### --- BB power spectra --- ###
fig, axes = plt.subplots(figsize=(12, 7), ncols=2, width_ratios=[20, 1], gridspec_kw={"wspace":0})

def animate(i):
    for ax in axes:
        ax.clear()
    index = i * len(Omega_cdm_vals)
    for j in range(len(Omega_cdm_vals)):
        axes[0].plot(l, results[index + j][4][2:], c=colours.to_rgba(Omega_cdm_vals[j]))
    axes[0].set(xlabel="Multipole Moment $l$", ylabel="Anisotropy Power ($\mu K^2$)",
           xscale='log',
           title=f'$\Omega_b = {round(Omega_b_vals[i], 3)}$',
           ylim=(0, 0.14))
    fig.colorbar(mappable=colours, cax=axes[1], label="$\Omega_{cdm}$")
    return fig,

axes[0].set_facecolor('k')   # black background
ani = animation.FuncAnimation(fig, animate, frames=frames, interval=1, cache_frame_data=False, blit=True)
ani.save("BBPowerSpectrum.gif", writer='pillow', fps=fps, dpi=300)



### --- Matter Power Spectrum --- ###
plt.style.use('default')
plt.rcParams.update({"text.usetex": True})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'

kvalues = np.logspace(-4.0, 0.0, 1000)
# zs = np.append([0], np.geomspace(0.001, 30, 6))
zs = np.linspace(0, 10, 4)
z_colours = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=min(zs), vmax=max(zs)), 
                                       cmap='winter')

fig, axes = plt.subplots(figsize=(12, 7), ncols=2, width_ratios=[20, 1], gridspec_kw={"wspace":0})

i = 0
for z in zs:
    if i == 0:
        lab1 = f"$\Omega_{{cdm}} = {min(Omega_cdm_vals)}$"
        lab2 = f"$\Omega_{{cdm}} = {max(Omega_cdm_vals)}$"
        i = 1
    else:
        lab1 = ''; lab2 = ''
    axes[0].plot(kvalues, results[14 * len(Omega_cdm_vals)][6].P(z, kvalues), c=z_colours.to_rgba(z), ls='--', label=lab1)
    axes[0].plot(kvalues, results[15 * len(Omega_cdm_vals) - 1][6].P(z, kvalues), c=z_colours.to_rgba(z), ls=':', label=lab2)

axes[0].set(xscale='log', yscale='log', 
       xlabel="Wavenumber $k$ ($h$ Mpc$^{-1}$)", ylabel="$P(k)$ ($h^{-3}$ Mpc$^3$)")
axes[0].legend()
fig.colorbar(mappable=z_colours, cax=axes[1], label='$z$')
fig.savefig("MatterPowerSpectrum.png", dpi=400, bbox_inches='tight')




### --- Matter Power Spectrum gif --- ###
length = 5
nt = len(Omega_cdm_vals) # number of timesteps
frames = np.arange(0, nt, 1)    # iterable for the animation function. Chooses which frames (indices) to animate.
fps = len(frames) // length  # fps for the final animation

fig, axes = plt.subplots(figsize=(12, 7), ncols=2, width_ratios=[20, 1], gridspec_kw={"wspace":0})

def animate(i):
    for ax in axes:
        ax.clear()
    for z in zs:
        axes[0].plot(kvalues, results[14 * len(Omega_cdm_vals) + i][6].P(z, kvalues), c=z_colours.to_rgba(z))
    axes[0].set(xscale='log', yscale='log', 
           xlabel="Wavenumber $k$ ($h$ Mpc$^{-1}$)", ylabel="$P(k)$ ($h^{-3}$ Mpc$^3$)",
           ylim=(0.1, 5e4),
           title=f"$\Omega_{{cdm}} = {round(Omega_cdm_vals[i], 3)}$")
    fig.colorbar(mappable=z_colours, cax=axes[1], label='$z$')
    return fig,
        
ani = animation.FuncAnimation(fig, animate, frames=frames, interval=1, cache_frame_data=False, blit=True)
ani.save("MatterPowerSpectrum.gif", writer='pillow', fps=fps, dpi=300)

plt.close('all')