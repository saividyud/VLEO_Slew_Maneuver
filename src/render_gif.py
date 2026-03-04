#!/usr/bin/env python3
# /// script
# dependencies = ["matplotlib", "scipy", "imageio"]
# ///
"""Render volcano observation animation to GIF from MATLAB .mat data.
3 views: along h, along volcano projection perp to h, and complementary axis."""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.io as sio
import imageio.v2 as imageio
import os, io

# Load data
script_dir = os.path.dirname(os.path.abspath(__file__))
data = sio.loadmat(os.path.join(script_dir, 'anim_data.mat'))

tspan = data['tspan'].flatten()
X_history = data['X_history']
visibility = data['visibility_history'].flatten()
pointing = data['pointing_eci_history']
r_volcano_0 = data['r_volcano_0'].flatten()
R_e = data['params']['R_e'][0, 0].item()
omega_earth = data['params']['omega_earth'][0, 0].item()
dt = tspan[1] - tspan[0]

# Find visibility window with padding
vis_indices = np.where(visibility > 0)[0]
if len(vis_indices) == 0:
    print('No visibility window found, using full time range.')
    i_start, i_end = 0, len(tspan) - 1
else:
    pad = int(60 / dt)
    i_start = max(0, vis_indices[0] - pad)
    i_end = min(len(tspan) - 1, vis_indices[-1] + pad)

# ~40 frames
window_len = i_end - i_start + 1
frame_skip = max(1, window_len // 40)
frame_indices = list(range(i_start, i_end + 1, frame_skip))

# Earth sphere
u = np.linspace(0, 2 * np.pi, 40)
v = np.linspace(0, np.pi, 25)
xs = R_e * np.outer(np.cos(u), np.sin(v))
ys = R_e * np.outer(np.sin(u), np.sin(v))
zs = R_e * np.outer(np.ones_like(u), np.cos(v))

# --- Compute 3 orthogonal view axes ---
mid_i = vis_indices[len(vis_indices) // 2] if len(vis_indices) > 0 else len(tspan) // 2
r_mid = X_history[mid_i, :3]
v_mid = X_history[mid_i, 3:6]

# e1: angular momentum direction (h = r x v)
h_vec = np.cross(r_mid, v_mid)
e1 = h_vec / np.linalg.norm(h_vec)

# Volcano at mid time
t_mid = tspan[mid_i]
cos_wt = np.cos(omega_earth * t_mid)
sin_wt = np.sin(omega_earth * t_mid)
r_volc_mid = np.array([
    cos_wt * r_volcano_0[0] - sin_wt * r_volcano_0[1],
    sin_wt * r_volcano_0[0] + cos_wt * r_volcano_0[1],
    r_volcano_0[2]
])

# e2: volcano direction projected perp to h
r_volc_unit = r_volc_mid / np.linalg.norm(r_volc_mid)
e2 = r_volc_unit - np.dot(r_volc_unit, e1) * e1
e2 = e2 / np.linalg.norm(e2)

# e3: complementary
e3 = np.cross(e1, e2)
e3 = e3 / np.linalg.norm(e3)


def vec_to_elev_azim(d):
    elev = np.degrees(np.arcsin(np.clip(d[2], -1, 1)))
    azim = np.degrees(np.arctan2(d[1], d[0]))
    return elev, azim


views = [
    (vec_to_elev_azim(e1), 'Along h'),
    (vec_to_elev_azim(e2), 'Along volcano proj.'),
    (vec_to_elev_azim(e3), 'Orthogonal'),
]

# Center on midpoint between satellite and volcano at mid-time, zoom to interaction
center = (r_mid + r_volc_mid) / 2
# Range: distance from center to satellite, with margin
dist_sat = np.linalg.norm(r_mid - center)
dist_volc = np.linalg.norm(r_volc_mid - center)
max_range = max(dist_sat, dist_volc, R_e * 0.15) * 2.5

gif_path = os.path.join(script_dir, 'volcano_observation.gif')
frames = []
n_total = len(frame_indices)

print(f'Rendering {n_total} frames x 3 views...')

for fi, i in enumerate(frame_indices):
    t = tspan[i]
    fig = plt.figure(figsize=(20, 7), dpi=120, facecolor='white')

    r_sat = X_history[i, :3]
    cos_wt = np.cos(omega_earth * t)
    sin_wt = np.sin(omega_earth * t)
    r_volc = np.array([
        cos_wt * r_volcano_0[0] - sin_wt * r_volcano_0[1],
        sin_wt * r_volcano_0[0] + cos_wt * r_volcano_0[1],
        r_volcano_0[2]
    ])

    vis_color = '#2ecc71' if visibility[i] > 0 else '#e74c3c'
    vis_text = 'VISIBLE' if visibility[i] > 0 else 'NOT VISIBLE'

    for vi, ((elev, azim), label) in enumerate(views):
        ax = fig.add_subplot(1, 3, vi + 1, projection='3d')

        # Earth (partial, clipped by view)
        ax.plot_surface(xs, ys, zs, alpha=0.25, color=[0.3, 0.5, 0.8],
                        edgecolor='none', rcount=25, ccount=25)

        # Satellite
        ax.scatter(*r_sat, c='red', s=80, zorder=5, depthshade=False)

        # Trajectory
        ax.plot(X_history[i_start:i+1, 0], X_history[i_start:i+1, 1],
                X_history[i_start:i+1, 2], 'r-', linewidth=1.5, alpha=0.7)

        # Volcano
        ax.scatter(*r_volc, c='#FFD700', edgecolors='black', s=120,
                   zorder=5, depthshade=False, linewidths=1.5)

        # Visibility line
        if visibility[i] > 0:
            ax.plot([r_volc[0], r_sat[0]], [r_volc[1], r_sat[1]],
                    [r_volc[2], r_sat[2]], color='#2ecc71', linewidth=2.5)
        else:
            ax.plot([r_volc[0], r_sat[0]], [r_volc[1], r_sat[1]],
                    [r_volc[2], r_sat[2]], 'r--', linewidth=1.5, alpha=0.5)

        # Camera pointing
        pv = pointing[i] * 2e6
        ax.quiver(r_sat[0], r_sat[1], r_sat[2], pv[0], pv[1], pv[2],
                  color='#3498db', linewidth=2, arrow_length_ratio=0.15)

        # Zoom to interaction area
        ax.set_xlim(center[0] - max_range, center[0] + max_range)
        ax.set_ylim(center[1] - max_range, center[1] + max_range)
        ax.set_zlim(center[2] - max_range, center[2] + max_range)
        ax.set_box_aspect([1, 1, 1])

        # Clean up - remove clutter
        ax.set_axis_off()
        ax.set_title(label, fontsize=13, fontweight='bold', pad=5)
        ax.view_init(elev=elev, azim=azim)

    fig.suptitle(f't = {t/60:.2f} min    {vis_text}',
                 fontsize=16, fontweight='bold', color=vis_color, y=0.97)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.01, wspace=0.02)

    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0.05)
    buf.seek(0)
    frames.append(imageio.imread(buf))
    buf.close()
    plt.close(fig)

    if (fi + 1) % 10 == 0:
        print(f'  {fi + 1}/{n_total} frames done')

print('Writing GIF...')
imageio.mimsave(gif_path, frames, duration=0.15, loop=0)
print(f'GIF saved to: {gif_path}')
