"""
Droplet density analysis from LAMMPS dump data - Plot fitted ellipsoid and Convex Hull points
"""

# Standart python packages
import numpy as np

# Packages to install with "pip3 install -U package-name"
#
#   Visualisation
#
#      Plotly static image generation requires Kaleido: pip3 install -U kaleido
import plotly.graph_objects as go

def fitted_ellipsoid(output_html, output_png, frame_number, dens_nonzero_xyz_full, hull_vertices_xyz, z_ref, z_ref_ellipse_params, el_center, el_radii, x_intersection_point, y_intersection_point, z_intersection_point):

    # The parametric equations of an ellipsoid of center (x0, y0, z0) and semi-axes A, B, C are:
    #   x = x0 + A*cos(u)*sin(v)
    #   y = y0 + B*sin(u)*sin(v)
    #   z = z0 + C*cos(v)
    # for u in [0,2pi] and v in [0,pi]  ===> whole ellipsoid
    # for u in [0,2pi] and v in [0,pi/2]  ===> top half ellipsoid: z ≥ 0
    # for u in [0,pi] and v in [0,pi/2]  ===> half of the top half ellipsoid: quarter of the ellipsoid
    #
    npts=100
    u = np.linspace(0.0, 2.0 * np.pi, npts)
    v = np.linspace(0.0, np.pi, npts)
    #
    u, v = np.meshgrid(u, v)
    #
    x_ellipsiod = el_center[0] + el_radii[0] * np.cos(u) * np.sin(v)
    y_ellipsiod = el_center[1] + el_radii[1] * np.sin(u) * np.sin(v)
    z_ellipsiod = el_center[2] + el_radii[2] * np.cos(v)
    #
    x_ellipsiod_flatten = x_ellipsiod.flatten()
    y_ellipsiod_flatten = y_ellipsiod.flatten()
    z_ellipsiod_flatten = z_ellipsiod.flatten()
    #
    zindx = np.where(z_ellipsiod_flatten >= z_ref)
    x_ellipsiod_filtered = x_ellipsiod_flatten[zindx]
    y_ellipsiod_filtered = y_ellipsiod_flatten[zindx]
    z_ellipsiod_filtered = z_ellipsiod_flatten[zindx]
    #print("x_ellipsiod_filtered.shape: ", x_ellipsiod_filtered.shape)
    #print("y_ellipsiod_filtered.shape: ", y_ellipsiod_filtered.shape)
    #print("z_ellipsiod_filtered.shape: ", z_ellipsiod_filtered.shape)
    
    # Use Plotly to plot 2D ellipse and Convex Hull points
    #
    fig = go.Figure()
    #
    fig.add_scatter3d(
        x=dens_nonzero_xyz_full[:, 0],
        y=dens_nonzero_xyz_full[:, 1],
        z=dens_nonzero_xyz_full[:, 2],
        mode='markers',
        marker=dict(size=7, color='rgba(135, 206, 250)',symbol='cross'),
        opacity=0.8
    )
    #
    fig.add_scatter3d(
        x=hull_vertices_xyz[:, 0],
        y=hull_vertices_xyz[:, 1],
        z=hull_vertices_xyz[:, 2],
        mode='markers',
        marker=dict(size=6, color='red')
    )
    #
    fig.add_scatter3d(
        x=x_intersection_point,
        y=y_intersection_point,
        z=z_intersection_point,
        mode="markers",
        marker=dict(size=6, color='black',symbol='diamond', opacity=0.7)
    )
    #
    fig.add_mesh3d(
        x=x_ellipsiod_filtered,
        y=y_ellipsiod_filtered,
        z=z_ellipsiod_filtered,
        alphahull=0,
        color='gold',
        opacity=0.4
    )
    #
    #   Customize created figure
    #
    fig.update_layout(
        autosize=False,
        width=2000,
        height=2000,
        #font=dict(family="Times", size=30),
        margin=dict(r=120, l=120, b=120, t=120),
        showlegend=False,
        #title={
        #    'text': "Frame " + frame_number,
        #    'x': 0.5,
        #    'y': 0.93,
        #    'xanchor': 'center',
        #    'yanchor': 'top'
        #}
    )
    #
    fig.update_layout(
        scene=dict(
            aspectmode='data',
            camera=dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=1.20, y=2.0, z=0.7)  # default: eye=dict(x=1.25, y=1.25, z=1.25)
            ),
            xaxis=dict(
                #backgroundcolor="rgb(200, 200, 230)",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="white",
                tickangle=-1,
                tickfont=dict(size=20),
                showaxeslabels=False,
                title=dict(text='X coordinate [Å]', font=dict(family="Times", size=60)),
                color="RebeccaPurple"
            ),
            yaxis=dict(
                #backgroundcolor="rgb(230, 200, 230)",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="white",
                tickangle=50,
                title=dict(text='Y [Å]', font=dict(family="Times", size=40)),
                color="RebeccaPurple"
            ),
            zaxis=dict(
                #backgroundcolor="rgb(230, 230, 200)",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="white",
                tickangle=-90,
                title=dict(text='Z [Å]', font=dict(family="Times", size=40)),
                color="RebeccaPurple"
            )
        )
    )
    #
    #   Interactive and Static Export
    #
    fig.write_html(output_html)
    #
    fig.write_image(output_png, engine="kaleido", scale=2)
    #
    #print("Created output files: ", output_html, output_png)
