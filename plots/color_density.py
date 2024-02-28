"""
Droplet density analysis from LAMMPS dump data - Plot 3D color density module
"""

# Standart python packages
import numpy as np

# Packages to install with "pip3 install -U package-name"
#
#   Visualisation
#
#      Plotly static image generation requires Kaleido: pip3 install -U kaleido
import kaleido #required
import plotly.graph_objects as go


def color_density(output_html, output_png, frame_number, mgrid_x, mgrid_y, mgrid_z, mgrid_density):

    #print("Plotly version: ",plotly.__version__)

    # Use Plotly to plot 3D xyz coordinates and density value as color (without "_nonzero" the plot is useless!)
    #
    #
    # Create a 3D volume trace with colors based on the grid values
    fig = go.Figure(data=go.Volume(
        x=mgrid_x.flatten(),
        y=mgrid_y.flatten(),
        z=mgrid_z.flatten(),
        value=mgrid_density.flatten(),
        opacity=0.7,   # needs to be small to see through all surfaces
        surface_count=60,  # needs to be a large number for good volume rendering
        #colorbar_title="Density (g/cm³)",
        lighting=dict(ambient=0.7),
        colorscale='geyser',
        colorbar=dict(thickness=100, xanchor="right", len=0.5, tickfont=dict(size=40), title=dict(text="[g/cm³]", font=dict(family="Times", size=40))),
        showscale=True,
        caps=dict(x_show=False, y_show=False, z_show=False)
    ))
    #
    #   Customize created figure
    #
    fig.update_layout(
        autosize=False,
        width=2000,
        height=2000,
        font=dict(family="Times", size=18),
        showlegend=False,
        margin=dict(r=120, l=120, b=120, t=120),
        #title={
        #    'text': "Frame " + frame_number,
        #    'x': 0.5,
        #    'y': 0.93,
        #    'xanchor': 'center',
        #    'yanchor': 'top'
        #    }
        )
    #
    fig.update_layout(
        scene=dict(
            aspectmode='data',
            camera=dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=2.0, y=2.0, z=0.7)
                ),
            xaxis=dict(
                #backgroundcolor="rgb(200, 200, 230)",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="white",
                tickangle=-1,
                title=dict(text='X [Å]', font=dict(family="Times", size=40)),
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
