"""
Droplet density analysis from LAMMPS dump data - Plot 2D most probable angle value per frame module
"""

# Packages to install with "pip3 install -U package-name"
#
#   Visualisation
#
#      Plotly static image generation requires Kaleido: pip3 install -U kaleido
import plotly.graph_objects as go

def tangent_angle_value_per_frame(frames, angles, output_html, output_png, mode):

    # Use Plotly to plot most probable angle value per frame
    #
    fig = go.Figure()
    #
    fig.add_trace(go.Scatter(
        x=frames,
        y=angles,
        mode=mode,
        line = dict(color='rgb(117,107,177)', width=6),
        marker = dict(color='black', size=10)
        )
    )
    #
    fig.update_layout(
        autosize=False,
        width=1000,
        height=1000,
        font=dict(family="Times", size=40),
        margin=dict(r=120, l=120, b=120, t=120),
        showlegend=False,
        #title={
        #    'text': "Water Contact Angle (Tangent)",
        #    'x': 0.5,
        #    'y': 0.93,
        #    'xanchor': 'center',
        #    'yanchor': 'top'
        #},
        xaxis_title='Frames',
        yaxis_title='Theta [°]'
        )
    fig.update_layout(
        xaxis=dict(
        showline=True,
        showgrid=True,
        showticklabels=True,
        linecolor='black',
        linewidth=4
        ),
        yaxis=dict(
            showline=True,
            showgrid=False,
            gridcolor='LightPink',
            showticklabels=True,
            linecolor='black',
            linewidth=4,
            ticks='outside'
        ),
        plot_bgcolor='white'
    )
    #   Interactive and Static Export
    #
    fig.write_html(output_html)
    #
    fig.write_image(output_png, engine="kaleido", scale=2)
    #
    #print("Created output files: ", output_html, output_png)


    
