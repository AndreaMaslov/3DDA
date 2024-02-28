"""
Droplet density analysis from LAMMPS dump data - Plot 2D Angles Esimated density for single frame module
"""

# Packages to install with "pip3 install -U package-name"
#
#   Visualisation
#
#      Plotly static image generation requires Kaleido: pip3 install -U kaleido
import plotly.graph_objects as go

def angles_dens(output_html, output_png, frame_number, frame_x_samples, frame_angles_dens):
    #print(frame_x_samples)
    #print(frame_angles_dens)
    # Use Plotly to plot Angles Esimated density for single frame
    #
    fig = go.Figure()
    #
    fig.add_trace(go.Scatter(
        x=frame_x_samples[:,0],
        y=frame_angles_dens,
        line = dict(color='rgb(117,107,177)', width=6)
        ))
    #
    fig.update_layout(
        autosize=False,
        width=1000,
        height=1000,
        font=dict(family="Times", size=40),
        margin=dict(r=120, l=120, b=120, t=120),
        showlegend=False,
        #title={
        #    'text': "Frame " + frame_number,
        #    'x': 0.5,
        #    'y': 0.93,
        #    'xanchor': 'center',
        #    'yanchor': 'top'
        #},
        xaxis_title='Contact angle [°]',
        yaxis_title='Esimated density'
        )
    fig.update_layout(
        xaxis=dict(
        showline=True,
        showgrid=True,
        showticklabels=True,
        range=[36, 48],
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

