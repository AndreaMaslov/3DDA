"""
Droplet density analysis from LAMMPS dump data - Plot 2D volume value per frame module
"""

# Packages to install with "pip3 install -U package-name"
#
#   Visualisation
#
#      Plotly static image generation requires Kaleido: pip3 install -U kaleido
import plotly.graph_objects as go

def volume_value_per_frame(frames, volume_values, output_html, output_png):

    # Use Plotly to plot most probable angle value per frame
    #
    fig = go.Figure()
    #
    fig.add_trace(go.Scatter(
        x=frames,
        y=volume_values,
        line = dict(color='red', width=2),
        marker = dict(size=0.1)
        ))
    #
    fig.update_layout(
        autosize=False,
        width=1000,
        height=1000,
        margin=dict(r=120, l=120, b=120, t=120),
        showlegend=False,
        title={
            'text': "Droplet Volume Values",
            'x': 0.5,
            'y': 0.93,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        xaxis_title='Frames',
        yaxis_title='Volume [Å^3]'
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
            showgrid=True,
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


    
