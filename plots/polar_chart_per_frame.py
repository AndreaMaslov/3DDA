"""
Droplet density analysis from LAMMPS dump data - Plot 2D polar chart module
"""

# Standart python packages
import numpy as np

# Packages to install with "pip3 install -U package-name"
#
#   Visualisation
#
#      Plotly static image generation requires Kaleido: pip3 install -U kaleido
import plotly.graph_objects as go

def polar_chart(output_html, output_png, frame_number, angle_per_frame_mprob_value, lst_tangent_angles_per_frame, angle_order_xy_intersection_points):

    # Create ndarray of tangent_angles and angle_order_xy_intersection_points
    #
    angles_intersection_order_per_frame = np.transpose([lst_tangent_angles_per_frame, angle_order_xy_intersection_points])
    #print("Statistics for angles_intersection_order_per_frame ndarray:")
    #print("    Number of axes (dimensions): ", angles_intersection_order_per_frame.ndim)
    #print("    Number of elements in each dimension: ", angles_intersection_order_per_frame.shape)
    #print("    All angles_intersection_order_per_frame elements: ", angles_intersection_order_per_frame)
    #
    angles_intersection_order_per_frame_sorted = angles_intersection_order_per_frame[np.argsort(angles_intersection_order_per_frame[:, 1])]
    #print("Statistics for angles_intersection_order_per_frame_sorted ndarray:")
    #print("    Number of axes (dimensions): ", angles_intersection_order_per_frame_sorted.ndim)
    #print("    Number of elements in each dimension: ", angles_intersection_order_per_frame_sorted.shape)
    #print("    All angles_intersection_order_per_frame_sorted elements: ", angles_intersection_order_per_frame_sorted)
    #
    # Append to the end of angles_intersection_order_per_frame_sorted ndarray it's first element: needed to close line in Polar Charts
    first_element = angles_intersection_order_per_frame_sorted[0]
    first_element = first_element.reshape(1,2)
    angles_intersection_order_per_frame_sorted = np.append(angles_intersection_order_per_frame_sorted, first_element, axis=0)
    #print("Statistics for angles_intersection_order_per_frame_sorted ndarray:")
    #print("    Number of axes (dimensions): ", angles_intersection_order_per_frame_sorted.ndim)
    #print("    Number of elements in each dimension: ", angles_intersection_order_per_frame_sorted.shape)
    #print("    All angles_intersection_order_per_frame_sorted elements: ", angles_intersection_order_per_frame_sorted)

    # Use Plotly to plot Polar Chart of most probable angle value
    #
    fig = go.Figure()
    #
    fig.add_trace(go.Scatterpolar(
        name = "CA",
        r = angles_intersection_order_per_frame_sorted[:,0],      # angles
        theta = angles_intersection_order_per_frame_sorted[:,1],  # intersection point's azimuths
        mode = 'markers+lines',   #  mode = 'markers'
        fill= 'toself',   #  fill= 'none',
        marker = dict(size=12, color="black", opacity=0.7),  #  marker = dict(size=6, color="royalblue"),
        line = dict(width=2, color="gold", dash="solid", shape="spline", smoothing=0)
    ))
    #
    # Define a circle of most probable angle value (from KDE)
    r, theta = circle_polar(0, 0, angle_per_frame_mprob_value)
    fig.add_trace(go.Scatterpolar(
        name = "Spherical CA",
        r=r,
        theta=theta,
        mode='lines',
        line = dict(width=10, color="MediumPurple", dash="solid", shape="spline", smoothing=1),
        opacity = 0.9
    ))
    
    fig.update_layout(
        polar = dict(
            radialaxis = dict(range = [10, 60],
                              showticklabels = True,
                              showline=False,
                              ticksuffix="°",
            tickfont=dict(
                family='Times',
                size=36,
                color='red'
                ),
            ),
            angularaxis = dict(showgrid = False,
                               linecolor='black',
                               linewidth= 6,
                               showticklabels = False,
                                ticks = '')
        ),
        showlegend = True,
        legend=dict(
        yanchor="top",
        y=1.2,
        xanchor="left",
        x=0.001
)
    )
    #
    fig.update_layout(
        polar_radialaxis_gridcolor="black",
        polar_angularaxis_gridcolor="white",
        polar_bgcolor='white',
        paper_bgcolor='white'
    )
    
    fig.update_layout(
        autosize=False,
        width=1000,
        height=1000,
        font=dict(family="Times", size=40),
        margin=dict(r=120, l=120, b=120, t=120),
        #title={
        #    'text': "Frame " + frame_number,
        #    'x': 0.5,
        #    'y': 0.93,
        #    'xanchor': 'center',
        #    'yanchor': 'top'
        #    }
        )
    #   Interactive and Static Export
    #
    fig.write_html(output_html)
    #
    fig.write_image(output_png, engine="kaleido", scale=2)
    #
    #print("Created output files: ", output_html, output_png)


def circle_polar(r0, theta0, rho, degrees = True):
    # compute the polar coordinates for 100 points on a circle of center (r0, theta0) and radius rho
    if degrees:
        theta0 = theta0 * np.pi / 180
    
    phi = np.linspace(0, 2*np.pi, 100)
    x = rho * np.cos(phi) + r0 * np.cos(theta0)
    y = rho * np.sin(phi) + r0 * np.sin(theta0)
    r = np.sqrt(x**2 + y**2)
    
    J = np.where(y<0)
    theta  = np.arctan2(y, x)
    theta[J]= theta[J] + 2*np.pi

    return (r, theta * 180 / np.pi) if degrees else (r, theta)

