import plotly.graph_objs as go


def get_black_layout():
    layout = go.Layout(title = "Observed versus Fitted Spectrum",
                       font = dict(
                           color = "white"
                       ),
                       yaxis = dict(
                           title = "Intensity",
                           color = "white",
                           showline = False,
                           zeroline = False
                       ),
                       xaxis = dict(
                           title = "mass/charge [Th]",
                           color = "white",
                       ),
                       plot_bgcolor = "black",
                       paper_bgcolor= "black")
    return layout