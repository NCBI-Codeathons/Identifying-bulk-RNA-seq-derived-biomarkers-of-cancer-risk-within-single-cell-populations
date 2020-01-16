import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

def heatmap(table, autosize=True, width=800, height=1000):
    if type(table.columns) == pd.MultiIndex:
        columns = table.columns.to_series().apply(lambda x: '{0}-{1}'.format(*x))
    else:
        columns = table.columns
    fig = go.Figure(data=go.Heatmap(
                    z=table,
                    x=columns,
                    y=table.index,
                    hoverongaps = False,))
    fig.update_layout(
        title="Gene set enrichment of DE genes by tumor_cluster",
        autosize=autosize,
        width=width,
        height=height,
    )
    fig.show()