from matplotlib.figure import Figure
from matplotlib.pyplot import axes
import numpy as np 
import seaborn as sns
from matplotlib.ticker import AutoMinorLocator

def ecdf(x):
    xs = np.sort(x)
    ys = np.linspace(0,1,len(xs))
    return xs, ys

def plot_ecdf(
    var:str,
    df,
    column_label,
    label_:str,
    ax:axes,
    group_colors:dict,
    label="",
    rasterized=False,
    linewidth=1
    ):
    y = df[(df[column_label] == label_)]
    if y.empty == False:
        xs, ys = ecdf(y[var])
        ax.plot(
            xs,
            ys,
            color=group_colors[y[column_label].iloc[0]],
            linewidth=linewidth,
            label=label,
            rasterized=rasterized
            )

def plot_box(
    df,
    x:str,
    var:str,
    ax:axes,
    fig:Figure,
    group_colors:dict,
    title='',
    x_offset=.2,
    y_offset=.03,
    width_ratio=2.5,
    height_ratio=1.5,
    boxwidth=.6,
    saturation=1,
    fliersize=.5,
    showfliers=True
):
    '''
    Function to overlay box plots to the right of my custom ecdf plots
    '''
    # find plotting position
    pos1 = ax.get_position()
    pos2 = [pos1.x0 + x_offset,
            pos1.y0 + y_offset,
            pos1.width/width_ratio,
            pos1.height/height_ratio]
    # create axes
    ax3 = fig.add_axes(pos2)
    # set up colors
    sns.set_palette(sns.color_palette(group_colors.values()))
    # plot boxes
    g=sns.boxplot(x=x,
                y=var,
                hue=x,
                data=df,
                width=boxwidth,
                ax=ax3,
                saturation=saturation,
                showfliers = showfliers,
                fliersize=fliersize,
                hue_order=group_colors.keys(),
                dodge=False
                )
    # make aesthetics adjustments
    ax3.axes.get_xaxis().set_ticks([])
    g.set(xlabel=None)
    g.set(ylabel=None)
    ax3.yaxis.set_minor_locator(AutoMinorLocator())
    ax3.spines["right"].set_visible(False)
    ax3.spines["left"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    ax3.spines["bottom"].set_visible(False)
    ax3.set_title(title,fontsize=7)
    ax3.get_legend().remove()
    return ax3

def set_equal_axis_range(ax1,ax2):
    axis_x_values = np.hstack(np.array((ax1.get_xlim(),ax2.get_xlim())))
    axis_y_values = np.hstack(np.array((ax1.get_ylim(),ax2.get_ylim())))
    ax1.set_xlim(axis_x_values.min(),axis_x_values.max())
    ax1.set_ylim(axis_y_values.min(),axis_y_values.max())
    ax2.set_xlim(axis_x_values.min(),axis_x_values.max())
    ax2.set_ylim(axis_y_values.min(),axis_y_values.max())