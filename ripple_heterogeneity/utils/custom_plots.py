from matplotlib.figure import Figure
from matplotlib.pyplot import axes
import numpy as np
import seaborn as sns
from matplotlib.ticker import AutoMinorLocator
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch


def ecdf(x):
    xs = np.sort(x)
    ys = np.linspace(0, 1, len(xs))
    return xs, ys


def plot_ecdf(
    var: str,
    df,
    column_label,
    label_: str,
    ax: axes,
    group_colors: dict,
    label="",
    rasterized=False,
    linewidth=1,
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
            rasterized=rasterized,
        )


def plot_box(
    df,
    x: str,
    var: str,
    ax: axes,
    fig: Figure,
    group_colors: dict,
    title="",
    x_offset=0.2,
    y_offset=0.03,
    width_ratio=2.5,
    height_ratio=1.5,
    boxwidth=0.6,
    saturation=1,
    fliersize=0.5,
    showfliers=True,
):
    """
    Function to overlay box plots to the right of my custom ecdf plots
    """
    # find plotting position
    pos1 = ax.get_position()
    pos2 = [
        pos1.x0 + x_offset,
        pos1.y0 + y_offset,
        pos1.width / width_ratio,
        pos1.height / height_ratio,
    ]
    # create axes
    ax3 = fig.add_axes(pos2)
    # set up colors
    sns.set_palette(sns.color_palette(group_colors.values()))
    # plot boxes
    g = sns.boxplot(
        x=x,
        y=var,
        hue=x,
        data=df,
        width=boxwidth,
        ax=ax3,
        saturation=saturation,
        showfliers=showfliers,
        fliersize=fliersize,
        hue_order=list(group_colors.keys()),
        order=list(group_colors.keys()),
        dodge=False,
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
    ax3.set_title(title, fontsize=7)
    ax3.get_legend().remove()
    return ax3


def set_equal_axis_range(ax1, ax2):
    """
    Makes x and y min and max the same between two plots
    """
    axis_x_values = np.hstack(np.array((ax1.get_xlim(), ax2.get_xlim())))
    axis_y_values = np.hstack(np.array((ax1.get_ylim(), ax2.get_ylim())))
    ax1.set_xlim(axis_x_values.min(), axis_x_values.max())
    ax1.set_ylim(axis_y_values.min(), axis_y_values.max())
    ax2.set_xlim(axis_x_values.min(), axis_x_values.max())
    ax2.set_ylim(axis_y_values.min(), axis_y_values.max())


def restore_natural_scale(ax, min_, max_, n_steps=4, x_axis=True, y_axis=True):
    """
    takes x and y ax that are in log10 and puts them into natural scale

    By default, it adjusts both x and y, but you can run this on a single
    axis or two times if you have different scales for x and y
    """
    ticks = np.linspace(min_, max_, n_steps)

    if x_axis:
        ax.set_xticks(ticks)
        ax.set_xticklabels(np.round(10**ticks, 3))

    if y_axis:
        ax.set_yticks(ticks)
        ax.set_yticklabels(np.round(10**ticks, 3))


def plot_events(events, labels, cmap="tab20", gridlines=True, alpha=0.75, ax=None):
    """
    events: nested list of nelpy EpochArrays
    labels: labels related to each event

    example:

        # load sleep states
        state_dict = loading.load_SleepState_states(basepath)

        # make nelpy epoch arrays
        nrem_epochs = nel.EpochArray(state_dict['NREMstate'])
        wake_epochs = nel.EpochArray(state_dict['WAKEstate'])
        rem_epochs = nel.EpochArray(state_dict['REMstate'])

        # add to list
        events = []
        events.append(nrem_epochs)
        events.append(wake_epochs)
        events.append(rem_epochs)

        # plot
        plt.figure(figsize=(20,5))
        plot_events(events,['nrem','wake','rem'])

    Ryan H 2022
    """
    # get colormap
    cmap = matplotlib.cm.get_cmap(cmap)

    # set up y axis
    y = np.linspace(0, 1, len(events) + 1)

    # set up ax if not provided
    if ax is None:
        ax = plt.gca()

    # iter over each event
    for i, evt in enumerate(events):

        # add horizontal line underneath
        if gridlines:
            ax.axhline(y[i] + np.diff(y)[0] / 2, color="k", zorder=-100, alpha=0.1)

        # plot events
        for pair in range(evt.n_intervals):
            ax.axvspan(
                evt.starts[pair],
                evt.stops[pair],
                y[i],
                y[i + 1],
                alpha=alpha,
                color=cmap(i * 0.1),
            )

    ax.set_yticks(y[:-1] + np.diff(y)[0] / 2)
    ax.set_yticklabels(labels)


def plot_ecdf_box(
    data=None,
    x=None,
    hue=None,
    hue_order=None,
    ax=None,
    fig=None,
    x_offset=0.2,
    y_offset=0.03,
    width_ratio=2.5,
    height_ratio=1.5,
    boxwidth=0.6,
    saturation=1,
    fliersize=0.5,
    showfliers=True,
    legend=False,
):
    """
    Code to make a cdf with a box plot off to the right side.

    Example:
    fig, axs = plt.subplots(1,4,
    figsize=functions.set_size("thesis", fraction=1.5, subplots=(1, 4)))
    fig.subplots_adjust(hspace=0.7, wspace=1)
    axs = axs.ravel()

    ax_ = custom_plots.plot_ecdf_box(
        data=temp_df_2,
        x="rankorder",
        hue="ca1_layer",
        hue_order=group_colors.keys(),
        x_offset=0.13,
        fig=fig,
        ax=axs[0]
    )
    """
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()

    sns.ecdfplot(data=data, x=x, hue=hue, hue_order=hue_order, ax=ax, legend=legend)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    pos1 = ax.get_position()
    pos2 = [
        pos1.x0 + x_offset,
        pos1.y0 + y_offset,
        pos1.width / width_ratio,
        pos1.height / height_ratio,
    ]

    # create axes
    ax3 = fig.add_axes(pos2)
    g = sns.boxplot(
        x=hue,
        y=x,
        hue=hue,
        data=data,
        width=boxwidth,
        ax=ax3,
        saturation=saturation,
        showfliers=showfliers,
        fliersize=fliersize,
        hue_order=hue_order,
        order=hue_order,
        dodge=False,
    )
    ax3.axes.get_xaxis().set_ticks([])
    g.set(xlabel=None)
    g.set(ylabel=None)
    ax3.yaxis.set_minor_locator(AutoMinorLocator())
    ax3.spines["right"].set_visible(False)
    ax3.spines["left"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    ax3.spines["bottom"].set_visible(False)
    try:
        ax3.get_legend().remove()
    except:
        pass

    return ax3


def lighten_color(color, amount=0.5):
    """
    Lightens a color by a certain percentage.
    This is useful for adjusting colors for a particular element of a page.
    :param color: The hex color code, e.g. #AABBCC
    :param amount: The amount to lighten the color by.
    :return: The lightened color code in hex, e.g. #FFFFFF
    """
    try:
        c = color.lstrip("#")
        c = tuple(int(c[i : i + 2], 16) for i in (0, 2, 4))
        c = (
            int((1 - amount) * c[0] + amount * 255),
            int((1 - amount) * c[1] + amount * 255),
            int((1 - amount) * c[2] + amount * 255),
        )
        return "#%02x%02x%02x" % c
    except ValueError:
        return color


def adjust_box_widths(g, fac):
    """
    Adjust the widths of a seaborn-generated boxplot.
    """

    # iterating through Axes instances
    for ax in g.axes:

        # iterating through axes artists:
        for c in ax.get_children():

            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5 * (xmin + xmax)
                xhalf = 0.5 * (xmax - xmin)

                # setting new width of box
                xmin_new = xmid - fac * xhalf
                xmax_new = xmid + fac * xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new

                # setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])
