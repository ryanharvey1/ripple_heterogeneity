import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import nelpy as nel
import nelpy.plotting as npl
import numpy as np
import functions

def plot_all_replay(
                    bst,
                    spiketrainarray,
                    tuningcurve,
                    tc_placecells,
                    idx=None,
                    title_str=None,
                    vmax=.1,
                    fraction=.25,
                    ratio=(3,1)
                    ):

    if idx is not None:
        bst = bst[idx]
    st = spiketrainarray
    tc = tuningcurve
    tc_placecells = tc_placecells

    no = tc_placecells.get_peak_firing_order_ids()
    st.reorder_units_by_ids(no, inplace=True)

    st_cut = st[bst.support]
    st_cut._support = bst.support # hacky fix so that we can plot events out of order
    st_cut = nel.utils.collapse_time(st_cut)

    # decode neural activity
    posterior, _, _, _ = nel.decoding.decode1D(
                                                bst=bst,
                                                ratemap=tc,
                                                xmin=tc.bins.min(),
                                                xmax=tc.bins.max()
                                                )
    width,height = functions.set_size('thesis', fraction=fraction, subplots=ratio)

    with npl.FigureManager(show=True, figsize=(width*bst.n_bins*.1,height)) as (fig, ax):
        npl.utils.skip_if_no_output(fig)

        pixel_width = 0.5
        if vmax == False:
            npl.imagesc(x=np.arange(bst.n_bins),
                        y=np.arange(np.round(tc.bins.max())),
                        data=posterior,
                        cmap=plt.cm.bone_r,
                        ax=ax)
        else:
            npl.imagesc(x=np.arange(bst.n_bins),
                        y=np.arange(np.round(tc.bins.max())),
                        data=posterior,
                        cmap=plt.cm.bone_r,
                        ax=ax,
                        vmax=vmax)

        npl.utils.no_yticks(ax)

        ax.vlines(np.arange(bst.lengths.sum())-pixel_width,
                    *ax.get_ylim(),
                    lw=1,
                    linestyle=':',
                    color='0.8')

        # ax.vlines(np.cumsum(bst.lengths)-pixel_width, *ax.get_ylim(), lw=1)

        ax.set_xlim(-pixel_width, bst.lengths.sum()-pixel_width)

        event_centers = np.insert(np.cumsum(bst.lengths),0,0)
        event_centers = event_centers[:-1] + bst.lengths/2 - 0.5

        ax.set_xticks(event_centers)
        if idx is not None:
            ax.set_xticklabels(idx)
        else:
            # ax.set_xticklabels(np.arange(bst.n_epochs))
            ax.set_xticklabels('')

        npl.utils.no_xticks(ax)
        npl.utils.clear_top_bottom(ax)

        # ax.set_ylim(np.round(tc.bins.min()),np.round(tc.bins.max()))

        divider = make_axes_locatable(ax)
        axRaster = divider.append_axes("top", size=1, pad=0)

        npl.rasterplot(st_cut, vertstack=True, ax=axRaster, lh=1.25)
        axRaster.set_xlim(st_cut.support.time.squeeze())
        
        bin_edges = np.linspace(st_cut.support.time[0,0],
                                st_cut.support.time[0,1],
                                bst.n_bins+1)

        axRaster.vlines(bin_edges, *ax.get_ylim(), lw=1, linestyle=':', color='0.8')
        # axRaster.vlines(bin_edges[np.cumsum(bst.lengths)], *ax.get_ylim(), lw=1, color='0.2')

        npl.utils.no_xticks(axRaster)
        npl.utils.no_xticklabels(axRaster)
        npl.utils.no_yticklabels(axRaster)
        npl.utils.no_yticks(axRaster)
        ax.set_ylabel('position [cm]')
        ax.set_xlabel('time bins (20 ms)')
        if title_str:
            fig.suptitle(title_str,fontsize=8)
        npl.utils.clear_left_right(axRaster)
        npl.utils.clear_right(axRaster)
        npl.utils.clear_top_bottom(axRaster)
    return ax