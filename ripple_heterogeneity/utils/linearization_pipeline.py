from ripple_heterogeneity.utils import loading
import numpy as np
import matplotlib.pyplot as plt
from track_linearization import make_track_graph
from track_linearization import get_linearized_position

plt.ion()

class NodePicker():
    """Interactive creation of track graph by looking at video frames."""

    def __init__(self, ax=None, video_filename=None, node_color="r", node_size=100):
        if ax is None:
            ax = plt.gca()
        self.ax = ax
        self.canvas = ax.get_figure().canvas
        self.cid = None
        self._nodes = []
        self.node_color = node_color
        self._nodes_plot = ax.scatter(
            [], [],
            zorder=5, s=node_size, color=node_color)
        self.edges = [[]]
        self.video_filename = video_filename

        ax.set_title('Left click to place node.\nRight click to remove node.'
                        '\nShift+Left click to clear nodes.\nCntrl+Left click two nodes to place an edge')
        self.canvas.draw()

        self.connect()

    @property
    def node_positions(self):
        return np.asarray(self._nodes)

    def connect(self):
        if self.cid is None:
            self.cid = self.canvas.mpl_connect(
                'button_press_event', self.click_event)

    def disconnect(self):
        if self.cid is not None:
            self.canvas.mpl_disconnect(self.cid)
            self.cid = None

    def click_event(self, event):
        if not event.inaxes:
            return
        if (event.key not in ['control', 'shift']) & (event.button == 1):  # left click
            self._nodes.append((event.xdata, event.ydata))
        if (event.key not in ['control', 'shift']) & (event.button == 3):  # right click
            self.remove_point((event.xdata, event.ydata))
        if (event.key == 'shift') & (event.button == 1):
            self.clear()
        if (event.key == 'control') & (event.button == 1):
            point = (event.xdata, event.ydata)
            distance_to_nodes = np.linalg.norm(
                self.node_positions - point, axis=1)
            closest_node_ind = np.argmin(distance_to_nodes)
            if len(self.edges[-1]) < 2:
                self.edges[-1].append(closest_node_ind)
            else:
                self.edges.append([closest_node_ind])
            # print(self.edges)
        if (event.key == 'enter'):
            self.disconnect()

        self.redraw()

    def redraw(self):
        # Draw Node Circles
        if len(self.node_positions) > 0:
            self._nodes_plot.set_offsets(self.node_positions)
        else:
            self._nodes_plot.set_offsets([])

        # Draw Node Numbers
        # self.ax.texts = []
        for ind, (x, y) in enumerate(self.node_positions):
            self.ax.text(x, y, ind,
                         zorder=6,
                         fontsize=12,
                         horizontalalignment='center',
                         verticalalignment='center',
                         clip_on=True,
                         bbox=None,
                         transform=self.ax.transData,
                         )
        # Draw Edges
        # self.ax.lines = []  # clears the existing lines
        for edge in self.edges:
            if len(edge) > 1:
                x1, y1 = self.node_positions[edge[0]]
                x2, y2 = self.node_positions[edge[1]]
                self.ax.plot([x1, x2], [y1, y2],
                            color=self.node_color,
                            linewidth=5, zorder=1000)

        self.canvas.draw_idle()
        # self.canvas.draw()

    def remove_point(self, point):
        if len(self._nodes) > 0:
            distance_to_nodes = np.linalg.norm(
                self.node_positions - point, axis=1)
            closest_node_ind = np.argmin(distance_to_nodes)
            self._nodes.pop(closest_node_ind)

    def clear(self):
        self._nodes = []
        self.edges = [[]]
        self.redraw()

def run(basepath):
# video_filename = r"Z:\Data\ORproject\OR15\day1"
    fig, ax = plt.subplots(figsize=(5, 5))

    behave_df = loading.load_animal_behavior(basepath)
    nan_idx = np.isnan(behave_df[["x","y"]]).any(axis=1)

    ax.scatter(behave_df.x, behave_df.y,color="grey",s=1)

    picker = NodePicker(ax=ax, video_filename=basepath)

    track_graph = make_track_graph(picker.node_positions, picker.edges)


    position = np.vstack([behave_df.x.values, behave_df.y.values])[:,~nan_idx].T
    position_df = get_linearized_position(
        position=position,
        track_graph=track_graph,
        edge_order=picker.edges,
    )

    behave_df.loc[~nan_idx, "linear_position"] = position_df.linear_position.values
    behave_df.loc[~nan_idx, "track_segment_id"] = position_df.track_segment_id.values
    behave_df.loc[~nan_idx, "projected_x_position"] = position_df.projected_x_position.values
    behave_df.loc[~nan_idx, "projected_y_position"] = position_df.projected_y_position.values
    
    return behave_df
