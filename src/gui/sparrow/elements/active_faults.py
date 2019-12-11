# https://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------

from __future__ import absolute_import, print_function, division

import numpy as num

from pyrocko import table, geometry, cake
from pyrocko.guts import Bool, Float
from pyrocko.gui.qt_compat import qw, qc, fnpatch

from pyrocko.dataset.active_faults import ActiveFaults
from pyrocko.gui.vtk_util import ScatterPipe

from .base import Element, ElementState

guts_prefix = 'sparrow'

km = 1e3
COLOR_FAULTS = (0.5, 0.3, .22)


def faults_to_points(faults):
    num_all_nodes = faults.nactive_faults_nodes()
    coords = num.zeros((num_all_nodes, 3))

    x = []
    y = []
    for f in faults.active_faults:
        num_nodes = len(f.lat)
        for i in range(0, num_nodes):
            x.append(f.lat[i])
            y.append(f.lon[i])
    for i in range(0, num_all_nodes):
        coords[i, :] = x[i], y[i], -10*km

    station_table = table.Table()

    station_table.add_col(('coords', '', ('lat', 'lon', 'depth')), coords)

    return geometry.latlondepth2xyz(
        station_table.get_col('coords'),
        planetradius=cake.earthradius)


def faults_to_color(faults):
    colors = []
    num_all_nodes = faults.nactive_faults_nodes()
    for f in faults.active_faults:
        num_nodes = len(f.lat)
        for i in range(0, num_nodes):
            colors.append(COLOR_FAULTS)
    return num.array(colors)


class ActiveFaultsState(ElementState):
    visible = Bool.T(default=True)
    size = Float.T(default=3.0)

    def create(self):
        element = ActiveFaultsElement()
        element.bind_state(self)
        return element


class ActiveFaultsElement(Element):

    def __init__(self):
        Element.__init__(self)
        self._parent = None
        self._state = None
        self._pipe = None
        self._controls = None
        self._active_faults = None
        self._listeners = []

    def bind_state(self, state):
        self._listeners.append(
            state.add_listener(self.update, 'visible'))
        self._listeners.append(
            state.add_listener(self.update, 'size'))
        self._state = state

    def get_name(self):
        return 'Active Faults'

    def set_parent(self, parent):
        self._parent = parent
        if not self._active_faults:
            self._active_faults = ActiveFaults()

        self._parent.add_panel(
            self.get_name(), self._get_controls(), visible=True)
        self.update()

    def update(self, *args):
        state = self._state

        if state.visible:
            if self._pipe is None:
                points = faults_to_points(self._active_faults)
                self._pipe = ScatterPipe(points)

                colors = faults_to_color(self._active_faults)
                self._pipe.set_colors(colors)

            self._pipe.set_size(state.size)
            self._pipe.set_symbol('sphere')
            self._parent.add_actor(self._pipe.actor)

        else:
            self._parent.remove_actor(self._pipe.actor)

        self._parent.update_view()

    def _get_controls(self):
        if self._controls is None:
            from ..state import state_bind_checkbox, state_bind_slider

            frame = qw.QFrame()
            layout = qw.QGridLayout()
            frame.setLayout(layout)

            layout.addWidget(qw.QLabel('Size'), 0, 0)

            slider = qw.QSlider(qc.Qt.Horizontal)
            slider.setSizePolicy(
                qw.QSizePolicy(
                    qw.QSizePolicy.Expanding, qw.QSizePolicy.Fixed))
            slider.setMinimum(0)
            slider.setMaximum(10)
            slider.setSingleStep(0.5)
            slider.setPageStep(1)
            layout.addWidget(slider, 0, 1)
            state_bind_slider(self, self._state, 'size', slider)

            cb = qw.QCheckBox('Show')
            layout.addWidget(cb, 1, 0)
            state_bind_checkbox(self, self._state, 'visible', cb)

            self._controls = frame

        return self._controls


__all__ = [
    'ActiveFaultsElement',
    'ActiveFaultsState'
]
