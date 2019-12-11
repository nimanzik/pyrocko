# https://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------

from __future__ import absolute_import, print_function, division

import logging

from pyrocko.table import Table, LocationRecipe
from pyrocko.model import Event
from pyrocko.guts_array import Array
from pyrocko.guts import Object, Bool, List, String, load
from pyrocko.geometry import arr_vertices, arr_faces
from pyrocko.gui.qt_compat import qw, qc, fnpatch
from pyrocko.gui.vtk_util import PolygonPipe

from .base import Element, ElementState
from .. import common


logger = logging.getLogger('geometry')

guts_prefix = 'sparrow'

km = 1e3


class Geometry(Object):

    properties = Table.T(default=Table.D(), optional=True)
    vertices = None
    faces = None
    event = Event.T(default=Event.D())
    times = Array.T(
        shape=(None,),
        dtype='float64',
        help='1d vector of times [s] wrt. event time for which '
             'properties have value',
        optional=True)
    stfs = Array.T(
        shape=(None, None),
        dtype='float64',
        help='2d array of source time functions for each patch, number of '
             'columns need tp be equal to number of times',
        optional=True)

    def setup(self, vert_c5, faces):

        self.vertices = Table()
        self.vertices.add_recipe(LocationRecipe())
        self.vertices.add_col((
            'c5', '',
            ('ref_lat', 'ref_lon', 'north_shift', 'east_shift', 'depth')),
            vert_c5)

        self.faces = Table()
        self.faces.add_col('faces', faces)

    def add_property(self, name, values):
        self.properties.add_col(name, values)

    def get_property(self, name):
        return self.properties.get_col(name)


class GeometryState(ElementState):

    visible = Bool.T(default=False)
    geometries = List.T(Geometry.T(), default=[])
    display_parameter = String.T(default='Slip [m]')

    def create(self):
        element = GeometryElement()
        element.bind_state(self)
        return element


class GeometryElement(Element):

    def __init__(self):
        self._listeners = []
        self._parent = None
        self._state = None
        self._controls = None
        self._pipe = []

    def remove(self):
        if self._parent and self._state:
            self._parent.state.elements.remove(self._state)

    def remove_pipes(self):
        for pipe in self._pipe:
            if isinstance(pipe.actor, list):
                for act in pipe.actor:
                    self._parent.remove_actor(act)
            else:
                self._parent.remove_actor(pipe.actor)
        self._pipe = []

    def set_parent(self, parent):
        self._parent = parent
        self._parent.add_panel(
            self.get_name(), self._get_controls(), visible=True)
        self.update()

    def unset_parent(self):
        self.unbind_state()
        if self._parent:
            if self._pipe:
                self.remove_pipes()

            if self._controls:
                self._parent.remove_panel(self._controls)
                self._controls = None

            self._parent.update_view()
            self._parent = None

    def bind_state(self, state):
        upd = self.update
        self._listeners.append(upd)
        state.add_listener(upd, 'visible')
        state.add_listener(upd, 'geometries')
        state.add_listener(upd, 'display_parameter')
        self._state = state

    def unbind_state(self):
        for listener in self._listeners:
            try:
                listener.release()
            except Exception as e:
                pass

        self._listeners = []
        self._state = None

    def get_name(self):
        return 'Geometry'

    def open_file_load_dialog(self):
        caption = 'Select one file containing a geometry to open'
        fns, _ = fnpatch(qw.QFileDialog.getOpenFileNames(
            self._parent, caption, options=common.qfiledialog_options))

        if fns:
            self.load_file(str(fns[0]))
        else:
            return

    def load_file(self, path):

        loaded_geometry = load(filename=path)

        self._parent.remove_panel(self._controls)
        self._controls = None
        self._state.geometries.append(loaded_geometry)
        self._parent.add_panel(
            self.get_name(), self._get_controls(), visible=True)

        self.update()

    def update(self, *args):

        state = self._state

        if state.visible:
            for geo in state.geometries:

                vertices = arr_vertices(geo.vertices.get_col('xyz'))
                faces = arr_faces(geo.faces.get_col('faces'))
                values = geo.get_property(state.display_parameter)

                self._pipe.append(
                    PolygonPipe(
                        vertices, faces,
                        values=values, contour=False,
                        cbar_title=state.display_parameter))

                if isinstance(self._pipe[-1].actor, list):
                    self._parent.add_actor_list(self._pipe[-1].actor)
                else:
                    self._parent.add_actor(self._pipe[-1].actor)

        else:
            if self._pipe:
                self.remove_pipes()

        self._parent.update_view()

    def _get_controls(self):
        if not self._controls:
            from ..state import state_bind_checkbox, state_bind_combobox

            frame = qw.QFrame()
            layout = qw.QGridLayout()
            frame.setLayout(layout)

            # load geometrie
            pb = qw.QPushButton('Load')
            layout.addWidget(pb, 0, 0)
            pb.clicked.connect(self.open_file_load_dialog)

            # visibility
            cb = qw.QCheckBox('Show')
            layout.addWidget(cb, 2, 0)
            state_bind_checkbox(self, self._state, 'visible', cb)

            # property choice
            layout.addWidget(qw.QLabel('Display parameter'), 1, 0)
            cb = qw.QComboBox()
            if self._state.geometries:
                props = []
                for geom in self._state.geometries:
                    for prop in geom.propertoes.get_col_names():
                        props.append(prop)

                unique_props = list(set(props))
                for i, s in enumerate(unique_props):
                    cb.insertItem(i, s)

            layout.addWidget(cb, 1, 1)
            state_bind_combobox(self, self._state, 'display_parameter', cb)

            self._controls = frame

        return self._controls


__all__ = [
    'GeometryElement',
    'GeometryState',
    'Geometry'
]