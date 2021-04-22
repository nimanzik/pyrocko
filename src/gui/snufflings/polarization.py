import logging
import numpy as num
import matplotlib
from matplotlib import patches
from pyrocko import util, trace, plot
from pyrocko.gui.snuffling import Snuffling, Param, Marker, Choice, Switch, \
    EventMarker

logger = logging.getLogger('pyrocko.gui.snufflings.polarization')

d2r = num.pi / 180.
r2d = 1.0 / d2r


def darken(color, f=0.5):
    return tuple(c*f for c in color[:3]) + color[3:]


class PCAError(Exception):
    pass


class LayoutError(Exception):
    pass


class Polarization(Snuffling):

    u'''
Polarization
============

Investigate patterns of ground motion during the passage of seismic waves.

This Snuffling can be used to analyze and visualize the polarization of seismic
waves from 3-component seismic recordings or to check component orientations of
a seismic sensors when used on signals with known directional properties. The
spatial pattern of ground movement is shown in horizontal and vertical
projections. Principal component analysis and rotation to radial/transverse
components are available as tools. Time window and filter settings can be
interactively adjusted.

Usage
-----

Select one or more normal/phase markers as anchor points for the extraction and
press *Run*. Multiple stations can be selected for direct comparison. Channels
matching the pattern-triplet given in the *Channels* setting are selected for
extraction of a time window around the anchor point. It is assumed that these
three channels correspond to sensor components in the order east, north,
vertical (upward), even if they are named differently.

The time window can be adjusted with the *Length*, *Length Factor*, and
*Offset* parameters. Extracted waveforms are filtered according the *Highpass*
and *Lowpass* parameters (Butterworth 4th order) and demeaned.

To rotate the horizontal components around a vertical axis, use the *\u0394
Azimuth* setting. When station coordinates and an "active" event are available,
the horizontal components can be rotated to radial (away from event) and
transverse (leftwards) orientations using the computed event back-azimuth
(dashed gray line).

If *Show 2D Eigensystems* is selected, a principal component analysis is
performed in each of the shown projects. Red lines are shown to depict the
eigenvectors and the eigenvalues are visualized using ellipse symbols. If *Show
3D Eigensystems* is selected, a principal component analysis is performed using
all three (non-rotated) components and the resulting eigenvectors are depicted
with purple lines.

By default the scaling is automatically adjusted to the maximum value (the
maximum vector length of the three-component signal within the selected time
window). The currently used scaling factor can be frozen by checking *Fix
Scale*.
'''

    def setup(self):
        self.set_name('Polarization')

        self.add_parameter(
            Choice(
                'Channels',
                's_channels',
                '*E, *N, *Z',
                ['*E, *N, *Z', '*1, *2, *Z', '*1, *2, *3', '*R, *T, *Z',
                 'p2, p1, p0', 'HJ1, HJ2, HJ3']))

        self.add_parameter(
            Param('Length [s]', 't_length', 1., 0.001, 1000.))

        self.add_parameter(
            Param('Length Factor', 'f_length', 1., 0., 2.))

        self.add_parameter(
            Param('Offset (relative)', 'ft_offset', 0., -2., 2.))

        self.add_parameter(
            Param(u'\u0394 Azimuth', 'azimuth', 0., -180., 180.))

        self.add_parameter(
            Param(
                'Highpass [Hz]', 'highpass', None, 0.001, 1000.,
                low_is_none=True))

        self.add_parameter(
            Param(
                'Lowpass [Hz]', 'lowpass', None, 0.001, 1000.,
                high_is_none=True))

        self.add_parameter(
            Param('Dot Position', 'dot_position', 0., 0., 1.))

        self.add_parameter(
            Switch(
                'Rotate to RT',
                'rotate_to_rt',
                False))

        self.add_parameter(
            Switch(
                'Show 2D Eigensystems',
                'show_eigensystem_2d',
                False))

        self.add_parameter(
            Switch(
                'Show 3D Eigensystem',
                'show_eigensystem_3d',
                False))

        self.add_parameter(
            Switch(
                'Fix Scale',
                'fix_scale',
                False))

        def new_figure():
            self.setup_figure_frame()
            self.call()

        self.add_trigger('New Figure', new_figure)

        self.fframe = None
        self.iframe = 0
        self.figure_key = None
        self.nsl_to_amax = {}
        self.last_rotate_to_rt = False

        self.font_size = float(matplotlib.rcParams['font.size'])

        self.colors = {
            'fg':  matplotlib.rcParams['axes.edgecolor'],
            'lines': 'black',
            'rot': plot.mpl_color('skyblue1'),
            'pca_2d': plot.mpl_color('scarletred1'),
            'pca_3d': plot.mpl_color('plum2'),
            'backazimuth': plot.mpl_color('aluminium3'),
            'time_window': plot.mpl_color('aluminium2')}

    def panel_visibility_changed(self, visible):
        viewer = self.get_viewer()
        if visible:
            viewer.pile_has_changed_signal.connect(self.adjust_controls)
            self.adjust_controls()

        else:
            viewer.pile_has_changed_signal.disconnect(self.adjust_controls)

    def adjust_controls(self):
        viewer = self.get_viewer()
        dtmin, dtmax = viewer.content_deltat_range()
        maxfreq = 0.5/dtmin
        minfreq = (0.5/dtmax)*0.001
        self.set_parameter_range('lowpass', minfreq, maxfreq)
        self.set_parameter_range('highpass', minfreq, maxfreq)

    def get_selected(self):
        markers = [
            marker for marker in self.get_selected_markers()
            if not isinstance(marker, EventMarker)]

        if not markers:
            self.fail(
                'No selected markers.\n\nCreate and select markers at points '
                'of interest. Normal and phase markers are accepted.')

        d = {}
        for marker in markers:

            tspan = self.transform_time_span(marker.tmin, marker.tmax)
            for nslc in marker.nslc_ids:
                nsl = nslc[:3]
                if nsl in d and d[nsl] != tspan:
                    self.fail(
                        'Inconsistent times for station %s.%s.%s (station '
                        'selected twice?)' % nsl)

                d[nsl] = tspan

        return d

    def transform_time_span(self, tmin, tmax):
        tmin = tmin + self.t_length * self.ft_offset
        tmax = tmin + self.t_length * self.f_length
        return (tmin, tmax)

    def make_selection_markers(self, nsl_to_tspan):
        selection_markers = []
        for nsl in sorted(nsl_to_tspan):
            tmin, tmax = nsl_to_tspan[nsl]
            patterns = self.get_patterns(nsl)

            marker = Marker(
                nslc_ids=patterns, tmin=tmin, tmax=tmax, kind=2)
            selection_markers.append(marker)

        return selection_markers

    def get_selected_channels(self):
        return [x.strip() for x in self.s_channels.split(',')]

    def get_selected_channels_rotated(self):
        chas = self.get_selected_channels()
        chas[0] += "'"
        chas[1] += "'"
        return chas

    def get_patterns(self, nsl):
        patterns = [nsl + (comp,) for comp in self.get_selected_channels()]
        return patterns

    def get_patterns_rotated(self, nsl):
        if self.azimuth == 0.0 and not self.rotate_to_rt:
            return self.get_patterns(nsl)
        else:
            patterns = [
                nsl + (comp,) for comp in self.get_selected_channels_rotated()]
            return patterns

    def get_traces(self, nsl_to_tspan, nsl_to_bazi, tpad):
        if self.highpass is not None:
            tpad_filter = 3.0 / self.highpass
        elif self.lowpass is not None:
            tpad_filter = 3.0 / self.lowpass
        else:
            tpad_filter = 0.0

        # prevent getting insanely long cutouts if e.g. if highpass is still
        # very low, e.g. while moving the slider
        tpad_filter = min(tpad_filter, 5.0 * self.t_length * self.f_length)

        d = {}
        for nsl in sorted(nsl_to_tspan):
            tmin, tmax = nsl_to_tspan[nsl]

            bazimuth = self.azimuth

            if self.rotate_to_rt:
                if nsl not in nsl_to_bazi:
                    self.fail(
                        'Cannot rotate to RT.\n\nStation coordinates must be '
                        'available and an event must be marked as the '
                        '"active" event (select event and press "e").')

                bazimuth += nsl_to_bazi[nsl] + 90.

            patterns = ['.'.join(t) for t in self.get_patterns(nsl)]
            group = []
            for trs in self.get_pile().chopper(
                    tmin=tmin, tmax=tmax, tpad=tpad + tpad_filter,
                    want_incomplete=False,
                    trace_selector=lambda tr: util.match_nslc(
                        patterns, tr.nslc_id)):

                for tr in trs:
                    tr = tr.copy()
                    if self.lowpass is not None \
                            and self.lowpass < 0.5/tr.deltat:
                        tr.lowpass(4, self.lowpass)

                    if self.highpass is not None \
                            and self.highpass < 0.5/tr.deltat:

                        tr.highpass(4, self.highpass)

                    tr.chop(tmin - tpad, tmax + tpad)
                    tr_chop = tr.chop(tmin, tmax, inplace=False)
                    y = tr.get_ydata()
                    tr.set_ydata(y - num.mean(tr_chop.get_ydata()))

                    group.append(tr)

            tr_e = self.get_trace(group, patterns[0])
            tr_n = self.get_trace(group, patterns[1])

            if tr_e and tr_n:
                cha_e = tr_e.channel
                cha_n = tr_n.channel
                group.extend(trace.rotate(
                    group, bazimuth,
                    in_channels=[cha_n, cha_e],
                    out_channels=[cha_n+"'", cha_e+"'"]))

            if group:
                d[nsl] = group

        return d

    def setup_figure_frame(self):
        self.iframe += 1
        self.fframe = self.figure_frame(
            'Particle Motion (%i)' % self.iframe)

        self.fframe.gcf().my_disconnect = None

    def get_figure(self):
        if not self.fframe or self.fframe.closed:
            self.setup_figure_frame()

        return self.fframe.gcf()

    def setup_figure(self, fig, nstations):
        fig.clf()
        new_axes = []

        def iwrap(iy, ix):
            return (ix + iy * 4) + 1

        for istation in range(nstations):
            axes_01 = fig.add_subplot(
                nstations, 4, iwrap(istation, 0), aspect=1.0)
            axes_02 = fig.add_subplot(
                nstations, 4, iwrap(istation, 1), aspect=1.0)
            axes_12 = fig.add_subplot(
                nstations, 4, iwrap(istation, 2), aspect=1.0)
            axes_tr = fig.add_subplot(
                nstations, 4, iwrap(istation, 3))

            for axes in (axes_01, axes_02, axes_12, axes_tr):
                axes.my_stuff = []

                axes.my_line, = axes.plot(
                    [], [], color=self.colors['lines'], lw=1.0)
                axes.my_dot, = axes.plot(
                    [], [], 'o', ms=4, color=self.colors['lines'])

            for axes in (axes_01, axes_02, axes_12):
                axes.get_xaxis().set_tick_params(
                    labelbottom=False, bottom=False)
                axes.get_yaxis().set_tick_params(
                    labelleft=False, left=False)

            axes_tr.get_yaxis().set_tick_params(
                left=False, labelleft=True, length=2.0)

            if istation != nstations - 1:
                axes_tr.get_xaxis().set_tick_params(
                    bottom=False, labelbottom=False)

            lines = []
            dots = []
            for i in range(3):
                lines.append(
                    axes_tr.plot(
                        [], [], color=self.colors['lines'], lw=1.0)[0])
                dots.append(
                    axes_tr.plot(
                        [], [], 'o', ms=4, color=self.colors['lines'])[0])

            axes_tr.my_lines = lines
            axes_tr.my_dots = dots
            axes_tr.my_stuff = []

            new_axes.append(
                (axes_01, axes_02, axes_12, axes_tr))

        def resize_handler(*args):
            self.layout(fig, new_axes)

        if fig.my_disconnect:
            fig.my_disconnect()

        cid_resize = fig.canvas.mpl_connect('resize_event', resize_handler)
        cid_dpi = fig.callbacks.connect('dpi_changed', resize_handler)

        def disconnect():
            fig.canvas.mpl_disconnect(cid_resize)
            fig.callbacks.disconnect(cid_dpi)

        fig.my_disconnect = disconnect

        self.axes = new_axes

    def get_trace(self, traces, pattern):
        trs = [tr for tr in traces if util.match_nslc([pattern], tr.nslc_id)]
        if len(trs) > 1:
            self.fail('Multiple traces matching pattern %s' % pattern)
        elif len(trs) == 0:
            return None
        else:
            return trs[0]

    def get_vector_abs_max(self, traces):
        tr_abs = None
        for tr in traces:
            if tr is not None:
                tr = tr.copy()
                tr.ydata **= 2
                if tr_abs is None:
                    tr_abs = tr
                else:
                    tr_abs.add(tr)

        tr_abs.set_ydata(num.sqrt(tr_abs.ydata))
        return num.max(tr_abs.ydata)

    def set_labels(
            self, istation, nstations, axes_01, axes_02, axes_12, axes_tr):

        chas = self.get_selected_channels()
        if self.azimuth != 0 or self.rotate_to_rt:
            chas_rot = self.get_selected_channels_rotated()
            rcolor = self.colors['rot']
        else:
            chas_rot = chas
            rcolor = self.colors['fg']

        axes_01.set_ylabel(chas[1])
        axes_02.set_ylabel(chas[2])
        axes_12.set_ylabel(chas[2])

        if istation == nstations - 1:
            axes_01.set_xlabel(chas[0])
            axes_02.set_xlabel(chas_rot[0])
            axes_12.set_xlabel(chas_rot[1])
            axes_02.get_xaxis().label.set_color(rcolor)
            axes_12.get_xaxis().label.set_color(rcolor)
            axes_tr.set_xlabel('Time [s]')

        axes_tr.set_yticks([0., 1., 2])
        axes_tr.set_yticklabels([chas_rot[0], chas_rot[1], chas[2]])
        for tlab in axes_tr.get_yticklabels()[:2]:
            tlab.set_color(rcolor)

    def pca(self, trs):

        if any(tr is None for tr in trs):
            raise PCAError('Missing component')

        nss = [tr.data_len() for tr in trs]
        if not all(ns == nss[0] for ns in nss):
            raise PCAError('Traces have different lengths.')

        ns = nss[0]

        if ns < 3:
            raise PCAError('Traces too short.')

        data = num.zeros((ns, len(trs)))
        for itr, tr in enumerate(trs):
            data[:, itr] = tr.ydata

        cov = num.cov(data, rowvar=False)

        evals, evecs = num.linalg.eigh(cov)

        azimuth = r2d*num.arctan2(evecs[1, 1], evecs[0, 1])
        azimuth = ((90. - azimuth) + 180) % 360. - 180.

        return cov, evals, evecs, azimuth

    def draw_cov_ellipse(self, evals, evecs, color, alpha=1.0):
        evals = num.sqrt(evals)
        ell = patches.Ellipse(
            xy=(0.0, 0.0),
            width=evals[0] * 2.,
            height=evals[1] * 2.,
            angle=r2d*num.arctan2(evecs[1, 0], evecs[0, 0]),
            zorder=-10,
            fc=color,
            ec=darken(color),
            alpha=alpha)

        return ell

    def draw(self, groups, nsl_to_tspan, tpad, nsl_to_bazi):

        for insl, nsl in enumerate(sorted(groups)):
            tmin, tmax = nsl_to_tspan[nsl]

            bazimuth = self.azimuth

            if self.rotate_to_rt:
                if nsl not in nsl_to_bazi:
                    self.fail(
                        'Cannot rotate to RT.\n\nActive event must '
                        'available (select event and press "e"). Station '
                        'coordinates must be available.')

                bazimuth += nsl_to_bazi[nsl] + 90.

            for axes in self.axes[insl]:
                while axes.my_stuff:
                    stuff = axes.my_stuff.pop()
                    stuff.remove()

            axes_01, axes_02, axes_12, axes_tr = self.axes[insl]

            axes_01.set_title('.'.join(nsl))

            trs_all = groups[nsl]

            patterns_orig = [
                '.'.join(t) for t in self.get_patterns(nsl)]

            trs_orig = [
                self.get_trace(trs_all, pattern) for pattern in patterns_orig]

            trs_orig_chopped = [
                (tr.chop(tmin, tmax, inplace=False) if tr else None)
                for tr in trs_orig]

            patterns_rot = [
                '.'.join(t) for t in self.get_patterns_rotated(nsl)]

            trs_rot = [
                self.get_trace(trs_all, pattern) for pattern in patterns_rot]

            trs_rot_chopped = [
                (tr.chop(tmin, tmax, inplace=False) if tr else None)
                for tr in trs_rot]

            if self.fix_scale and nsl in self.nsl_to_amax:
                amax = self.nsl_to_amax[nsl]
            else:
                amax = self.get_vector_abs_max(trs_orig_chopped)
                self.nsl_to_amax[nsl] = amax

            for ix, iy, axes, trs in [
                    (0, 1, axes_01, trs_orig_chopped),
                    (0, 2, axes_02, trs_rot_chopped),
                    (1, 2, axes_12, trs_rot_chopped)]:

                axes.set_xlim(-amax*1.05, amax*1.05)
                axes.set_ylim(-amax*1.05, amax*1.05)

                if not (trs[ix] and trs[iy]):
                    continue

                x = trs[ix].get_ydata()
                y = trs[iy].get_ydata()

                axes.my_line.set_data(x, y)
                ipos = int(round(self.dot_position * (x.size-1)))
                axes.my_dot.set_data(x[ipos], y[ipos])

            tref = tmin
            for itr, (tr, tr_chopped) in enumerate(zip(
                    trs_rot, trs_rot_chopped)):

                if tr is None or tr_chopped is None:
                    axes_tr.my_lines[itr].set_data([], [])
                    axes_tr.my_dots[itr].set_data([], [])

                else:
                    y = tr.get_ydata() / (2.*amax) + itr
                    t = tr.get_xdata()
                    t = t - tref

                    ipos = int(round(
                        self.dot_position * (tr_chopped.data_len()-1)))

                    yp = tr_chopped.ydata[ipos] / (2.*amax) + itr
                    tp = tr_chopped.tmin - tref + tr_chopped.deltat*ipos

                    axes_tr.my_lines[itr].set_data(t, y)
                    axes_tr.my_dots[itr].set_data(tp, yp)

            if self.azimuth != 0.0 or self.rotate_to_rt:
                color = self.colors['rot']

                xn = num.sin(bazimuth*d2r)
                yn = num.cos(bazimuth*d2r)
                xe = num.sin(bazimuth*d2r + 0.5*num.pi)
                ye = num.cos(bazimuth*d2r + 0.5*num.pi)

                l1, = axes_01.plot(
                    [0., amax*xn],
                    [0., amax*yn],
                    color=color)

                chas_rot = self.get_selected_channels_rotated()

                font_size = self.font_size

                a1 = axes_01.annotate(
                    chas_rot[1],
                    xy=(amax*xn, amax*yn),
                    xycoords='data',
                    xytext=(-font_size*(xe+.5*xn), -font_size*(ye+.5*yn)),
                    textcoords='offset points',
                    va='center',
                    ha='center',
                    color=color)

                l2, = axes_01.plot(
                    [0., amax*xe],
                    [0., amax*ye],
                    color=color)

                a2 = axes_01.annotate(
                    chas_rot[0],
                    xy=(amax*xe, amax*ye),
                    xycoords='data',
                    xytext=(-font_size*(xn+.5*xe), -font_size*(yn+.5*ye)),
                    textcoords='offset points',
                    va='center',
                    ha='center',
                    color=color)

                axes_01.my_stuff.extend([l1, a1, l2, a2])

            axes_tr.my_stuff.append(axes_tr.axvspan(
                tmin - tref, tmax - tref, color=self.colors['time_window']))

            axes_tr.set_ylim(-1, 3)
            axes_tr.set_xlim(tmin - tref - tpad, tmax - tref + tpad)

            self.set_labels(insl, len(groups), *self.axes[insl])

            if self.show_eigensystem_2d:

                for (ix, iy, axes, trs) in [
                        (0, 1, axes_01, trs_orig_chopped),
                        (0, 2, axes_02, trs_rot_chopped),
                        (1, 2, axes_12, trs_rot_chopped)]:

                    try:
                        cov, evals, evecs, azimuth = self.pca(
                            [trs[ix], trs[iy]])

                        ell = self.draw_cov_ellipse(
                            evals[:2], evecs[:2, :2],
                            color=self.colors['pca_2d'], alpha=0.5)

                        axes.add_artist(ell)
                        axes.my_stuff.append(ell)

                        l1, = axes.plot(
                            [-amax*evecs[0, -1], amax*evecs[0, -1]],
                            [-amax*evecs[1, -1], amax*evecs[1, -1]],
                            color=self.colors['pca_2d'], alpha=0.5)

                        l2, = axes.plot(
                            [-amax*evecs[0, -2], amax*evecs[0, -2]],
                            [-amax*evecs[1, -2], amax*evecs[1, -2]],
                            color=self.colors['pca_2d'], alpha=0.2)

                        axes.my_stuff.extend([l1, l2])

                    except PCAError as e:
                        logger.warn('PCA failed: %s' % e)

            if self.show_eigensystem_3d:
                try:
                    cov, evals, evecs, azimuth = self.pca(trs_orig_chopped)
                    cosa = num.cos(bazimuth*d2r)
                    sina = num.sin(bazimuth*d2r)
                    rot = num.array(
                        [[cosa, -sina, 0.0],
                         [sina, cosa, 0.0],
                         [0.0, 0.0, 1.0]], dtype=num.float)

                    evecs_rot = num.dot(rot, evecs)

                    for (ix, iy, axes, evecs_) in [
                            (0, 1, axes_01, evecs),
                            (0, 2, axes_02, evecs_rot),
                            (1, 2, axes_12, evecs_rot)]:

                        # ell = self.draw_cov_ellipse(
                        #     evals[[ix, iy]], evecs[[ix, iy], [ix, iy]],
                        #     color=self.colors['pca_3d'], alpha=0.5)

                        # axes.add_artist(ell)
                        # axes.my_stuff.append(ell)

                        for (ie, alpha) in [
                                (-1, 0.8),
                                (-2, 0.4),
                                (-3, 0.2)]:

                            lv, = axes.plot(
                                [-amax*evecs_[ix, ie], amax*evecs_[ix, ie]],
                                [-amax*evecs_[iy, ie], amax*evecs_[iy, ie]],
                                color=self.colors['pca_3d'], alpha=alpha)

                            axes.my_stuff.extend([lv])

                except PCAError as e:
                    logger.warn('PCA failed: %s' % e)

            if nsl in nsl_to_bazi:
                l1, = axes_01.plot(
                    [0., amax*num.cos((90. - nsl_to_bazi[nsl])*d2r)],
                    [0., amax*num.sin((90. - nsl_to_bazi[nsl])*d2r)],
                    '--',
                    color=self.colors['backazimuth'])

                axes_01.my_stuff.extend([l1])

    def get_bazis(self):
        event = self.get_viewer().get_active_event()
        if not event:
            return {}

        nsl_to_bazi = dict(
            (station.nsl(), event.azibazi_to(station)[1])
            for station in self.get_stations())

        return nsl_to_bazi

    def layout(self, fig, axes):

        # Do not access self in here. Called from resize in finished plots.

        def get_pixels_factor(fig):
            try:
                r = fig.canvas.get_renderer()
                return 1.0 / r.points_to_pixels(1.0)
            except AttributeError:
                return 1.0

        def rect_to_figure_coords(rect):
            l, b, w, h = rect
            return (l / width, b / height, w / width, h / height)

        ny = len(axes)
        if ny == 0:
            raise LayoutError('No axes given.')

        nx = len(axes[0])

        width, height = fig.canvas.get_width_height()
        pixels = get_pixels_factor(fig)

        margin_left = margin_right = 4. * self.font_size / pixels
        margin_top = margin_bottom = 4. * self.font_size / pixels

        spacing_width = 3. * self.font_size / pixels
        spacing_height = 4. * self.font_size / pixels

        axes_height_avail = height - (ny - 1) * spacing_height \
            - margin_top - margin_bottom

        if axes_height_avail <= 0.0:
            raise LayoutError('Not enough space vertically.')

        axes_width_avail = width - (nx - 1) * spacing_width \
            - margin_left - margin_right

        a_height = axes_height_avail / ny
        a_width = axes_width_avail / (nx + 2)

        a = min(a_height, a_width)

        pad_height = (a_height - a) * ny
        pad_width = (a_width - a) * (nx + 2)

        if axes_width_avail <= 0.0:
            raise LayoutError('Not enough space horizontally.')

        for iy in range(ny):
            y = height - 0.5 * pad_height - margin_top \
                - (iy + 1) * a - iy * spacing_height
            h = a
            for ix in range(nx):
                x = margin_right + 0.5 * pad_width + ix * (a + spacing_width)
                w = a if ix != (nx - 1) else a * 3.0
                axes[iy][ix].set_position(
                    rect_to_figure_coords((x, y, w, h)), which='both')

    def call(self):

        self.cleanup()

        if self.rotate_to_rt != self.last_rotate_to_rt:
            # reset delta azimuth to avoid confusion

            self.set_parameter('azimuth', 0.0)

        self.last_rotate_to_rt = self.rotate_to_rt

        nsl_to_tspan = self.get_selected()
        selection_markers = self.make_selection_markers(nsl_to_tspan)

        self.add_markers(selection_markers)

        nsl_to_bazi = self.get_bazis()

        tpad = self.t_length * self.f_length
        groups = self.get_traces(nsl_to_tspan, nsl_to_bazi, tpad)
        if not groups:
            self.fail(
                'No matching traces. Check time and channel settings. Traces '
                'may not contain gaps within the extracted time window and in '
                'the padding areas left and right. Traces are extracted with '
                'additional padding of 3 x filter corner period to eliminate '
                'artifacts.')

        fig = self.get_figure()

        figure_key = (len(groups), self.iframe)

        if not self.figure_key or self.figure_key != figure_key:
            self.setup_figure(fig, len(groups))
            self.figure_key = figure_key

        self.draw(groups, nsl_to_tspan, tpad, nsl_to_bazi)

        self.layout(fig, self.axes)

        self.fframe.draw()
        tabs = self.fframe.parent().parent()
        # bring plot to front if we are not looking at the markers
        from pyrocko.gui.pile_viewer import PileViewer
        if not isinstance(tabs.currentWidget(), PileViewer):
            tabs.setCurrentWidget(self.fframe)


def __snufflings__():
    return [Polarization()]
