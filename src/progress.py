from __future__ import absolute_import, print_function

import sys
import time

from .get_terminal_size import get_terminal_size


spinner = u'\u25dc\u25dd\u25de\u25df'
skull = u'\u2620'
check = u'\u2714'
bar = u'[- ]'
blocks = u'\u2588\u2589\u258a\u258b\u258c\u258d\u258e\u258f '

ansi_up = u'\033[%iA'
ansi_down = u'\033[%iB'
ansi_left = u'\033[%iC'
ansi_right = u'\033[%iD'
ansi_next_line = u'\033E'

ansi_erase_display = u'\033[2J'
ansi_window = u'\033[%i;%ir'
ansi_move_to = u'\033[%i;%iH'

ansi_clear_down = u'\033[0J'
ansi_clear_up = u'\033[1J'
ansi_clear = u'\033[2J'

ansi_clear_right = u'\033[0K'

ansi_scroll_up = u'\033D'
ansi_scroll_down = u'\033M'

ansi_reset = u'\033c'


g_force_viewer_off = False


class TerminalStatusWindow(object):
    def __init__(self, parent=None):
        self._terminal_size = get_terminal_size()
        self._height = 0
        self._state = 0
        self._parent = parent

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.stop()

    def print(self, s):
        print(s, end='', file=sys.stderr)

    def flush(self):
        print('', end='', flush=True, file=sys.stderr)

    def start(self):
        sx, sy = self._terminal_size
        self._state = 1

    def stop(self):
        if self._state == 1:
            sx, sy = self._terminal_size
            self._resize(0)
            self.print(ansi_move_to % (sy-self._height, 1))
            self.flush()

        self._state = 2
        if self._parent:
            self._parent.hide()

    def _start_show(self):
        sx, sy = self._terminal_size
        self.print(ansi_move_to % (sy-self._height+1, 1))

    def _end_show(self):
        sx, sy = self._terminal_size
        self.print(ansi_move_to % (sy-self._height, 1))
        self.print(ansi_clear_right)

    def _resize(self, height):
        sx, sy = self._terminal_size
        k = height - self._height
        if k > 0:
            self.print(ansi_scroll_up * k)
            self.print(ansi_window % (1, sy-height))
        if k < 0:
            self.print(ansi_window % (1, sy-height))
            self.print(ansi_scroll_down * abs(k))

        self._height = height

    def draw(self, lines):
        if self._state == 0:
            self.start()

        if self._state != 1:
            return

        self._terminal_size = get_terminal_size()
        sx, sy = self._terminal_size
        nlines = len(lines)
        self._resize(nlines)
        self._start_show()

        for iline, line in enumerate(lines):
            if len(line) > sx - 1:
                line = line[:sx-1]

            self.print(ansi_clear_right + line)
            if iline != nlines - 1:
                self.print(ansi_next_line)

        self._end_show()
        self.flush()


class DummyStatusWindow(object):

    def __init__(self, parent=None):
        self._parent = parent

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.stop()

    def stop(self):
        if self._parent:
            self._parent.hide()

    def draw(self, lines):
        pass


class Task(object):
    def __init__(
            self, progress, id, name, n, state='working', logger=None,
            group=None):

        self._id = id
        self._name = name
        self._condition = ''
        self._ispin = 0
        self._i = None
        self._n = n
        self._done = False
        assert state in ('waiting', 'working')
        self._state = state
        self._progress = progress
        self._logger = logger
        self._tcreate = time.time()
        self._group = group

    def __call__(self, it):
        err = False
        try:
            n = 0
            for obj in it:
                self.update(n)
                yield obj
                n += 1

            self.update(n)

        except Exception:
            err = True
            raise

        finally:
            if err:
                self.fail()
            else:
                self.done()

    def log(self, s):
        if self._logger is not None:
            self._logger.info(s)

    def get_group_time_start(self):
        if self._group:
            return self._group.get_group_time_start()
        else:
            return self._tcreate

    def task(self, *args, **kwargs):
        kwargs['group'] = self
        return self._progress.task(*args, **kwargs)

    def update(self, i=None, condition=''):
        self._state = 'working'
        self._condition = condition
        if i is not None:
            if self._n is not None:
                i = min(i, self._n)

            self._i = i

        self._progress._update()

    def done(self, condition=''):
        self.duration = time.time() - self._tcreate

        if self._state in ('done', 'failed'):
            return

        self._condition = condition
        self._state = 'done'
        self._progress._end(self)

    def fail(self, condition=''):
        self.duration = time.time() - self._tcreate

        self._condition = condition
        self._state = 'failed'
        self._progress._end(self)

    def _str_state(self):
        s = self._state
        if s == 'waiting':
            return '  '
        elif s == 'working':
            self._ispin += 1
            return spinner[self._ispin % len(spinner)] + ' '
        elif s == 'done':
            return ''  # check
        elif s == 'failed':
            return skull + ' '
        else:
            return '? '

    def _str_progress(self):
        if self._i is None:
            return self._state
        elif self._n is None:
            return '%i' % self._i
        else:
            nw = len(str(self._n))
            return ('%' + str(nw) + 'i / %i%s') % (
                self._i, self._n,
                '  %3.0f%%' % ((100. * self._i) / self._n)
                if self._state == 'working' else '')

    def _str_condition(self):
        if self._condition:
            return '%s' % self._condition
        else:
            return ''

    def _str_bar(self):
        if self._state == 'working' and self._n and self._i is not None:
            nb = 20
            fb = nb * float(self._i) / self._n
            ib = int(fb)
            ip = int((fb - ib) * (len(blocks)-1))
            s = blocks[0] * ib
            if ib < nb:
                s += blocks[-1-ip] + (nb - ib - 1) * blocks[-1] + blocks[-2]

            # s = ' ' + bar[0] + bar[1] * ib + bar[2] * (nb - ib) + bar[3]
            return s
        else:
            return ''

    def __str__(self):
        return '%s%s: %s' % (
            self._str_state(),
            self._name,
            ' '.join([
                self._str_progress(),
                self._str_bar(),
                self._str_condition(),
                ]))


class Progress(object):

    def __init__(self):
        self._current_id = 0
        self._current_group_id = 0
        self._tasks = {}
        self._tasks_done = []
        self._last_update = 0.0
        self._term = None

    def view(self, viewer='terminal'):
        if g_force_viewer_off:
            viewer = 'off'

        if viewer == 'terminal':
            self._term = TerminalStatusWindow(self)
        elif viewer == 'off':
            self._term = DummyStatusWindow(self)
        else:
            raise ValueError('Invalid viewer choice: %s' % viewer)

        return self._term

    def hide(self):
        self._update(force=True)
        self._term = None

    def task(self, name, n=None, logger=None, group=None):
        self._current_id += 1
        task = Task(
            self, self._current_id, name, n, logger=logger, group=group)
        self._tasks[task._id] = task
        self._update(force=True)
        return task

    def _end(self, task):
        del self._tasks[task._id]
        self._tasks_done.append(task)
        self._update(force=True)

    def _update(self, force=False):
        now = time.time()
        if self._last_update + 0.1 < now or force:
            tasks_done = self._tasks_done
            self._tasks_done = []
            if self._term:
                for task in tasks_done:
                    task.log(str(task))

            lines = self._lines()
            if self._term:
                self._term.draw(lines)

            self._last_update = now

    def _lines(self):
        task_ids = sorted(self._tasks)
        lines = []
        for task_id in task_ids:
            task = self._tasks[task_id]
            lines.extend(str(task).splitlines())

        return lines


progress = Progress()
