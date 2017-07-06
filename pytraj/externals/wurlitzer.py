# Note: This code is copied from
# https://github.com/minrk/wurlitzer
# git tag: e9549c29ead108cbe4e7cf5bf0cbe0dcbc1ce9e0

# The MIT License (MIT)
#
# Copyright (c) 2016 Min RK
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Capture C-level FD output on pipes

Use `wurlitzer.pipes` or `wurlitzer.sys_pipes` as context managers.
"""
from __future__ import print_function

__version__ = '0.2.1.dev'

__all__ = [
    'pipes',
    'sys_pipes',
    'sys_pipes_forever',
    'stop_sys_pipes',
    'Wurlitzer',
]

from contextlib import contextmanager
import ctypes
from fcntl import fcntl, F_GETFL, F_SETFL
import io
import os
import select
import sys
import threading
import warnings

libc = ctypes.CDLL(None)

try:
    c_stdout_p = ctypes.c_void_p.in_dll(libc, 'stdout')
    c_stderr_p = ctypes.c_void_p.in_dll(libc, 'stderr')
except ValueError:  # pragma: no cover
    # libc.stdout is has a funny name on OS X
    c_stdout_p = ctypes.c_void_p.in_dll(libc, '__stdoutp')  # pragma: no cover
    c_stderr_p = ctypes.c_void_p.in_dll(libc, '__stderrp')  # pragma: no cover

STDOUT = 2
PIPE = 3

_default_encoding = getattr(sys.stdin, 'encoding', None) or 'utf8'
if _default_encoding.lower() == 'ascii':
    # don't respect ascii
    _default_encoding = 'utf8'  # pragma: no cover


class Wurlitzer(object):
    """Class for Capturing Process-level FD output via dup2
    
    Typically used via `wurlitzer.capture`
    """
    flush_interval = 0.2

    def __init__(self, stdout=None, stderr=None, encoding=_default_encoding):
        """
        Parameters
        ----------
        stdout: stream or None
            The stream for forwarding stdout.
        stderr = stream or None
            The stream for forwarding stderr.
        encoding: str or None
            The encoding to use, if streams should be interpreted as text.
        """
        self._stdout = stdout
        if stderr == STDOUT:
            self._stderr = self._stdout
        else:
            self._stderr = stderr
        self.encoding = encoding
        self._save_fds = {}
        self._real_fds = {}
        self._handlers = {}
        self._handlers['stderr'] = self._handle_stderr
        self._handlers['stdout'] = self._handle_stdout

    def _setup_pipe(self, name):
        real_fd = getattr(sys, '__%s__' % name).fileno()
        save_fd = os.dup(real_fd)
        self._save_fds[name] = save_fd

        pipe_out, pipe_in = os.pipe()
        os.dup2(pipe_in, real_fd)
        os.close(pipe_in)
        self._real_fds[name] = real_fd

        # make pipe_out non-blocking
        flags = fcntl(pipe_out, F_GETFL)
        fcntl(pipe_out, F_SETFL, flags | os.O_NONBLOCK)
        return pipe_out

    def _decode(self, data):
        """Decode data, if any
        
        Called before pasing to stdout/stderr streams
        """
        if self.encoding:
            data = data.decode(self.encoding, 'replace')
        return data

    def _handle_stdout(self, data):
        if self._stdout:
            self._stdout.write(self._decode(data))

    def _handle_stderr(self, data):
        if self._stderr:
            self._stderr.write(self._decode(data))

    def _setup_handle(self):
        """Setup handle for output, if any"""
        self.handle = (self._stdout, self._stderr)

    def _finish_handle(self):
        """Finish handle, if anything should be done when it's all wrapped up."""
        pass

    def __enter__(self):
        # flush anything out before starting
        libc.fflush(c_stdout_p)
        libc.fflush(c_stderr_p)
        # setup handle
        self._setup_handle()

        # create pipe for stdout
        pipes = []
        names = {}
        if self._stdout:
            pipe = self._setup_pipe('stdout')
            pipes.append(pipe)
            names[pipe] = 'stdout'
        if self._stderr:
            pipe = self._setup_pipe('stderr')
            pipes.append(pipe)
            names[pipe] = 'stderr'

        def forwarder():
            """Forward bytes on a pipe to stream messages"""
            while True:
                # flush libc's buffers before calling select
                libc.fflush(c_stdout_p)
                libc.fflush(c_stderr_p)
                r, w, x = select.select(pipes, [], [], self.flush_interval)
                if not r:
                    # nothing to read, next iteration will flush and check again
                    continue
                for pipe in r:
                    name = names[pipe]
                    data = os.read(pipe, 1024)
                    if not data:
                        # pipe closed, stop polling
                        pipes.remove(pipe)
                    else:
                        handler = getattr(self, '_handle_%s' % name)
                        handler(data)
                if not pipes:
                    # pipes closed, we are done
                    break

        self.thread = threading.Thread(target=forwarder)
        self.thread.daemon = True
        self.thread.start()

        return self.handle

    def __exit__(self, exc_type, exc_value, traceback):
        # flush the underlying C buffers
        libc.fflush(c_stdout_p)
        libc.fflush(c_stderr_p)
        # close FDs, signaling output is complete
        for real_fd in self._real_fds.values():
            os.close(real_fd)
        self.thread.join()

        # restore original state
        for name, real_fd in self._real_fds.items():
            save_fd = self._save_fds[name]
            os.dup2(save_fd, real_fd)
            os.close(save_fd)
        # finalize handle
        self._finish_handle()


@contextmanager
def pipes(stdout=PIPE, stderr=PIPE, encoding=_default_encoding):
    """Capture C-level stdout/stderr in a context manager.
    
    The return value for the context manager is (stdout, stderr).
    
    Examples
    --------
    
    >>> with capture() as (stdout, stderr):
    ...     printf("C-level stdout")
    ... output = stdout.read()
    """
    stdout_pipe = stderr_pipe = False
    # setup stdout
    if stdout == PIPE:
        stdout_r, stdout_w = os.pipe()
        stdout_w = os.fdopen(stdout_w, 'wb')
        if encoding:
            stdout_r = io.open(stdout_r, 'r', encoding=encoding)
        else:
            stdout_r = os.fdopen(stdout_r, 'rb')
        stdout_pipe = True
    else:
        stdout_r = stdout_w = stdout
    # setup stderr
    if stderr == STDOUT:
        stderr_r = None
        stderr_w = stdout_w
    elif stderr == PIPE:
        stderr_r, stderr_w = os.pipe()
        stderr_w = os.fdopen(stderr_w, 'wb')
        if encoding:
            stderr_r = io.open(stderr_r, 'r', encoding=encoding)
        else:
            stderr_r = os.fdopen(stderr_r, 'rb')
        stderr_pipe = True
    else:
        stderr_r = stderr_w = stderr
    if stdout_pipe or stderr_pipe:
        capture_encoding = None
    else:
        capture_encoding = encoding
    w = Wurlitzer(stdout=stdout_w, stderr=stderr_w, encoding=capture_encoding)
    try:
        with w:
            yield stdout_r, stderr_r
    finally:
        # close pipes
        if stdout_pipe:
            stdout_w.close()
        if stderr_pipe:
            stderr_w.close()


def sys_pipes(encoding=_default_encoding):
    """Redirect C-level stdout/stderr to sys.stdout/stderr
    
    This is useful of sys.sdout/stderr are already being forwarded somewhere.
    
    DO NOT USE THIS if sys.stdout and sys.stderr are not already being forwarded.
    """
    return pipes(sys.stdout, sys.stderr, encoding=encoding)


_mighty_wurlitzer = None


def sys_pipes_forever(encoding=_default_encoding):
    """Redirect all C output to sys.stdout/err
    
    This is not a context manager; it turns on C-forwarding permanently.
    """
    global _mighty_wurlitzer
    if _mighty_wurlitzer is None:
        _mighty_wurlitzer = sys_pipes(encoding)
    _mighty_wurlitzer.__enter__()


def stop_sys_pipes():
    """Stop permanent redirection started by sys_pipes_forever"""
    global _mighty_wurlitzer
    if _mighty_wurlitzer is not None:
        _mighty_wurlitzer.__exit__(None, None, None)
        _mighty_wurlitzer = None


def load_ipython_extension(ip):
    """Register me as an IPython extension
    
    Captures all C output during execution and forwards to sys.
    
    Does nothing on terminal IPython.
    
    Use: %load_ext wurlitzer
    """
    if not getattr(ip, 'kernel'):
        warnings.warn(
            "wurlitzer extension doesn't do anything in terminal IPython")
        return
    ip.events.register('pre_execute', sys_pipes_forever)
    ip.events.register('post_execute', stop_sys_pipes)


def unload_ipython_extension(ip):
    """Unload me as an IPython extension
    
    Use: %unload_ext wurlitzer
    """
    if not getattr(ip, 'kernel'):
        return
    ip.events.unregister('pre_execute', sys_pipes_forever)
    ip.events.unregister('post_execute', stop_sys_pipes)
