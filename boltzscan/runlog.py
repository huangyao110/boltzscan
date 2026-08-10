"""Small PLINK-style command logs for BoltzScan CLI runs."""

from __future__ import annotations

from datetime import datetime
import logging
from pathlib import Path
import platform
import shlex
import socket
import sys
import time

from boltzscan import __version__


class _Tee:
    """Write console text to its original stream and to a command log."""

    def __init__(self, stream, log_handle):
        self.stream = stream
        self.log_handle = log_handle

    def write(self, text):
        self.stream.write(text)
        self.log_handle.write(text)
        self.log_handle.flush()
        return len(text)

    def flush(self):
        self.stream.flush()
        self.log_handle.flush()

    def isatty(self):
        return self.stream.isatty()

    @property
    def encoding(self):
        return self.stream.encoding

    def __getattr__(self, name):
        return getattr(self.stream, name)


class CommandLog:
    """Write one concise, self-contained CLI invocation to ``path``."""

    def __init__(self, path, command, arguments):
        self.path = Path(path)
        self.command = tuple(str(part) for part in command)
        self.arguments = arguments
        self._started = None
        self._handle = None
        self._stdout = None
        self._stderr = None
        self._logging_handler = None

    def _line(self, text=""):
        self._handle.write(text + "\n")
        self._handle.flush()

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.path.open("w", encoding="utf-8")
        self._started = time.monotonic()
        now = datetime.now().astimezone()
        self._line("BoltzScan".center(72, "="))
        self._line(f"Version: {__version__}")
        self._line(f"Start:   {now.isoformat(timespec='seconds')}")
        self._line(f"Host:    {socket.gethostname()} ({platform.system()} {platform.release()})")
        self._line("Command: " + shlex.join(self.command))
        self._line("Options:")
        for name, value in sorted(vars(self.arguments).items()):
            if name != "command":
                self._line(f"  {name} = {value}")
        self._line("Output:")

        self._stdout, self._stderr = sys.stdout, sys.stderr
        sys.stdout = _Tee(self._stdout, self._handle)
        sys.stderr = _Tee(self._stderr, self._handle)

        # ``logging.basicConfig`` may have captured the original stderr before
        # the tee was installed, so attach the same file explicitly.
        self._logging_handler = logging.FileHandler(self.path, mode="a", encoding="utf-8")
        self._logging_handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        logging.getLogger().addHandler(self._logging_handler)
        return self

    def __exit__(self, exc_type, exc, traceback):
        logging.getLogger().removeHandler(self._logging_handler)
        self._logging_handler.close()
        sys.stdout, sys.stderr = self._stdout, self._stderr
        elapsed = time.monotonic() - self._started
        self._line(f"Status:  {'FAILED' if exc_type else 'OK'}")
        if exc is not None:
            self._line(f"Error:   {exc}")
        self._line(f"End:     {datetime.now().astimezone().isoformat(timespec='seconds')}")
        self._line(f"Elapsed: {elapsed:.2f} seconds")
        self._line()
        self._handle.close()
        return False
