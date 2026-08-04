# f90wrap conf-file: read '!>' and '!!' comments as documentation.
#
# f90wrap only recognises its own '!%' marker and Doxygen commands (\brief,
# \param, ...).  RABE documents in native style instead, so rewrite the !!
# markers to '!%' before the parser sees them.
#
# Used via:  f90wrap ... --conf-file python/docs_conf.py
#
# Note on the default arguments below: f90wrap exec()s this file from inside a
# function (see f90wrap/scripts/main.py), so a nested def cannot look names up
# in this file's namespace.  Binding them as defaults captures them at def time.

import re

import f90wrap.parser

if not hasattr(f90wrap.parser, "F90File"):
    raise RuntimeError(
        "docs_conf.py expects f90wrap.parser.F90File; this f90wrap version is "
        "not supported. Documentation would be silently dropped, so failing here."
    )


def _init(
    self, fname, _orig=f90wrap.parser.F90File.__init__, _mark=re.compile(r"!(!|>)")
):
    _orig(self, fname)
    self.lines = [_mark.sub("!%", line, count=1) for line in self.lines]


f90wrap.parser.F90File.__init__ = _init
