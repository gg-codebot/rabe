import inspect
import re

from rabe.fieldline_mod import FlockOfFieldlines

DOC = FlockOfFieldlines.__init__.__doc__
PARAM_HEADER = re.compile(r"^(\w+) : \w")

# for testing
FAULTY_DOC = """
self = Flock_Of_Fieldlines_T(iota, nfp, third_argument)
Defined at ...

Parameters
----------
iota : float64
  rotational transform of field lines

nfp : float64
  number of field periods; must be positive integer

third_argument : float64

Returns
-------
flock : Flock_Of_Fieldlines_T
"""


def split_sections(doc):
    """Return (summary, parameter block) of a numpydoc-style docstring.

    f90wrap indents the docstring to the level of the generated `def`, so
    strip that common indentation before matching line starts.
    """
    summary, _, rest = inspect.cleandoc(doc).partition("Parameters")
    return summary, rest.partition("Returns")[0]


def test_summary_has_prose(doc, min_words=20):
    has_failed = False
    summary, _ = split_sections(doc)
    prose = [
        s
        for line in summary.splitlines()
        if (s := line.strip())
        and not s.startswith("Defined at")
        and not s.startswith("self = ")
    ]
    n_words = len(" ".join(prose).split())
    if not n_words >= min_words:
        print(
            f"summary has {n_words} words of prose, expected >= {min_words}:\n{summary}"
        )
        has_failed = True
    return has_failed


def test_every_parameter_is_described(doc, n_params=7, min_words=3):
    has_failed = False
    _, params = split_sections(doc)
    described = {}
    name = None
    for line in params.splitlines():
        if match := PARAM_HEADER.match(line):
            name = match.group(1)
            described[name] = ""
        elif name and line.startswith(" "):
            described[name] += " " + line.strip()

    if not len(described) == n_params:
        print(
            f"found {len(described)} parameters {list(described)}, expected {n_params}"
        )
        has_failed = True

    undescribed = [n for n, desc in described.items() if len(desc.split()) < min_words]
    if undescribed:
        print(f"parameters described in fewer than {min_words} words: {undescribed}")
        has_failed = True

    return has_failed


if __name__ == "__main__":
    doc = DOC
    print("--------------------------------------")
    print("Sample testing generated docstring ...")
    has_no_prose = test_summary_has_prose(doc)
    has_no_parameter = test_every_parameter_is_described(doc)
    assert not (
        has_no_prose or has_no_parameter
    ), "Docstring test failed - see output above"
    print("Docstring test passed!")
    print("--------------------------------------")
