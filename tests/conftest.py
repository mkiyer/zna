"""Repo-level pytest options.

``--merge-backend={auto,python}`` selects which merge backend ``auto`` resolves to for
the session -- the suite's second configuration (the merge test-environment notes list
three: compiled, reference, no extension). ``auto`` is the default and changes nothing:
the compiled backend wherever it is built. ``python`` narrows the backend PREFERENCE to
the reference kernel rather than just selecting it once, because ``zna merge`` resolves
``--backend auto`` again on every run and would otherwise re-select the compiled one
mid-session. An explicit ``accel`` -- the cross-backend differentials, and the tests that
shell out with ``--backend accel`` -- still loads the compiled backend, so those keep
comparing the two.

    python -m pytest tests/test_merge.py tests/test_merge_encode.py --merge-backend=python

The third configuration, no extension at all, is a different install (an environment
where the merge extension did not build), not an option: an option cannot un-import a
module that is present.
"""
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--merge-backend", choices=("auto", "python"), default="auto",
        help="merge backend 'auto' resolves to for the session: auto (the compiled one "
             "where built) or python (the reference kernel; explicit 'accel' still "
             "loads the compiled one for the cross-backend tests)")


def pytest_configure(config):
    if config.getoption("--merge-backend") != "python":
        return
    from zna.merge import backend
    backend._PREFERENCE = ("python",)
    backend._default = backend._default_name = None
    backend.use()


def pytest_report_header(config):
    from zna.merge import backend
    return (f"merge backend: --merge-backend={config.getoption('--merge-backend')} -> "
            f"auto resolves to {backend.get_merge_backend_name()}; available "
            f"{backend.available_merge_backends()}")


@pytest.fixture
def merge_backend_option(request):
    """The session's ``--merge-backend`` value."""
    return request.config.getoption("--merge-backend")
