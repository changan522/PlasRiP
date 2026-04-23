from __future__ import annotations

import plasrip


def test_package_version_exists():
    assert isinstance(plasrip.__version__, str)
    assert plasrip.__version__


def test_public_api_exports():
    assert callable(plasrip.run_host_prediction)
    assert callable(plasrip.run_ptu_assignment)
    assert callable(plasrip.run_core_analysis)
    assert callable(plasrip.run_annotation)
    assert callable(plasrip.run_pipeline)
