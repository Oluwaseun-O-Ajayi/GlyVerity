"""Tests for the evidence-layer registry."""

import pytest

from glyverity import evidence as ev


class _Dummy(ev.EvidenceLayer):
    name = "dummy_test_layer"
    version = "1"

    def evaluate(self, features):
        return [ev.Evidence(f.feature_id, self.name, "test", 1.0) for f in features]


def teardown_function():
    ev._unregister("dummy_test_layer")


def test_register_and_evaluate():
    ev.register(_Dummy)
    assert "dummy_test_layer" in ev.available_layers()
    layer = ev.get_layer("dummy_test_layer")()
    feats = [ev.Feature("f1", 500.0, 2, 10.0, 1e5)]
    out = layer.evaluate(feats)
    assert out[0].feature_id == "f1"
    assert out[0].layer == "dummy_test_layer"


def test_duplicate_name_rejected():
    ev.register(_Dummy)

    class _Other(_Dummy):
        pass

    with pytest.raises(ValueError):
        ev.register(_Other)


def test_empty_name_rejected():
    class _NoName(ev.EvidenceLayer):
        def evaluate(self, features):
            return []

    with pytest.raises(ValueError):
        ev.register(_NoName)


def test_unknown_layer():
    with pytest.raises(KeyError):
        ev.get_layer("does_not_exist")
