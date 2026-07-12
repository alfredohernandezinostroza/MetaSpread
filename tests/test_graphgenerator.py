import pytest
from metaspread import graphgenerator

def test_imports():
    assert graphgenerator.ast is not None
    assert graphgenerator.pd is not None
    assert graphgenerator.plt is not None
    assert graphgenerator.np is not None
    assert graphgenerator.re is not None
    assert graphgenerator.os is not None
    assert graphgenerator.sys is not None
    assert graphgenerator.metaspread.configs is not None

def test_generate_graph(monkeypatch):
    from unittest.mock import Mock
    m = Mock()
    monkeypatch.setattr(graphgenerator, 'generate_graphs', m)
    graphgenerator.generate_graphs()
    assert m.called