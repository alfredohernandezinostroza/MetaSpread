import pytest
from metaspread import interactive

def test_main_menu(monkeypatch):
    from unittest.mock import Mock
    m = Mock()
    monkeypatch.setattr(interactive, 'main_menu', m)
    interactive.main_menu()
    assert m.called