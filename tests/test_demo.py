import runpy
from pathlib import Path


def test_demo_runs(capsys):
    demo = Path(__file__).resolve().parents[1] / "examples" / "demo.py"
    runpy.run_path(str(demo), run_name="__main__")
    out = capsys.readouterr().out
    assert "CHRNA7" in out
    assert "无法核验" in out
