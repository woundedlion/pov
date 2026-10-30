"""Shared fake SCons environment for PlatformIO hook tests."""

from pathlib import Path
import re
import types

TOOLS = Path(__file__).resolve().parent


class FakeEnv(dict):
    """The subset of the SCons construction env the hooks touch."""

    def __init__(self, **kw):
        super().__init__(kw)
        self.methods = {}
        self.post_actions = []
        self.pre_actions = []
        self.lib_builders = []
        self.dependencies = []

    def subst(self, s):
        return re.sub(r"\$(\w+)",
                      lambda m: str(self.get(m.group(1), m.group(0))), str(s))

    def Append(self, **kw):
        for key, val in kw.items():
            self.setdefault(key, []).extend(val)

    def Replace(self, **kw):
        self.update(kw)

    def AddMethod(self, fn, name):
        self.methods[name] = fn

    def AddPostAction(self, target, action):
        self.post_actions.append((target, action))

    def AddPreAction(self, target, action):
        self.pre_actions.append((target, action))

    def Depends(self, target, dependency):
        self.dependencies.append((target, dependency))

    def File(self, path):
        return path

    def ConvertInoToCpp(self):
        """PlatformIO's converter: writes CONVERTED_SKETCH beside the sketch."""
        if "CONVERTED_SKETCH" in self:
            sketch = self.methods["FindInoNodes"](self)[0]
            with open(sketch + ".cpp", "w", encoding="utf-8", newline="") as cpp:
                cpp.write(self["CONVERTED_SKETCH"])

    def GetLibBuilders(self):
        return self.lib_builders


def load_hook(name, **scons_globals):
    """Exec a hook the way PlatformIO does: SCons globals injected, not imported."""
    path = TOOLS / name
    mod = types.ModuleType(path.stem)
    mod.__dict__["Import"] = lambda *a, **k: None
    mod.__dict__.update(scons_globals)
    exec(compile(path.read_text(encoding="utf-8"), str(path), "exec"), mod.__dict__)
    return mod
