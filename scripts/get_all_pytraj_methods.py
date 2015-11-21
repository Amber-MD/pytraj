import inspect
import pytraj as pt

x = dir(pt)
xd = pt.__dict__
y = [_ for _ in x if callable(xd[_]) and
     not _.startswith('_') and
     not inspect.isclass(xd[_])]
print(y)
