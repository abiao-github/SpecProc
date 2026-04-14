import src.gui.main_window
cls = src.gui.main_window.MainWindow
for name, method in vars(cls).items():
    if hasattr(method, '__code__'):
        if 'logger' in method.__code__.co_varnames:
            print(f"WARNING: logger is local in method: {name}")

for name, obj in vars(src.gui.main_window).items():
    if callable(obj) and hasattr(obj, '__code__'):
        if 'logger' in obj.__code__.co_varnames:
            print(f"WARNING: logger is local in function: {name}")
