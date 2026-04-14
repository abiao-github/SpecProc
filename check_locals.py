import src.gui.main_window
func = src.gui.main_window.MainWindow._execute_selected_stages
print("Local variables:", func.__code__.co_varnames)
if 'logger' in func.__code__.co_varnames:
    print("WARNING: logger is a local variable!")
else:
    print("logger is NOT a local variable.")
