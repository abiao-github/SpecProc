import src.core.extraction
func = src.core.extraction.process_extraction_stage
print("Local variables in extraction:", func.__code__.co_varnames)
if 'logger' in func.__code__.co_varnames:
    print("WARNING: logger is a local variable in extraction!")
else:
    print("logger is NOT a local variable in extraction.")
    
for name, obj in vars(src.core.extraction).items():
    if callable(obj) and hasattr(obj, '__code__'):
        if 'logger' in obj.__code__.co_varnames:
            print(f"WARNING: logger is local in {name}")
