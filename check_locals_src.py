import os
import ast

def check_file(filepath):
    try:
        with open(filepath, 'r') as f:
            tree = ast.parse(f.read())
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                local_names = set()
                for child in ast.walk(node):
                    if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Store):
                        local_names.add(child.id)
                    elif isinstance(child, ast.ImportFrom):
                        for name in child.names:
                            local_names.add(name.asname or name.name)
                    elif isinstance(child, ast.Import):
                        for name in child.names:
                            local_names.add(name.asname or name.name)
                
                # Check if it accesses 'logger' BEFORE assignment, or just if 'logger' is in local_names
                if 'logger' in local_names:
                    print(f"WARNING: logger is a local variable in {filepath} : function {node.name}")
    except Exception as e:
        pass

for root, dirs, files in os.walk('src'):
    for file in files:
        if file.endswith('.py'):
            check_file(os.path.join(root, file))
