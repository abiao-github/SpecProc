import traceback
import sys

def process():
    try:
        from src.gui.main_window import SpecProcMainWindow
    except ImportError as e:
        print(f"ImportError: {e}")
        return
        
with open("src/gui/main_window.py", "r") as f:
    content = f.read()
    
# Replace string in _on_execution_complete
new_content = content.replace('self.log_text.append(f"Error: {message}")', 'import traceback\n        self.log_text.append(f"Error: {message}\\n{traceback.format_exc()}")')

with open("src/gui/main_window.py", "w") as f:
    f.write(new_content)
