with open("src/core/extraction.py", "r") as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    new_lines.append(line)
    if "err = table_hdu.data[err_col] if err_col in table_hdu.data.columns.names else None" in line:
        new_lines.append("                \n")
        new_lines.append("                if len(wav.shape) > 1 and wav.shape[0] == 1:\n")
        new_lines.append("                    wav = wav[0]\n")
        new_lines.append("                if len(flux.shape) > 1 and flux.shape[0] == 1:\n")
        new_lines.append("                    flux = flux[0]\n")
        new_lines.append("                if err is not None and len(err.shape) > 1 and err.shape[0] == 1:\n")
        new_lines.append("                    err = err[0]\n")
        
with open("src/core/extraction.py", "w") as f:
    f.writelines(new_lines)

