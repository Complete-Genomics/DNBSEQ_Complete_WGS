import subprocess, sys

template, reports = sys.argv[1], sys.argv[2:]
reports = " ".join(reports)

cmd = f"""
paste -d "," {template} {reports}
"""

process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

stdout, stderr = process.communicate()

if process.returncode != 0:
    print(f"Error running command:\n{stderr.decode().strip()}")
else:
    # Get the number of heterozygous SNPs from the command output
    counts = stdout.decode()
    print(counts)