import subprocess

output = subprocess.run(["python","learn_subprocess.py"], capture_output=True)
print()
print(output.stdout)
print()
