import re
import os

with open('git_diff.txt', 'r', encoding='utf-8') as f:
    diff = f.read()

deleted_funcs = []
deleted_classes = []

current_file = None
for line in diff.split('\n'):
    if line.startswith('--- a/'):
        current_file = line[6:]
    elif line.startswith('-def '):
        func_name = re.match(r'-def (\w+)\(', line)
        if func_name:
            deleted_funcs.append((current_file, func_name.group(1)))
    elif line.startswith('-class '):
        class_name = re.match(r'-class (\w+)', line)
        if class_name:
            deleted_classes.append((current_file, class_name.group(1)))

print("Deleted Functions:")
for f in deleted_funcs:
    print(f)
print("\nDeleted Classes:")
for c in deleted_classes:
    print(c)
