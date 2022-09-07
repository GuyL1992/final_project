import os
from subprocess import Popen, PIPE, STDOUT

def read_file(file_name):

    with open(file_name,'r') as file:

            def parse_row(row):
                return list(map(lambda point: float(point),row.split(',')))

            observations = list(map(lambda row: parse_row(row), file))
            return observations

os.system("python3 setup.py build_ext --inplace")



BASE_PATH = "tests/"
goal = input("goal: ")
num = input("test number: ")
test = f"spk_{num}.txt" if goal != "jacobi" else f"jacobi_{num}.txt"

os.system(f"python3 spkmeans.py 3 {goal} {BASE_PATH}{test}")

input("copy result")

path = f"{BASE_PATH}output.txt"
test_path = f"{BASE_PATH}{goal}/{test}"

result = read_file(path)
test = read_file(test_path)

print (result == test)


