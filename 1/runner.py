from subprocess import check_output
import re


def get_output(cmd):
    return check_output(cmd.split()).decode()

types = [1, 2]
cores = [1, 2, 4]
t1 = None
for t in types:
    print(f'Type {t}')
    for c in cores:
        cmd = f'./run.sh {t} {c}'
        output = get_output(cmd)
        result = float(re.findall('DONE: \d: (.+)\n', output)[0])
        if c == 1:
            t1 = result
        Sp = t1 / result
        Ep = Sp / c * 100
        print(f"{result:.2f} {Sp:.2f} {Ep:.2f}")
