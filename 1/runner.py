from subprocess import check_output
import re
import sys

print('I am alive')

def get_output(cmd):
    return check_output(cmd.split()).decode()

types = [1]
cores = [8]
t1 = None
for t in types:
    print(f'P\tRes\tSp\tEp')
    for c in cores:
        cmd = f'./run.sh {t} {c}'
        output = get_output(cmd)
        print(output, file=sys.stderr)
        result = float(re.findall('DONE: \d+: (.+)\n', output)[0])
        if c == 1:
            t1 = result
        Sp = t1 / result
        Ep = Sp / c * 100
        print(f"{c}\t{result:.2f}\t{Sp:.2f}\t{Ep:.2f}")
