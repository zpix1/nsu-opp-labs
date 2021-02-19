import glob
import re
import pprint
import collections

def get_time(filename):
    data = open(filename).read()
    res = re.findall(r'DONE: \d+: (\d+.\d+)', data)[0]
    return float(res)

first_times = [
    None,
    72.2521,
    72.3259
]

total = {}

all_files = glob.glob('run_mul*.sh')

for runner in all_files:
    for result in glob.glob('run_mul*.o*'):
        if result.startswith(runner):
            break
    else:
        pass
        # print('NOT FOUND:', runner)

for f in glob.glob('run_mul*.o*'):
    # print(f)
    res = re.match(r'run_mul_version_(\d+)_(\d+)_\d+-\d+.sh.o(\d+)', f)
    version = int(res.group(1))
    ncp = int(res.group(2))
    result = get_time(f)
    # print(version, ncp, result)
    t1 = first_times[version]
    Sp = t1 / result
    Ep = Sp / ncp * 100

    total[(version, ncp)] = (result, Sp, Ep)

ans = collections.defaultdict(list)

for v in sorted(total.keys()):
    version = v[0]
    ncp = v[1]
    ans['{}_{}'.format(version, 'PC')].append(ncp)
    ans['{}_{}'.format(version, 'T')].append(total[v][0])
    ans['{}_{}'.format(version, 'Sp')].append(total[v][1])
    ans['{}_{}'.format(version, 'Ep')].append(total[v][2])

print(dict(ans))