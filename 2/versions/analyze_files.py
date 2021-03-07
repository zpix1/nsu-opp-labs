import glob
import re
import pprint
import collections
import matplotlib.pyplot as plt

def get_time(filename):
    data = open(filename).read()
    try:
        res = re.findall(r'DONE: (\d+.\d+)', data)[0]
    except:
        print('Not found {}'.format(filename))
        exit(1)
    return float(res)

first_times = [
    None,
    get_time(glob.glob('rmv_1_1_1-1.sh.o*')[0]),
    get_time(glob.glob('rmv_2_1_1-1.sh.o*')[0])
]

total = {}

all_files = glob.glob('rmv_*.sh')

for runner in all_files:
    for result in glob.glob('rmv_*.o*'):
        if result.startswith(runner):
            break
    else:
        pass
        print('NOT FOUND:', runner)

for f in glob.glob('rmv_*.o*'):
    # print(f)
    res = re.match(r'rmv_(\d+)_(\d+)_\d+-\d+.sh.o(\d+)', f)
    version = int(res.group(1))
    ncp = int(res.group(2))
    result = get_time(f)
    # print(version, ncp, result)
    t1 = first_times[version]
    Sp = t1 / result
    Ep = Sp / ncp * 100

    total[(version, ncp)] = (result, Sp, Ep)

ans = collections.defaultdict(list)

print(f'version\tC\tT\tA\tE')
for v in sorted(total.keys()):
    version = v[0]
    ncp = v[1]
    print(f'{version}\t{ncp}\t{total[v][0]:.2f}\t{total[v][1]:.2f}\t{total[v][2]:.2f}')
    ans['{}_{}'.format(version, 'PC')].append(ncp)
    ans['{}_{}'.format(version, 'T')].append(total[v][0])
    ans['{}_{}'.format(version, 'Sp')].append(total[v][1])
    ans['{}_{}'.format(version, 'Ep')].append(total[v][2])

if False:
    data = dict(ans)

    def gen_plot(x, y, plot, name, v):
        if name == 'T':
            plot.set_title('Время выполнения, с.' + ' (в. {})'.format(v))
            plot.set_ylim([0, 70])
        if name == 'Sp':
            plot.set_title('Ускорение' + ' (в. {})'.format(v))
            plot.set_ylim([0, 10])
        if name == 'Ep':
            plot.set_title('Эффективность, %' + ' (в. {})'.format(v))
            plot.set_ylim([0, 100])
        plot.plot(x, y)
        plot.set_xticks(x)
        return plot

    fig, plots = plt.subplots(nrows=2, ncols=3)

    i = 0
    for v in ['1', '2']:
        for t in ['T', 'Sp', 'Ep']:
            x = data['{}_{}'.format(v, 'PC')]
            y = data['{}_{}'.format(v, t)]
            gen_plot(x, y, plots.flatten()[i % 3], t, v)
            i += 1
    fig.tight_layout()
    fig.savefig('results_plots.png', dpi=900)