''' input: file with list of generator buses, one per line
    output: new file with list of (bus, cost) pairs, one per line
    defaults to uniform random cost in range [0.5, 2.0]'''

from random import uniform
from sys import argv

def main(filename):
    fin = open(filename)
    fout = open(filename + '.new', 'w')
    for line in fin:
        fout.write('%s %.2f\n' % (line[:-1], uniform(.5, 2)))
    fout.write('END')
    fin.close()
    fout.close()

if __name__ == '__main__':
    if len(argv) != 2:
        print 'usage: python create_costs filename'
    main(argv[1])
