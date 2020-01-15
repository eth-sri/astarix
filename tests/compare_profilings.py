import sys
import pandas as pd

def read_tsv(fn):
    df = pd.read_csv(fn, sep='\t')
    df = df.sort_values(by=['readname'])
    df = df.set_index('readname')
    return df

def equal(c1, c2):
    if not c1.equals(c2):
        #set(c1) ^ set(c2)
        print(c1, c2)
        print("Difference: ")
        L = [ '{}: {} != {}'.format(k, v, c2[k]) for k, v in c1.items() if v != c2[k] ]
        print('\n'.join(L))
        return False
    return True

if __name__ == "__main__":
    assert(len(sys.argv) == 3)  # 'two profiling filenames expected'
    a = read_tsv(sys.argv[1])
    b = read_tsv(sys.argv[2])

    print(a.columns)
    if equal(a.index.to_series(), b.index.to_series()):
        if equal(a['spell'], b['spell']):
            if equal(a['cost'], b['cost']):
                print('All good.')
                exit(0)
    exit(1)
