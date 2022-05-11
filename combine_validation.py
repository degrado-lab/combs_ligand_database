import os
import pickle
import argparse
import pandas as pd

def parse_args():
    argp = argparse.ArgumentParser()
    argp.add_argument('validation_dir', type=os.path.realpath, 
                      help="Path to directory at which pickled validation "
                      "dataframes are stored.")
    argp.add_argument('outpath', type=os.path.realpath, 
                      help="Path at which to output the final combined "
                      "validation dataframe.")
    return argp.parse_args()

args = parse_args()
df_list = []
for _dir in sorted(os.listdir(args.validation_dir)):
    for pkl in os.listdir(args.validation_dir + '/' + _dir):
        try:
            with open(args.validation_dir + '/' + _dir + '/' + pkl, 'rb') as f:
                df_list.append(pickle.load(f))
        except Exception as e:
            print(e)
            print('Could not open ' + args.validation_dir + '/' + _dir + '/' + pkl)
df = pd.concat(df_list)
with open(args.outpath, 'wb') as f:
    pickle.dump(df, f, protocol=5)
