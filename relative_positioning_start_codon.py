import numpy as np
import os,sys,argparse
import pandas as pd
import time
import glob
import matplotlib as plt

def alignSamGenic(SamFile, utr3Folder, utr5Folder, cdsFolder, output, normalized_sora):
    df100 = pd.read_csv(SamFile, delim_whitespace=True, names=['Qname', 'flag', 'transcript_id', 'start', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'lol', 'a7a'], usecols=['Qname', 'transcript_id', 'start', 'SEQ'])
    df100['length'] = df100['SEQ'].str.len()
    df100['end'] = df100['start'] + df100['length']
    dfs=[]
    all_files = glob.glob(utr3Folder + "*")
    for filename in all_files:
        df = pd.read_csv(filename, delim_whitespace=True)
        dfs.append(df)
        frame = pd.concat(dfs, axis=0)
    frame1 = frame.sort_values(by='chr')
    frame = frame1.reset_index()
    utr3 = frame[['transcript_id', 'total_length']]
    utr3['transcript_id'] = utr3['transcript_id'].str.replace('_3utr', '')
    utr3.columns = ['transcript_id', 'utr3_length']
    dfs=[]
    all_files = glob.glob(utr5Folder + "*")
    for filename in all_files:
        df = pd.read_csv(filename, delim_whitespace=True)
        dfs.append(df)
        frame = pd.concat(dfs, axis=0)
    frame1 = frame.sort_values(by='chr')
    frame = frame1.reset_index()
    utr5 = frame[['transcript_id', 'total_length']]
    utr5['transcript_id'] = utr5['transcript_id'].str.replace('_5utr', '')
    utr5.columns = ['transcript_id', 'utr5_length']
    dfs=[]
    all_files = glob.glob(cdsFolder + "*")
    for filename in all_files:
        df = pd.read_csv(filename, delim_whitespace=True)
        dfs.append(df)
        frame = pd.concat(dfs, axis=0)
    frame1 = frame.sort_values(by='chr')
    frame = frame1.reset_index()
    cds = frame[['transcript_id', 'total_length']]
    cds['transcript_id'] = cds['transcript_id'].str.replace('_cds', '')
    cds.columns = ['transcript_id', 'cds_length']
    df1 = df100.merge(utr5, on=['transcript_id'], how='left')
    df2 = df1.merge(cds, on=['transcript_id'], how='left')
    df3 = df2.merge(utr3, on=['transcript_id'], how='left')
    df4 = df3[df3['cds_length'].notna()]
    for i in range(1,101):
        df4[str(i) + 'b'] = df4['utr5_length'] + i
        df4['-' + (str(i)) + 'b'] = df4['utr5_length'] - i
        df4[-i] = np.where((df4['-' + str(i) + 'b'] > df4['start']) & (df4['-' + str(i) + 'b'] < df4['end']), 1, 0)
        df4[i] = np.where((df4[str(i) + 'b'] > df4['start']) & (df4[str(i) + 'b'] < df4['end']), 1, 0)

    names = list(range(-100, 100 + 1))
    names += ['transcript_id']

    lol = pd.DataFrame(df4, columns=names)
    lol1 = lol.set_index('transcript_id')
    lol2 = lol1.T
    lol2['sum'] = lol2.sum(axis=1)
    lol2['normalized'] = lol2['sum']/lol2['sum'].sum()
    lol2.to_csv(output, index=False, sep='\t')
    ax = lol2['sum'].plot(style='.-', figsize=(20,10), color='red')
    plt.pyplot.xticks(fontsize=22)
    plt.pyplot.yticks(fontsize=22)
    plt.pyplot.xlabel('\nbins', fontsize=24)
    plt.pyplot.ylabel('reads', fontsize = 24)
    fig = ax.get_figure()
    fig.savefig(normalized_sora, dpi=600)

def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-SamFile', required=True)
    parser.add_argument('-utr3Folder', required=True)
    parser.add_argument('-utr5Folder', required=True)
    parser.add_argument('-cdsFolder', required=True)
    parser.add_argument('-output', required=True)
    parser.add_argument('-normalized_sora', required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = getArgs()
    start = time.time()
    filename = alignSamGenic(args.SamFile, args.utr3Folder, args.utr5Folder, args.cdsFolder, args.output, args.normalized_sora)
    end = time.time()
    print ('time elapsed:' + str(end - start))
