#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import json
import re

def load_list(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]

def load_missingness(path):
    df = pd.read_csv(path, sep='\t')
    if 'Sample' not in df.columns or 'Percent_Missingness' not in df.columns:
        raise ValueError(f"{path} must contain 'Sample' and 'Percent_Missingness'")
    df['Sample'] = df['Sample'].astype(str).str.strip()
    df['Missing'] = df['Percent_Missingness']
    return df.drop(columns=['Percent_Missingness'])

def parse_html_header(path):
    meta = {}
    with open(path) as f:
        for line in f:
            text = line.strip(" \n\t-<!>")
            m = re.match(r'^([A-Za-z0-9_]+):\s*"(.*)"$', text)
            if m:
                meta[m.group(1)] = m.group(2)
    return meta

def write_mqc_json(df, metadata, output):
    data = {
        row['Sample']: {k: v for k, v in row.items() if k != 'Sample'}
        for _, row in df.iterrows()
    }
    pconfig = {
        'id':        metadata.get('id', metadata.get('section_name')),
        'ylab':      'Value',
        'xlab':      'Metric',
        'xDecimals': False,
        'tt_label':  'Metric',
        'min':       0.0,
        'max':       100.0,
        'scale':     'YlGnBu'
    }
    out = {'data': data, 'pconfig': pconfig}
    for k, v in metadata.items():
        if k not in ('id'):
            out[k] = v
    with open(output, 'w') as f:
        json.dump(out, f, indent=2)

def main(args):
    inds_pre  = load_list(args.inds_pre)
    inds_post = load_list(args.inds_post)
    miss_pre  = load_missingness(args.miss_pre)[['Sample', 'Missing']].rename(columns={'Missing':'Missing_Pre'})
    miss_post = load_missingness(args.miss_post)[['Sample', 'Missing']].rename(columns={'Missing':'Missing_Post'})
    df = pd.DataFrame({'Sample': inds_pre})
    df = df.merge(miss_pre,  on='Sample', how='left') \
           .merge(miss_post, on='Sample', how='left')
    df.loc[~df['Sample'].isin(inds_post), 'Missing_Post'] = np.nan
    if args.header:
        meta = parse_html_header(args.header)
        write_mqc_json(df, meta, args.output)
    else:
        df.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--inds-pre',  required=True)
    p.add_argument('--inds-post', required=True)
    p.add_argument('--miss-pre',  required=True)
    p.add_argument('--miss-post', required=True)
    p.add_argument('--output',    required=True)
    p.add_argument('--header',    required=False)
    main(p.parse_args())
