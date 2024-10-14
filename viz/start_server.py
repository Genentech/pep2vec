import panel as pn
import holoviews as hv
import dashboard_v3 as dashboard
import pandas as pd
import numpy as np
import datashader as ds
import warnings
warnings.filterwarnings("ignore")
pn.extension('terminal')
pn.extension('tabulator')
hv.extension('bokeh')
print('Loading Data')
df_temp = pd.read_parquet('./df_for_vizualization.parquet').reset_index(drop=True)
df_temp['peptide_length'] = df_temp['peptide'].apply(len)
if 'cFlank' in df_temp.columns: df_temp = df_temp.drop(columns=['cFlank'])
if 'nFlank' in df_temp.columns: df_temp = df_temp.drop(columns=['nFlank'])
df_temp = df_temp.rename(columns={'n_flank':'nFlank', 'c_flank':'cFlank'})
if 't_attn3_max' not in df_temp.columns: df_temp['t_attn3_max'] = 0
if 'chosen_allele' not in df_temp.columns: df_temp['chosen_allele'] = df['allele']
df_temp[['EL_pred', 'emb_0', 'emb_1', 'peptide_length']] = df_temp[['EL_pred', 'emb_0', 'emb_1', 'peptide_length']].astype(float)
df_temp[['peptide', 'peptide_core', 'nFlank', 'cFlank']] = df_temp[['peptide', 'peptide_core', 'nFlank', 'cFlank']].astype(str)
colorby_cols = [ 'split','mhctype', 
            'chosen_allele_type', 'EL_pred', 'peptide_length', 
            'dataset', 'chosen_allele', 't_attn3_max', 't_attn3_argmax']
agg_func =     [ds.count_cat,ds.count_cat, 
            ds.count_cat, ds.mean, ds.mean, 
            ds.count_cat, ds.count_cat, ds.mean, ds.mean]
dash = dashboard.get_dashboard(df=df_temp, agg_func=agg_func, colorby_cols=colorby_cols, use_gpu=True)




pn.template.VanillaTemplate(
    site="pep2vec",
    title="vizualization",
    sidebar=None,
    main=dash.show(),#
).servable()

