import datashader as ds
import datashader.transfer_functions as tf
from holoviews.operation.datashader import rasterize
import logging
import hvplot.pandas
import param
import panel as pn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import holoviews as hv
from holoviews import opts
from holoviews.operation.datashader import datashade, dynspread
from motif_plots import motif_plot_depletion
import colorcet as cc
long_key = list(set(cc.glasbey_warm + cc.glasbey_cool + cc.glasbey_dark))

ROW_HEIGHT = 20
COL_WIDTHS = {'count':50, 'chosen_allele':200}
stylesheet = """
.tabulator-cell {
    font-size: 12px;
}
"""

def get_dashboard(df, colorby_cols=None, agg_func=None, use_gpu=True):
    if colorby_cols is None: colorby_cols = [ 'split','mhctype', 'chosen_allele_type', 'EL_pred', 'peptide_length', 'dataset', 'chosen_allele', 't_attn3_max', 't_attn3_argmax']
    if agg_func is None: agg_func =     [ds.count_cat,ds.count_cat, ds.count_cat, ds.mean, ds.count_cat, ds.count_cat, ds.count_cat, ds.mean, ds.mean]
    
   
    
    class Dashboard(param.Parameterized):
        trigger_return = param.Integer(default=0)
        trigger_update_scatter = param.Integer(default=0)
        cluster_label_option = param.ObjectSelector(default=0, objects=[0,1,2,3,4,5,6,7,8,9,10])

        def __init__(self, df, agg_func, colorby_cols, use_gpu=False):
            super().__init__()
            
            if 'index' not in colorby_cols: colorby_cols = colorby_cols + ['index']
            self.agg_func, self.colorby_cols, self.use_gpu = agg_func, colorby_cols, use_gpu
            self.setup_logging()
            self.debug = 10
            
            df['cluster_label'] = 0
            df['chosen_allele_type'] = df['chosen_allele'].str[:2]
            for col in colorby_cols:
                if col not in df.columns:
                    df[col] = 0

            for ix, val in enumerate(self.agg_func):
                if val is ds.count_cat:
                    df[colorby_cols[ix]] = df[colorby_cols[ix]].astype('category')
                else:
                    df[colorby_cols[ix]] = df[colorby_cols[ix]].fillna(0)
            df['index'] = df.index.values
            datacols = ['peptide', 'peptide_core', 'nFlank', 'cFlank', 'emb_0', 'emb_1']
            self.df = df[colorby_cols+datacols].copy()
            self.color_option = pn.widgets.Select(name='Color Option',value='chosen_allele', options=colorby_cols)
            self.color_option.link(self, callbacks={'value': self.query_update_func})
            self.motif_option = pn.widgets.Select(name='Motif Option',value='MHC2_BC', options=['MHC2_BC','MHC1_8mer', 'MHC1_9mer', 'MHC1_10mer', 'MHC1_11mer', 'MHC1_12mer'])
            self.motif_option.link(self, callbacks={'value': self.motif_update_func})
            
            self.density_option = pn.widgets.Select(name='Density Option',value='sparse', options=['sparse', 'dense'])
            self.density_option.link(self, callbacks={'value': self.query_update_func})
            self.p_scatter = None
            plt.ioff()
            self.fig, self.axes = plt.subplots(1, 1, figsize=(2,1.5), squeeze=True, tight_layout=True, )
            self.fig_2, self.axes_2 = plt.subplots(1, 1, figsize=(2,1.5), squeeze=True, tight_layout=True)
            self.fig_3, self.axes_3 = plt.subplots(1, 1, figsize=(2,1.5), squeeze=True, tight_layout=True)

            self.fig_4, self.axes_4 = plt.subplots(1, 1, figsize=(2,1.5), squeeze=True, tight_layout=True)
            self.fig_5, self.axes_5 = plt.subplots(1, 1, figsize=(2,1.5), squeeze=True, tight_layout=True)
            plt.ion()
            self.link = hv.link_selections.instance(unselected_alpha=1)
            self.link.param.watch(self.trigger_updates, 'selection_expr')       
            self.gspec = pn.GridSpec(width=1600,height=1000)
            self.user_query = pn.widgets.TextInput(name='Query Value')
            self.user_query.link(self, callbacks={'value': self.user_query_func})                
            self.update_df_query_indexs()
            self.update_df_selected()
            self.update_scatter()
            self.update_plots()


        def setup_logging(self):
            LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
            file_handler = logging.FileHandler(filename='test.log', mode='w')
            file_handler.setFormatter(logging.Formatter(LOG_FORMAT))
            self.logger = logging.getLogger(__name__)
            self.logger.addHandler(file_handler)
            self.logger.setLevel(logging.INFO)

            self.debug_info = pn.widgets.Debugger(name='Debugger info level', level=logging.INFO, sizing_mode='stretch_both',
                                    logger_names=[__name__],#comment this line out to get all panel errors
                                    )
            self.view_debug = pn.Column(self.debug_info, sizing_mode='stretch_both')
            self.logger.info('Logger Created')

            
        def motif_update_func(self, event=None, event2=None):
            self.logger.info('motif_update_func')
            self.trigger_updates(None)
            
        def user_query_func(self, event=None, event2=None):
            self.logger.info(self.user_query.value)
            self.query_update_func()
            
        def query_update_func(self, event=None, event2=None):
            self.logger.info('query_update_func')
            self.link.selection_expr = None
            self.update_df_query_indexs()
            self.update_scatter()

        def update_df_query_indexs(self):
            try:
                if self.user_query.value:
                    self.logger.info(f'Filtering by {self.user_query.value}')
                    self.df_query_indexs = self.df.query(self.user_query.value).index.values
                else:
                    self.df_query_indexs = None
            except:                
                self.logger.exception('Dataframe query or update failed')
        
        def update_df_selected(self):    
            if self.df_query_indexs is not None:
                self.df_selected_indexs = self.link.filter(self.df.loc[self.df_query_indexs,:]).index.values
            else:
                self.df_selected_indexs = self.link.filter(self.df).index.values
            self.df_selected = self.df.loc[self.df_selected_indexs,:]
            result = []
            for col in ['mhctype', 'chosen_allele', 'dataset']:
                summary = self.df_selected[col].value_counts().head(5).reset_index()
                summary['group'] = col
                result.append(summary.rename(columns={col:'name'}))
            self.df_selected_summary = pd.concat(result)
            
        def trigger_updates(self, event):
            self.logger.info('trigger_updates')
            self.update_df_selected()
            self.update_plots()
            self.param.trigger('trigger_return')

        def update_plots(self):
            try:
                self.axes.clear() 
                
                if self.motif_option.value == 'MHC2_BC':
                    df_local = self.df_selected.query('mhctype=="mhc2"')
                    if df_local.shape[0] > 4: motif_plot_depletion(self.axes, df_local['peptide_core'].str.pad(9, side='right', fillchar='*').values, special_chars=False)
                elif self.motif_option.value == 'MHC1_9mer':
                    df_local = self.df_selected.query('mhctype=="mhc1" and peptide_length==9')
                    if df_local.shape[0] > 4: motif_plot_depletion(self.axes, df_local['peptide'].values)
                elif self.motif_option.value == 'MHC1_10mer':
                    df_local = self.df_selected.query('mhctype=="mhc1" and peptide_length==10')
                    if df_local.shape[0] > 4: motif_plot_depletion(self.axes, df_local['peptide'].values)
                elif self.motif_option.value == 'MHC1_11mer':
                    df_local = self.df_selected.query('mhctype=="mhc1" and peptide_length==11')
                    if df_local.shape[0] > 4: motif_plot_depletion(self.axes, df_local['peptide'].values)
                elif self.motif_option.value == 'MHC1_12mer':
                    df_local = self.df_selected.query('mhctype=="mhc1" and peptide_length==12')
                    if df_local.shape[0] > 4: motif_plot_depletion(self.axes, df_local['peptide'].values)
                elif self.motif_option.value == 'MHC1_8mer':
                    df_local = self.df_selected.query('mhctype=="mhc1" and peptide_length==8')
                    if df_local.shape[0] > 4: motif_plot_depletion(self.axes, df_local['peptide'].values)


                if self.motif_option.value == 'MHC2_BC':
                    self.axes_2.clear()
                    self.axes_3.clear()
                    self.axes_4.clear()
                    self.axes_5.clear()
                    df_local = self.df_selected.query('mhctype=="mhc2"')
                    if df_local.shape[0] > 4:
                        motif_plot_depletion(self.axes_2, df_local['nFlank'].str.pad(5, side='left', fillchar='*').values)
                        motif_plot_depletion(self.axes_3, df_local['cFlank'].str.pad(5, side='right', fillchar='*').values)
                        
                else:
                    self.axes_2.clear()
                    self.axes_3.clear()
                    df_local = self.df_selected.query('mhctype=="mhc1"')
                    if df_local.shape[0] > 4:
                        motif_plot_depletion(self.axes_2, df_local['nFlank'].str.pad(5, side='left', fillchar='*').values)
                        motif_plot_depletion(self.axes_3, df_local['cFlank'].str.pad(5, side='right', fillchar='*').values)

            except:
                self.logger.exception('Plot generation failed')


        def update_scatter(self):
            self.logger.info('update_scatter')
            try:
                #self.ds = hv.Dataset(self.df_scatter)
                #self.link = hv.link_selections.instance(selected_color=None, unselected_color=None, unselected_alpha=1)
                selected_color = self.color_option.value
                selected_color_idx = np.argwhere(np.array(self.color_option.options) == selected_color)[0][0]
                aggregator = self.agg_func[selected_color_idx](selected_color)
                if self.df[selected_color].dtype.name == 'category':
                    ck = long_key
                    cm = long_key
                else:
                    ck = "viridis"
                    cm = "viridis"
                if self.df_query_indexs is not None:
                    points = self.df.loc[self.df_query_indexs,:].hvplot(x='emb_0', y='emb_1', kind='points', hover_cols=self.colorby_cols)
                else:
                    points = self.df.hvplot(x='emb_0', y='emb_1', kind='points', hover_cols=self.colorby_cols)

                shaded = datashade(points, expand=False, aggregator=aggregator, color_key=ck, min_alpha=0, x_sampling=0.25, y_sampling=0.25, vdim_prefix='',
                                    cnorm="eq_hist", precompute=True, dynamic=True).opts(responsive=True)
                                        
                self.p_scatter = self.link(shaded)     
                self.param.trigger('trigger_return')
            except:
                self.logger.exception('Scatter plot failed')


        @param.depends('cluster_label_option')
        def save_label(self):
            self.df.loc[self.link.filter(self.df).index, 'cluster_label'] = self.cluster_label_option
            self.temp = 'updated'
            return self.param.cluster_label_option
        
        

        def view_scatter(self): return self.p_scatter
        
        def view_binding_core_confidence(self): return self.df_selected['t_attn3_max'].hvplot.hist(responsive=True, shared_axes=False)
        
        def view_binding_core_location(self): return self.df_selected['t_attn3_argmax'].hvplot.hist(responsive=True, shared_axes=False)
        
        def view_el_pred(self): return self.df_selected['EL_pred'].hvplot.hist(responsive=True,  shared_axes=False)
        
        def view_p_bar(self, event=None): return self.df_selected['peptide_length'].hvplot.hist(responsive=True, shared_axes=False, bins=np.linspace(self.df_selected['peptide_length'].min()-0.5, self.df_selected['peptide_length'].max()+0.5, int(self.df_selected['peptide_length'].max()-self.df_selected['peptide_length'].min()+2)))
    
        def view_split_df(self): return pn.widgets.Tabulator(self.df_selected['split'].value_counts().reset_index().head(),show_index=False, configuration={'rowHeight':ROW_HEIGHT},stylesheets=[stylesheet], text_align='left', sizing_mode='stretch_width')

        def view_mhctype_df(self): return pn.widgets.Tabulator(self.df_selected['mhctype'].value_counts().reset_index().head(), show_index=False, configuration={'rowHeight':ROW_HEIGHT},stylesheets=[stylesheet], text_align='left', sizing_mode='stretch_width')
            
        def view_chosen_allele_df(self): return pn.widgets.Tabulator(self.df_selected['chosen_allele'].value_counts().head(n=10).reset_index(), show_index=False, configuration={'rowHeight':ROW_HEIGHT},stylesheets=[stylesheet], text_align='left', sizing_mode='stretch_width')

        def view_dataset_df(self): return pn.widgets.Tabulator(self.df_selected['dataset'].value_counts().head().reset_index(), show_index=False, configuration={'rowHeight':ROW_HEIGHT},stylesheets=[stylesheet], text_align='left', sizing_mode='stretch_width')
    
        def view_motif(self): return pn.Tabs(('Binding Core', self.fig), ('nFlank', self.fig_2), ('cFlank', self.fig_3))

        def view_summary_df(self): return pn.widgets.Tabulator(self.df_selected_summary, show_index=False, configuration={'rowHeight':ROW_HEIGHT},hidden_columns=['group'], stylesheets=[stylesheet], text_align='left', sizing_mode='stretch_width', groupby=['group'])

        def view_bars(self): return pn.Row((self.view_binding_core_location() + self.view_binding_core_confidence() + self.view_p_bar() + self.view_el_pred()).cols(4).opts(toolbar='right'))

        def view_save_selected(self):
            button = pn.widgets.Button(name='Save Selected points to ./df_selected.parquet', button_type='primary', disabled=False)
            button.on_click(lambda event: self.df_selected.to_parquet('df_selected.parquet')) 
            return button

        def show(self):
            self.gspec[0,0] = pn.Column(self.color_option)
            self.gspec[1,0] = self.user_query
            self.gspec[0,2] = self.view_save_selected
            self.gspec[2:13,0] = self.view_summary_df
            self.gspec[13:20,0] = pn.Column(self.motif_option, self.view_motif)
            self.gspec[1:16,   1:4] = self.view_scatter
            self.gspec[16:20,   1:4] =  self.view_bars 
            return pn.panel(self.gspec)
    
    dash = Dashboard(df, agg_func=agg_func, colorby_cols=colorby_cols, use_gpu=use_gpu)
    return dash

