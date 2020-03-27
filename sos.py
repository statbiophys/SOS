'''
    Copyright (C) 2020 Isacchini Giulio
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.'''
import numpy as np
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import flask
import os
import olga.load_model as olga_load_model
import olga.generation_probability as pgen
import pandas as pd
import dash_daq as daq
from utils import *
from minimal_sonia import MinimalSonia
import urllib.parse as urllib
from textwrap import dedent
data_tables=[pd.read_csv(local_directory+'data_'+i+'.csv') for i in options_of]
app = dash.Dash(__name__,meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}])
import plotly.graph_objs as go
import pathlib
import time
import base64

test_base64 = base64.b64encode(open(local_directory+'plots.png', 'rb').read()).decode('ascii')
pobs = base64.b64encode(open(local_directory+'pobs.png', 'rb').read()).decode('ascii')

BASE_PATH = pathlib.Path(__file__).parent.resolve()

############################
# HERE I DEFINE THE LAYOUT #
############################

app.layout = html.Div([
                html.Div(
                html.Div(
                id='olga-control-tabs',
                className='control-tabs',
                children=[
                    dcc.Tabs(id='olga-tabs', value='what-is', children=[
                        dcc.Tab(
                            label='About',
                            value='what-is',
                            children=html.Div(children=[
                                html.P(id='ciao',children='801',style={'width': '50%', 'display': 'none'}),
                                html.P(style={'padding-top':10}),
                                html.H4(className='what-is', children='Simple Olga Sonia (SOS)',style={'text-align':'center'}),
                                html.P('SOS evaluates the probability of generation (Pgen) and probability in the periphery (Ppost) of  '
                                    'specific T and B cell receptor sequences of humans and mice.'),
                                html.P('SOS also generates synthetic generated and peripheral repertoires of naive cells.'),
                                html.P('SOS is built on top of the  ',style={'display':'inline-block'}),
                                html.A('IGOR',href='https://github.com/qmarcou/IGoR', target='_blank',style={'display':'inline-block','padding-left':'5px'}),
                                html.P('[2], ',style={'display':'inline-block'}),
                                html.A('OLGA ',href='https://github.com/zsethna/OLGA', target='_blank',style={'display':'inline-block','padding-left':'5px'}),
                                html.P('[3] and',style={'display':'inline-block'}),
                                html.A(' SONIA ',href='https://github.com/statbiophys/SONIA/', target='_blank',style={'display':'inline-block','padding-left':'5px'}),
                                html.P('[4,5] software.',style={'display':'inline-block'}),
                                html.P('Pgen is conditioned on productivity of the sequence. It only takes into account VDJ '
                                        'recombination in the absence of selection. Ppost takes selection into account, and should '
                                        'better reflect frequency of occurrence in the periphery. However note that it is insensitive to '
                                        'HLA type. The models were learned on a possibly biased panel of individuals and may '
                                        'depend on the sequencing technologies used. For more accurate dataset-specific comparisons '
                                        'it is advised to infer new models using IGOR and SONIA.'),
                                html.P('Evaluate command: '),
                                html.P( 'sonia-evaluate --humanTRB CASSTGNYGAFF --v_mask TRBV9 --j_mask TRBJ-1',
                                    style={'fontFamily':'Courier New, monospace'}),
                                html.P('Generate command: '),
                                html.P('sonia-generate --humanTRB -n 1000 --pgen (or --ppost)',
                                    style={'fontFamily':'Courier New, monospace'}),
                                html.P('To cite this tool: '),
                                html.P('[1] Isacchini et al (2020) SOS: A web server for SONIA and OLGA, in preparation.'),
                                html.P('References:'),
                                html.A('[2] Marcou Q, Mora T, Walczak AM (2018) High-throughput immune repertoire analysis with IGoR. Nat Commun 9, 561 (2018) ', target='_blank',href='https://www.nature.com/articles/s41467-018-02832-w'),
                                html.P(),
                                html.A('[3] Sethna Z, Elhanati Y, Callan CG, Walczak AM, Mora T(2019) OLGA: fast computation of'
                                        ' generation probabilities of B- and T-cell receptor amino acid sequences and motifs. '
                                        'Bioinformatics, Volume 35, Issue 17, 1 September 2019, Pages 2974–2981', target='_blank',
                                        href='https://academic.oup.com/bioinformatics/article/35/17/2974/5292315',),
                                html.P(),
                                html.A('[4] Sethna Z, Isacchini G, Dupic T, Mora T, Walczak AM, Elhanati Y, '
                                    'Population variability in the generation and thymic selection of T-cell repertoires, (2020) bioRxiv, https://doi.org/10.1101/2020.01.08.899682 ', target='_blank',href= ' https://www.biorxiv.org/content/10.1101/2020.01.08.899682v1'),
                                html.P(),
                                html.A('[5]  Isacchini G,Sethna Z, Elhanati Y ,Nourmohammad A, Mora T, Walczak AM, '
                                    'On generative models of T-cell receptor sequences,(2019) bioRxiv, https://doi.org/10.1101/857722 ', target='_blank',href=' https://www.biorxiv.org/content/10.1101/857722v2')

                            ],style={'maxHeight': '60vh', 'overflow': 'scroll'}), selected_style=tab_selected_style
                        ),

                        dcc.Tab(
                            label='Evaluate',
                            value='what-is-not',
                            children=html.Div(className='control-tab',children=[
                                html.P(style={'padding-top':10}),
                                html.H4(className='what-is-not', children='Evaluate the probability',style={'text-align':'center'}),

                                html.P(children=''' Select the model'''),
                                dcc.Dropdown(
                                    id='name-dropdown',
                                    options=[{'label':name, 'value':i} for i,name in enumerate(names)],
                                    value = 0,clearable=False),
                                html.P(children='''Input CDR3,V and J genes ''',style={'padding-top':10}),
                                html.Div([dcc.Dropdown(id='cdr3-type',options=[{'label': 'aminoacid', 'value': 1},
                                {'label': 'nucleotide', 'value': 0}],value=1,clearable=False)],style={'width': '50%','display': 'inline-block','vertical-align': 'middle'}),
                                dcc.Input(id='cdr3',value='CASSDAQGRNRGTEAFF', type='text',style={'width': '50%','backgroundColor': '#ffffff','display': 'inline-block','vertical-align': 'middle''inline-block'}),
                                html.Div([dcc.Dropdown(id='v_option',clearable=True,value='TRBV9')],style={'width': '50%', 'display': 'inline-block'}),
                                html.Div([dcc.Dropdown(id='j_option',clearable=True,value='TRBJ1-1')],style={'width': '50%', 'display': 'inline-block'}),
                                html.P(children=' '),
                                html.Button(id='submit-button', children='PGEN',style={'width': '50%'},title='olga-compute_pgen --humanTRB CASSVGGNTEAFF --v_mask TRBV9 --j_mask TRBJ-1'),
                                html.Button(id='fake-button', children='Compute',style={'display': 'none'}),
                                html.Div(id='model-type',children='human_T_beta',style={'display':'none'}),
                                html.P(style={'padding-top':10}),
                                dcc.Loading(id="loading-2", children=[
                                    html.Div(id='pgen-nucleotide',children=1e-80,style={'display':'none'}),
                                    html.Div(id='pgen-value',children=1e-80,style={'display':'none'}),
                                    html.P(id='print-pgen-nucleotide',children='',style={'width': '100%'}),
                                    html.P(id='print-pgen',children='',style={'width': '100%'}),
                                    ]),
                                html.Button(id='submit-ppost', children='PPOST',style={'width': '50%'},title='sonia-evaluate --humanTRB CASSVGGNTEAFF --v_mask TRBV9 --j_mask TRBJ-1'),
                                html.P(style={'padding-top':10}),
                                dcc.Loading(id="loading-3", children=[
                                    html.Div(id='ppost-value',children=1e-80,style={'display':'none'}),
                                    html.Div(id='pgen-for-ppost',children=1e-80,style={'display':'none'}),
                                    html.P(id='print-ppost',children='',style={'width': '100%'})
                                     ]),

                                html.Div([
                                html.P(children='''How many unique productive nucleotide sequences do you have?'''),
                                dcc.Input(id='n-productive',value=100000, type='number',style={'width': '100%','backgroundColor': '#ffffff',})
                                 ],style={'width': '100%','display': 'inline-block','vertical-align':'middle','padding-bottom': 10}),
                                html.Br(),
                                html.Div([
                                html.P(children='Probability to find the aminoacid sequence by pure chance (p-value):', style={'width': '100%','display':'inline-block','vertical-align':'middle','padding-left':'1'}),
                                html.Img(src='data:image/png;base64,{}'.format(pobs),style={'width': '50%', 'display': 'inline-block'}),

                                html.P(id='probability-to-observe',children='', style={'width': '30%','display':'inline-block','vertical-align':'middle'})],

#                                daq.GraduatedBar(id='my-graduated-bar',value=0,min=0,max=100,step=5,showCurrentValue=True,style={})],

                                style={'width': '100%'})],style={'width': '100%','display': 'inline-block','vertical-align': 'middle'}),style=tab_style, selected_style=tab_selected_style),
                        dcc.Tab(
                            label='Generate',
                            value='generate-seq',
                            children=html.Div([
                                html.P(style={'padding-top':10}),
                                html.H4(className='what-is', children='Generate synthetic repertoires',style={'text-align':'center'}),
                                html.P(children=''' Select the model''',style={'vertical-align': 'middle'}),
                                dcc.Dropdown(
                                    id='name-dropdown-copy',
                                    options=[{'label':name, 'value':i} for i,name in enumerate(names)],
                                    value = 0,clearable=False,
                                    style={'width': '100%','vertical-align': 'middle'}),
                                html.P(style={'padding-top':10}),
                                html.P(children='Generate synthetic sequences from the model (max 1000)', style={'vertical-align': 'middle'}),
                                html.Div([dcc.Dropdown(id='sampling-origin',options=[{'label': 'pre-selection', 'value': 0},
                                {'label': 'post-selection', 'value': 1}],value=0,clearable=False)],style={'width': '33%','display': 'inline-block','vertical-align': 'middle'}),
                                dcc.Input(id='n-generation', placeholder="n sequences", type='number',style={'width': '33%','backgroundColor': '#ffffff','display': 'inline-block','vertical-align': 'middle'}),
                                dcc.Input(id="seed", type="number", placeholder="seed",style={'width': '33%','backgroundColor': '#ffffff','display': 'inline-block','vertical-align': 'middle'}),
                                html.P(style={'padding-top':10}),
                                html.Button(id='generate-button', children='GENERATE',style={'width': '50%'},title='sonia-generate --humanTRB -n 1000 --pgen (or --ppost)'),
                                html.P(style={'padding-top':10}),
                                dcc.Loading(id="loading-1", children=
                                    [ html.P(id='file', children='',style={'display': 'none'}),
                                    html.Div(html.A('download', id='my-link'),style={'width': '100%','height': '100%','vertical-align': 'middle','textAlign':'center'})],
                                    type="default",style={'width': '100%','height': '100%','textAlign':'center'})
                            ],style={'width': '100%','display': 'inline-block','vertical-align':'middle','position': 'relative'})
                        ,style=tab_style, selected_style=tab_selected_style),

                        ]
                        ,style=tabs_styles,colors={'primary':'gold',"background": "cornsilk"})],style={'width': '80%','padding-left':'10%', 'padding-right':'0%',
                                    'height': '80%','padding-bottom':'10%', 'padding-top':'10%'
                                            }),id='mobile',style={'width': '50%','display': 'inline-block','vertical-align':'top','height': '90%'}),
        # PLOTS
        html.Div([
        #variables (NOT SHOWN BUT STORED IN THE LAYOUT)
        html.Div(html.Img(src='data:image/png;base64,{}'.format(test_base64),style={'width': '90%','height': '100%','padding-left':'20%','padding-top':'20%'}),
            style={'width': '50%','height': '50%','display': 'inline-block','vertical-align':'top'}) ,
        html.Div(dcc.Graph(
        id='Q-graph',
        config={ 'staticPlot': True},
        style={'height': '35vh','width': '120%',},

        figure={'data': [
                {'x': data_tables[0].binning_Q.values, 'y': data_tables[0].Q_gen.values, 'type': 'line', 'name': 'pre-selection'},
                {'x': data_tables[0].binning_Q.values, 'y': data_tables[0].Q_data.values, 'type': 'line', 'name': 'post-selection'}],
        "layout": go.Layout(
                    paper_bgcolor = 'rgba(0,0,0,0)',
                    plot_bgcolor = 'rgba(0,0,0,0)',
                    xaxis = {
                    "side": "bottom",
                    "type": "linear",
                    "range": [data_tables[0].binning_Q.values.min(),data_tables[0].binning_Q.values.max()],
                    "title": " log10 Q",
                    },
                    showlegend=False,
                    yaxis= {
                        "side": "left",
                        "type": "linear",
                        "range": [0,data_tables[0].Q_data.values.max()+0.8]},
                    margin = {'t': 0,'b':100},
                    legend={'x':0,'y':1})}),style={'width': '50%','height': '50%','display': 'inline-block','vertical-align':'top'}) ,
html.Div(
        dcc.Graph(
        id='ppost-graph',
        config={ 'staticPlot': True},
        style={'height': '35vh','width': '120%',},
        figure={'data': [
                {'x': data_tables[0].binning_ppost.values, 'y': data_tables[0].ppost_gen.values, 'type': 'line', 'name': 'pre-selection'},
                {'x': data_tables[0].binning_ppost.values, 'y': data_tables[0].ppost_data.values, 'type': 'line', 'name': 'post-selection'}],
        "layout": go.Layout(
                    paper_bgcolor = 'rgba(0,0,0,0)',
                    plot_bgcolor = 'rgba(0,0,0,0)',
                    xaxis = {
                    "side": "bottom",
                    "type": "linear",
                    "range": [data_tables[0].binning_ppost.values.min(),data_tables[0].binning_ppost.values.max()],
                    "title": " log10 Ppost",
                    },
                    showlegend=False,
                    yaxis= {
                        "side": "left",
                        "type": "linear",
                        "range": [0,data_tables[0].ppost_data.values.max()+0.1]},
                    margin = {'t': 0,'b':100},
                    legend={'x':0,'y':1})}),style={'width': '50%','height': '50%','display': 'inline-block','vertical-align':'top'}) ,
        html.Div(dcc.Graph(
        id='pgen-graph',
        config={ 'staticPlot': True},
        style={'height': '35vh','width': '120%',},
        figure={'data': [
                {'x': data_tables[0].binning_pgen.values, 'y': data_tables[0].pgen_gen.values, 'type': 'line', 'name': 'pre-selection'},
                {'x': data_tables[0].binning_pgen.values, 'y': data_tables[0].pgen_data.values, 'type': 'line', 'name': 'post-selection'}],
        "layout": go.Layout(
                    paper_bgcolor = 'rgba(0,0,0,0)',
                    plot_bgcolor = 'rgba(0,0,0,0)',
                    xaxis = {
                    "side": "bottom",
                    "type": "linear",
                    "range": [data_tables[0].binning_pgen.values.min(),data_tables[0].binning_pgen.values.max()],
                    "title": " log10 Pgen",
                    },
                    yaxis= {
                        "side": "left",
                        "type": "linear",
                        "range": [0,data_tables[0].pgen_data.values.max()+0.1]},
                    margin = {'t': 0,'b':100},
                    legend={'x':0,'y':1})}),style={'width': '50%','height': '50%','display': 'inline-block','vertical-align':'top'}) ,


        ],id='mobile2',
        style={'width': '45%','height': '80%','display': 'inline-block','vertical-align': 'top',
                'padding-left':'0%', 'padding-right':'5%',
                'padding-top':'5%', 'padding-bottom':'5%',}),

])
default_seqs=[['CASSDAQGRNRGTEAFF','TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT','TRBV9','TRBJ1-1'],
            ['CAYRLIQGAQKLVF','TGTGCTTATAGATTAATTCAGGGAGCCCAGAAGCTGGTATTT','TRAV38-2/DV8','TRAJ54'],
            ['CARHYGSGSEWFDPW','TGTGCGAGACACTATGGTTCGGGGAGTGAATGGTTCGACCCCTGG','IGHV4-34','IGHJ5'],
            ['CMQALQTPQITF','TGCATGCAAGCTCTACAAACTCCTCAGATCACCTTC','IGKV2-28','IGKJ5'],
            ['CQSADSSGQVF','TGTCAATCAGCAGACAGCAGTGGACAGGTCTTC','IGLV3-25','IGLJ1'],
            ['CTCSLTSNERLFF','TGCACCTGCAGTTTGACATCCAACGAAAGATTATTTTTC','TRBV1','TRBJ1-4']]
###############################
# HERE I DEFINE THE FUNCTIONS #
###############################
app.clientside_callback(
    """
    function(s) {
        return screen.width;
    }
    """,
    Output('ciao', 'children'),
    [Input('fake-button', 'n_clicks')]
)

@app.callback(
    Output('cdr3', 'value'),
    [Input('cdr3-type', 'value'),Input('name-dropdown', 'value')])
def update_bar(value,index):
    if value>0: return default_seqs[index][0]
    else: return default_seqs[index][1]

@app.callback(
    Output('mobile', 'style'),
    [Input('ciao', 'children')])
def update_mobile(children):
    if int(children)>800: return {'width': '50%','display': 'inline-block','vertical-align':'top','height': '90%'}
    else:return {'width': '100%','display': 'inline-block','vertical-align':'top','height': '90%'}

@app.callback(
    Output('mobile2', 'style'),
    [Input('ciao', 'children')])
def update_mobile(children):
    if int(children)>800:  return {'width': '45%','height': '80%','display': 'inline-block','vertical-align': 'top',
                'padding-left':'0%', 'padding-right':'5%',
                'padding-top':'5%', 'padding-bottom':'5%',}
    else: return {'width': '90%','display': 'none','vertical-align': 'top',
                'padding-left':'5%', 'padding-right':'5%',
                'padding-top':'5%', 'padding-bottom':'5%','height': '90%'}

#@app.callback(
#    Output('my-graduated-bar', 'value'),
#    [Input('probability-to-observe', 'children')])
#ef update_bar(value):
#    return value

@app.callback(Output('pgen-graph', 'figure'), [Input('pgen-value', 'children'),Input('name-dropdown', 'value'),Input('v_option', 'value'),Input('j_option', 'value')])
def update_graph_pgen(pgen_value,index_,v_gene,j_gene):
    if v_gene is None or j_gene is None:
       return {'data': [],
            "layout":
                go.Layout(
                    paper_bgcolor = 'rgba(0,0,0,0)',
                    plot_bgcolor = 'rgba(0,0,0,0)',
                    xaxis = {
                    "side": "bottom",
                    "type": "linear",
                    "range": [data_tables[index_].binning_pgen.values.min(),data_tables[index_].binning_pgen.values.max()],
                    "title": " log10 Pgen",
                    },

                    yaxis= {
                        "side": "left",
                        "type": "linear",
                        "range": [0,data_tables[index_].pgen_data.values.max()+0.1]},
                    margin = {'t': 0,'b':100},
                    legend={'x':0,'y':1})
        }
    return {'data': [
                {'x': data_tables[index_].binning_pgen.values, 'y': data_tables[index_].pgen_gen.values, 'type': 'line', 'name': 'pre-selection'},
                {'x': data_tables[index_].binning_pgen.values, 'y': data_tables[index_].pgen_data.values, 'type': 'line', 'name': 'post-selection'},

                {'x': [np.log10(pgen_value+1e-80), np.log10(pgen_value+1e-80)], 'y': [-1,2], 'type': 'line', 'name': 'your seq'},
            ],
            "layout":
                go.Layout(
                    paper_bgcolor = 'rgba(0,0,0,0)',
                    plot_bgcolor = 'rgba(0,0,0,0)',
                    xaxis = {
                    "side": "bottom",
                    "type": "linear",
                    "range": [data_tables[index_].binning_pgen.values.min(),data_tables[index_].binning_pgen.values.max()],
                    "title": " log10 Pgen",
                    },

                    yaxis= {
                        "side": "left",
                        "type": "linear",
                        "range": [0,data_tables[index_].pgen_data.values.max()+0.1]},
                    margin = {'t': 0,'b':100},
                    legend={'x':0,'y':1})
        }

@app.callback(Output('ppost-graph', 'figure'), [Input('ppost-value', 'children'),Input('pgen-for-ppost', 'children'),Input('name-dropdown', 'value'),Input('v_option', 'value'),Input('j_option', 'value')])
def update_graph_ppost_(ppost_value,pgen_value,index_,v_gene,j_gene):
    if v_gene is None or j_gene is None:
       return{
       'data': [],
            "layout": go.Layout(
                    paper_bgcolor = 'rgba(0,0,0,0)',
                    plot_bgcolor = 'rgba(0,0,0,0)',
                    xaxis = {
                    "side": "bottom",
                    "type": "linear",
                    "range": [data_tables[index_].binning_ppost.values.min(),data_tables[index_].binning_ppost.values.max()],
                    "title": " log10 Ppost",
                    },
                    showlegend=False,
                    yaxis= {
                        "side": "left",
                        "type": "linear",
                        "range": [0,data_tables[index_].ppost_data.values.max()+0.1]},
                    margin = {'t': 0,'b':100},
                    legend={'x':0,'y':1})
            }
    return {'data': [
                 {'x': data_tables[index_].binning_ppost.values, 'y': data_tables[index_].ppost_gen.values, 'type': 'line', 'name': 'pre-selection'},
                {'x': data_tables[index_].binning_ppost.values, 'y': data_tables[index_].ppost_data.values, 'type': 'line', 'name': 'post-selection'},
                 {'x': [np.log10((ppost_value+1e-50)), np.log10((ppost_value+1e-50))], 'y': [-1,4], 'type': 'line', 'name': 'your seq'},
            ],
            "layout": go.Layout(
                    paper_bgcolor = 'rgba(0,0,0,0)',
                    plot_bgcolor = 'rgba(0,0,0,0)',
                    xaxis = {
                    "side": "bottom",
                    "type": "linear",
                    "range": [data_tables[index_].binning_ppost.values.min(),data_tables[index_].binning_ppost.values.max()],
                    "title": " log10 Ppost",
                    },
                    showlegend=False,
                    yaxis= {
                        "side": "left",
                        "type": "linear",
                        "range": [0,data_tables[index_].ppost_data.values.max()+0.1]},
                    margin = {'t': 0,'b':100},
                    legend={'x':0,'y':1})
            }

@app.callback(Output('Q-graph', 'figure'), [Input('ppost-value', 'children'),Input('pgen-for-ppost', 'children'),Input('name-dropdown', 'value'),Input('v_option', 'value'),Input('j_option', 'value')])
def update_graph_ppost(ppost_value,pgen_value,index_,v_gene,j_gene):
    if pgen_value<1e-70 and ppost_value <1e-70 : ppost_value=1.
    if v_gene is None or j_gene is None:
        return {'data': [],
            "layout": go.Layout(
                    paper_bgcolor = 'rgba(0,0,0,0)',
                    plot_bgcolor = 'rgba(0,0,0,0)',
                    xaxis = {
                    "side": "bottom",
                    "type": "linear",
                    "range": [data_tables[index_].binning_Q.values.min(),data_tables[index_].binning_Q.values.max()],
                    "title": " log10 Q",
                    },
                    showlegend=False,
                    yaxis= {
                        "side": "left",
                        "type": "linear",
                        "range": [0,data_tables[index_].Q_data.values.max()+0.8]},
                    margin = {'t': 0,'b':100},
                    legend={'x':0,'y':1})
            }
    return {'data': [
                {'x': data_tables[index_].binning_Q.values, 'y': data_tables[index_].Q_gen.values, 'type': 'line', 'name': 'pre-selection'},
                {'x': data_tables[index_].binning_Q.values, 'y': data_tables[index_].Q_data.values, 'type': 'line', 'name': 'post-selection'},
                {'x': [np.log10((ppost_value+1e-50)/(pgen_value+1e-50)+1e-30), np.log10((ppost_value+1e-50)/(pgen_value+1e-50)+1e-30)], 'y': [-1,4], 'type': 'line', 'name': 'your seq'},
            ],
            "layout": go.Layout(
                    paper_bgcolor = 'rgba(0,0,0,0)',
                    plot_bgcolor = 'rgba(0,0,0,0)',
                    xaxis = {
                    "side": "bottom",
                    "type": "linear",
                    "range": [data_tables[index_].binning_Q.values.min(),data_tables[index_].binning_Q.values.max()],
                    "title": " log10 Q",
                    },
                    showlegend=False,
                    yaxis= {
                        "side": "left",
                        "type": "linear",
                        "range": [0,data_tables[index_].Q_data.values.max()+0.8]},
                    margin = {'t': 0,'b':100},
                    legend={'x':0,'y':1})
            }
@app.callback(Output('pgen-value', 'children'),
              [Input('submit-button', 'n_clicks')],[State('v_option', 'value'),State('j_option', 'value'),State('name-dropdown', 'value'),State('cdr3', 'value'),State('cdr3-type', 'value')])
def update_pgen(n_clicks,v_gene, j_gene, index,cdr3,type_cdr3):
    if n_clicks is None: return 1e-80
    else:
        time.sleep(1)
        qm=MinimalSonia(qfiles[index],norms[index][1])
        cdr3=str(cdr3).upper()
        if not bool(type_cdr3): cdr3=nt2aa(cdr3)
        return compute_pgen(index,[cdr3,v_gene,j_gene])

@app.callback([Output('ppost-value', 'children'),Output('pgen-for-ppost', 'children')],
              [Input('submit-ppost', 'n_clicks')],[State('v_option', 'value'),State('j_option', 'value'),State('name-dropdown', 'value'),State('cdr3', 'value'),State('cdr3-type', 'value')])
def update_ppost(n_clicks,v_gene, j_gene, index,cdr3,type_cdr3):
    if n_clicks is None: return 1e-80,1e-80
    else:
        #time.sleep(1)
        qm=MinimalSonia(qfiles[index],norms[index][1])
        cdr3=str(cdr3).upper()
        if not bool(type_cdr3): cdr3=nt2aa(cdr3)
        list_seqs=[]
        if v_gene is None and j_gene is None: #both
            genV,genJ=return_genes(int(index))
            for vgene in genV:
                for jgene in genJ:
                    list_seqs.append([str(cdr3).upper(), vgene, jgene])
        elif v_gene is None:
            genV,genJ=return_genes(int(index))
            for vgene in genV:
                list_seqs.append([str(cdr3).upper(), vgene, j_gene])
        elif j_gene is None:
            genV,genJ=return_genes(int(index))
            for jgene in genJ:
                list_seqs.append([str(cdr3).upper(), v_gene, jgene])
        else:
            list_seqs.append([str(cdr3).upper(), v_gene, j_gene])
        pgens= compute_pgen_for_ppost(index,list_seqs)
        qs=qm.compute_sel_factor(list_seqs)
        pposts=pgens*qs
        return pposts.sum(),pgens.sum()

@app.callback(Output('pgen-nucleotide', 'children'),
              [Input('submit-button', 'n_clicks')],[State('v_option', 'value'),State('j_option', 'value'),State('name-dropdown', 'value'),State('cdr3', 'value'),State('cdr3-type', 'value')])
def update_pgen_nt(n_clicks,v_gene, j_gene, index,cdr3,type_cdr3):
    if n_clicks is None: return 1e-80
    else:
        return compute_pgen_nt(index,str(cdr3).upper(), v_gene, j_gene,type_cdr3)

@app.callback(Output('print-pgen-nucleotide', 'children'),[Input('pgen-nucleotide', 'children')])
def update_pgen_nucleotide(pgen_val):
        if pgen_val < 1e-70 and pgen_val > 1e-90: return u''
        return u'''The Pgen of the nucleotide sequence is: {:.2e}'''.format(pgen_val)

@app.callback(Output('print-pgen-nucleotide', 'style'),[Input('pgen-nucleotide', 'children')],[State('cdr3-type', 'value')])
def update_pgen_nucloetide_display(n_clicks,cdr3_type):
        if cdr3_type<1: return  {'width': '100%','text-align': 'left'}
        else: return  {'display': 'none'}

@app.callback(Output('print-pgen', 'children'),[Input('pgen-value', 'children')])
def update_pgen_amin(pgen_val):
        if pgen_val < 1e-70 and pgen_val > 1e-90: return u''
        return u'''The Pgen of the aminoacid sequence is: {:.2e}'''.format(pgen_val)

@app.callback(Output('print-ppost', 'children'),[Input('ppost-value','children')])
def update_ppost_amin(ppost_val):
        if ppost_val < 1e-70 and ppost_val > 1e-90: return u''
        return u'''The Ppost of the aminoacid sequence is: {:.2e} '''.format(ppost_val)

@app.callback(Output('probability-to-observe', 'children'),
              [Input('n-productive','value'),Input('ppost-value', 'children')])
def update_output(n_prod,ppost_value):
    try:
        n=float(n_prod)
        final_value=(1-np.exp(-n*ppost_value)) # TO CHECK
    except:
        return ''
    if final_value< 1e-50: return ''
    return '= {:.1e}'.format(final_value)

@app.callback(Output('v_option', 'options'),[Input('name-dropdown', 'value')])
def update_date_dropdown(name):
    genV,genJ=return_genes(int(name))
    return [{'label': i, 'value': i} for i in genV]

@app.callback(Output('v_option', 'value'),[Input('name-dropdown', 'value'),Input('v_option', 'options')])
def update_value_v_dropdown(name,option):
    return default_seqs[name][2]

@app.callback(Output('j_option', 'value'),[Input('name-dropdown', 'value'),Input('j_option', 'options')])
def update_value_j_dropdown(name,option):
    return default_seqs[name][3]

@app.callback(Output('j_option', 'options'),[Input('name-dropdown', 'value')])
def update_date_dropdown_again(name):
    genV,genJ=return_genes(int(name))
    return [{'label': i, 'value': i} for i in genJ]


@app.callback(Output('my-link', 'style'), [Input('file', 'children')])
def show_download(value):
    if value is '': return { 'display': 'none'}
    return {'display': 'inline-block','text-align':'middle'}

@app.callback(Output('my-link', 'href'), [Input('generate-button', 'n_clicks'),Input('file', 'children')])
def update_download_link(n,file_csv):
    return file_csv

@app.callback(Output('file', 'children'), [Input('generate-button', 'n_clicks')],[State('n-generation','value'),State('name-dropdown-copy', 'value'),State('sampling-origin', 'value'),State('seed', 'value')])
def generate_seqs(n_clicks,n,index_,ppost,seed):
    if n_clicks is None or n is None: return ''
    n=int(n)
    time.sleep(1)
    seqs=sample_olga(n,index_,ppost,seed)
    df=pd.DataFrame(seqs,columns=['nucleotide','amino_acid','v_gene','j_gene'])
    csv_string = df.to_csv(index=False, encoding='utf-8')
    csv_string = "data:text/csv;charset=utf-8," + urllib.quote(csv_string)
    return csv_string

if __name__ == '__main__':
    app.run_server(debug=True)