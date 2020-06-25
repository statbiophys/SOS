import dash_katex
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import pandas as pd
import urllib
from utils import *
import base64


def generate_graph(bins, gen, data, vline=None, title='', showlegend=False):
    """ 
    Plot a graph comparing two curves + potentially a vertical line
    Return a dictionary with keys "data" and "layout", that can be inserted in a dcc.Graph object.
    """
    figure = {'data': [], 'layout': None}
    figure['data'] += [{'x': bins, 'y': gen, 'type': 'line', 'name': 'Pre-selection'}]
    figure['data'] += [{'x': bins, 'y': data, 'type': 'line', 'name': 'Post-selection'}]
    if vline is not None:
        figure['data'] += [{'x': [np.log10(vline)]*2, 'y': [-1, 2*data.max()],
                            'type': 'name', 'name': 'Your sequence'}]
    figure['layout'] = go.Layout(
        paper_bgcolor = 'rgba(0,0,0,0)',
        plot_bgcolor = 'rgba(0,0,0,0)',
        xaxis = {
            "side": "bottom",
            "type": "linear",
            "range": [bins.min(), bins.max()],
            "title": title,
            "zeroline": False,
        },
        showlegend=showlegend,
        yaxis= {
            "side": "left",
            "type": "linear",
            "range": [0, 1.1 * data.max()],
        },
        margin = {'t': 0,'b':100},
        legend={'x':0,'y':1, 'bgcolor':'white'})
    return figure


def generate_graph_immprint(x1, x2, y1, y2, thresh, val, title=''):
    """ Plot one of the S/I immprint graph and return the figure object """
    def scale(x):
        return np.log10(x + 1)
    
    figure = {'data': [], 'layout': None}
    figure['data'] += [{'x': scale(x1), 'y': y1/np.max(y1), 'type': 'line',
                        'name': 'Same donor', 'line': {'width': 3} }]
    figure['data'] += [{'x': scale(x2), 'y': y2/np.max(y2), 'type': 'line',
                        'name': 'Different donors', 'line': {'width': 3} }]
    figure['data'] += [{'x': [scale(thresh)]*2, 'y': [-1, 2],
                        'type': 'name', 'name': 'Threshold', 'line': {'width': 3} }]
    figure['data'] += [{'x': [scale(val)]*2, 'y': [-1, 2],
                        'type': 'name', 'name': 'Your samples', 'line': {'width': 3} }]
        
    figure['layout'] = go.Layout(
        paper_bgcolor = 'rgba(0,0,0,0)',
        plot_bgcolor = 'rgba(0,0,0,0)',
        xaxis = {
            "tickmode":'array',
            "tickvals":[scale(1), scale(10), scale(100)],
            "ticktext":["1", "10", "100"],
            "side": "bottom",
            "type": "linear",
            "range": [0, max(scale(x1).max(), scale(x2).max())],
            "title": title,
            "zeroline": False,
        },
        showlegend=True,
        yaxis= {
            "side": "left",
            "type": "linear",
            "showticklabels": False,
            "range": [0, 1.1],
        },
        margin = {'t': 0,'b':100},
        legend={'x':0.5,'y':1, 'bgcolor': 'white'})
    return figure


test_base64 = base64.b64encode(open(local_directory+'assets/plots.png', 'rb').read()).decode('ascii')

# load the sample cdr3s for Olga/Sonia
data_tables=[pd.read_csv(local_directory+'examples/data_'+i+'.csv') for i in options_of]

# load the test samples for Immprint
# the encoding horror that follows seems to be necessary
df_sampleA = pd.read_csv('examples/sampleA.csv')
csv_string_sampleA = "data:text/csv;base64," +  base64.b64encode(df_sampleA.to_csv(index=False, encoding='utf-8').encode('utf-8')).decode("UTF-8")
df_sampleA_full = pd.read_csv('examples/sampleA_full.csv')
csv_string_sampleA_full = "data:text/csv;base64," +  base64.b64encode(df_sampleA_full.to_csv(index=False, encoding='utf-8').encode('utf-8')).decode("utf-8")

df_sampleB = pd.read_csv('examples/sampleB.csv')
csv_string_sampleB = "data:text/csv;base64," +  base64.b64encode(df_sampleB.to_csv(index=False, encoding='utf-8').encode('utf-8')).decode("UTF-8")
df_sampleB_full = pd.read_csv('examples/sampleB_full.csv')
csv_string_sampleB_full = "data:text/csv;base64," +  base64.b64encode(df_sampleB_full.to_csv(index=False, encoding='utf-8').encode('utf-8')).decode("UTF-8")

# padding p
padP = html.P(style={"height": "1%"})

# store values to share between callback
# (bad idea to use global python values)
hidden_layout = html.Div([html.Div(id='hidden-pgen-value'),
                          html.Div(id='hidden-ppost-value')], className='hidden')


## Layout of the tab where pgen/ppost are computed ##

layout_evaluate_tab = [html.H4('Evaluate the probability'),

                       html.P(''' Select the model'''),
                       dcc.Dropdown(
                           id='name-dropdown',
                           className="dropdown-menu",
                           options=[{'label':name, 'value':i} for i, name in enumerate(names)],
                           value=0, clearable=False),
                       padP,
                       html.P('''Input CDR3,V and J genes '''),
                       html.Div([dcc.Dropdown(id='cdr3-type',
                                              options=[{'label': 'Amino-Acid', 'value': "aa"},
                                                       {'label': 'Nucleotide', 'value': "nt"}],
                                              value='aa', clearable=False)],
                                className="dropdown-menu", style={'width':'48%'}),
                       html.P(style={'width': '4%'}),
                       dcc.Input(id='cdr3',
                                 value='CASSDAQGRNRGTEAFF',
                                 type='text', className="input", style={"width":'48%'}),
                       html.P(),
                       html.Div([dcc.Dropdown(id='v_option', clearable=True, value='TRBV9')],
                                className="dropdown-menu", style={"width": "48%"}),
                       html.P(style={'width': '4%'}),
                       html.Div([dcc.Dropdown(id='j_option', clearable=True, value='TRBJ1-1')],
                                className="dropdown-menu", style={"width": "48%"}),
                       html.Div(
                           [html.Button(id='submit-button-pgen',
                                       children='PGEN',
                                       title=('olga-compute_pgen')),
                            html.P(style={"width": "10%"}),
                            html.Button(id='submit-button-ppost',
                                        children='PPOST',
                                        title=('sonia-evaluate'))],
                           className='centered'),
                       html.Button(className="hidden", id='fake-button', children='Compute'),
                       padP,
                       html.Div(id='model-type',children='human_T_beta', className='hidden'),
                       dcc.Loading(id="loading-2",
                                   children=[
                                       html.P(id='print-pgen', children=''),
                                       html.P(id='print-ppost', children='')
                                       ]
                                   ),
                       html.Div([
                           html.P('''How many unique productive nucleotide sequences do you have?'''),
                           dcc.Input(id='n-productive', value=100000, type='number', className="input")
                       ]),
                       html.Br(),
                       html.Div([
                           html.P('Probability to find the aminoacid '
                                  'sequence by pure chance (p-value):'),
                           dash_katex.DashKatex(
                               id='probability-to-observe',
                               expression=r'P_{\mathrm{obs}} = 1 - \left( 1 - P_{\mathrm{post}}\right)^N \approx 1 - e^{-N P_{\mathrm{post}}}'),
                       #    html.P(id='probability-to-observe',children='')
                       ])]


## Layout of the tab to generate sequences ##

layout_generate_tab = [html.H4(children='Generate synthetic repertoires'),
                       html.P(children=''' Select the model'''),
                       dcc.Dropdown(
                           id='name-dropdown-copy',
                           options=[{'label':name, 'value': i} for i,name in enumerate(names)],
                           value=0, clearable=False),
                       html.P(children='Generate synthetic sequences from the model (max 1000)'),
                       
                       html.Div([
                           html.Div(dcc.Dropdown(id='sampling-origin',
                                              options=[{'label': 'pre-selection', 'value': 0},
                                                       {'label': 'post-selection', 'value': 1}],
                                                 value=0, clearable=False),
                                    style={'width': '32%'},
                                    className='dropdown-menu'),
                           html.P(' ', style={'width': '2%'}),
                           dcc.Input(id='n-generation',
                                     placeholder="n sequences",
                                     value=10,
                                     type='number',
                                     style={'width': '32%'},
                                     className="input"),
                           html.P(' ', style={'width': '2%'}),
                           dcc.Input(id="seed", type="number",
                                     placeholder="seed",
                                     style={'width': '32%'},
                                     className="input")],
                                style={'display': 'flex'}),
                       html.Div(
                       html.Button(id='generate-button',
                                   children='GENERATE',
                                   style={'width': '50%'},
                                   title='sonia-generate --humanTRB -n 1000 --pgen (or --ppost)'), className='centered'),
                       dcc.Loading(id="loading-1",
                                   children= [html.P(id='file', children='', className="hidden"),
                                              html.Div(html.A('download', id='dl-link',
                                                              className='hidden'),
                                                       className="centered")],
                                   type="default")]

## Layout of the immprint tab ##

layout_immprit_tab = [
    dcc.Markdown("""
    #### ImmPrint (TCR only)

    ImmPrint answers the question "**were these two samples extracted from the same donor ?**" by looking at the numbers and sequences of shared unique TCR between the two samples.
        """),
    padP,
    html.P('Receptor type', style={'width': '45%'}),
    html.P('Estimator', style={'width': '45%'}),
    html.Div(dcc.RadioItems(options=[
        {'label': 'TRB only', 'value': 'trb'},
        {'label': 'Full receptor', 'value': 'full'}],
                   value='trb', id='radio-full'), style={'width': '45%', 'display':'inline-block'}),
    html.Div(dcc.RadioItems(
        options=[{'label': 'Number of shared sequences 𝓢', 'value': 's'},
                 {'label': 'ImmPrint pgen-based estimator 𝓘 (can be slow)', 'value': 'i'}],
        value='s', id='radio-estimator'), style={'width':'45%', 'display':'inline-block'}),
    padP,
    html.Div([
        html.P('Sample A'),
        html.Div(
            dcc.Upload(id='upload-dataset-A',
                         className='upload',
                         children=html.Div(
                             ['Drag and Drop or ',
                             html.A(id='link_select_A',
                                    children=' select a csv file')]
                         ), style={ # can't put this in a css it seems
                             'width': '80%',
                             'height': '50px',
                             'lineHeight': '60px',
                             'borderWidth': '1px',
                             'borderStyle': 'dashed',
                             'borderRadius': '5px',
                             'textAlign': 'center',
                             'margin': '10px'
                         }),
            style={'width': '75%', 'display': 'inline-block'}),
        html.Div(html.A('example_P1.csv', id='dl-A'),
                 style={'width': '21%', 'display': 'inline-block'}),
    ]),
    html.Div([
        html.P('Sample B'),
        html.Div(
            dcc.Upload(id='upload-dataset-B',
                         className='upload',
                         children=html.Div(
                             ['Drag and Drop or ',
                             html.A(id='link_select_B',
                                    children=' select a csv file')]
                         ), style={ # can't put this in a css it seems
                             'width': '80%',
                             'height': '50px',
                             'lineHeight': '60px',
                             'borderWidth': '1px',
                             'borderStyle': 'dashed',
                             'borderRadius': '5px',
                             'textAlign': 'center',
                             'margin': '10px'
                         }), className='draganddrop', style={'width': '75%', 'display': 'inline-block'}
        ),
        html.Div(html.A('example_P2.csv', id='dl-B'),
                 style={'width': '21%', 'display': 'inline-block'})
    ]),
    dcc.Loading(
        html.Div(html.Button(id='immprint-submit-button',
                            children='Immprint',
                            title=('immprint sampleA.csv sampleB.csv')),
                className='centered')
    ),
    padP,
    html.Div(id='immprint-answer')
]


## generate the default plots for the tabs ##

figure_Q = generate_graph(data_tables[0].binning_Q.values,
                          data_tables[0].Q_gen.values,
                          data_tables[0].Q_data.values,
                          title=r"$\log_{10}(Q)$")
figure_ppost = generate_graph(data_tables[0].binning_ppost.values,
                              data_tables[0].ppost_gen.values,
                              data_tables[0].ppost_data.values,
                          title=r"$\log_{10}(P_{\mathrm{post}})$")
figure_pgen = generate_graph(data_tables[0].binning_pgen,
                              data_tables[0].pgen_gen.values,
                             data_tables[0].pgen_data.values,
                             title=r"$\log_{10}(P_{\mathrm{gen}})$",
                             showlegend=True)

plots = [html.Div(html.Img(src='assets/description_plots.png',
                           id='descr-graph'),
                  className='plot-div')]
plots += [html.Div(dcc.Graph(id=idd,
                            config={'staticPlot': True},
                            figure=fig,
                            className='graph'),
                  className='plot-div')
         for idd, fig in [('Q-graph', figure_Q), ('ppost-graph', figure_ppost), ('pgen-graph', figure_pgen)]]

pgen_results = [html.Div(plots, id='pgen-results')]
generate_results = [html.Div('', id='generate-results')]
immprint_results = [dcc.Loading(html.Div( '', id='immprint-results'))]


