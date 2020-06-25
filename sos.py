import dash_katex
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table 
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State
import io

import pandas as pd
from utils import *
import base64
import time
import urllib
import immprint

# import all the html structure
from layout import *

app = dash.Dash(__name__,
                meta_tags=[{"name": "viewport",
                            "content": "width=device-width, initial-scale=1"}],
                external_scripts=[
                    'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML',
                ])
app.title = "SOS"
server = app.server


@app.callback(
    [Output('cdr3', 'value'),
     Output('v_option', 'value'),
     Output('j_option', 'value')],
    [Input('cdr3-type', 'value'), Input('name-dropdown', 'value')])
def update_bar(cdr3_type, index):
    """
    Return the default cdr3/V/J either as nucleotide or amino-acid,
    depending on the user-choice in the form 'cdr3-type'.
    """
    default_seqs=[['CASSDAQGRNRGTEAFF','TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT','TRBV9','TRBJ1-1'],
            ['CAYRLIQGAQKLVF','TGTGCTTATAGATTAATTCAGGGAGCCCAGAAGCTGGTATTT','TRAV38-2/DV8','TRAJ54'],
            ['CARHYGSGSEWFDPW','TGTGCGAGACACTATGGTTCGGGGAGTGAATGGTTCGACCCCTGG','IGHV4-34','IGHJ5'],
            ['CMQALQTPQITF','TGCATGCAAGCTCTACAAACTCCTCAGATCACCTTC','IGKV2-28','IGKJ5'],
            ['CQSADSSGQVF','TGTCAATCAGCAGACAGCAGTGGACAGGTCTTC','IGLV3-25','IGLJ1'],
            ['CTCSLTSNERLFF','TGCACCTGCAGTTTGACATCCAACGAAAGATTATTTTTC','TRBV1','TRBJ1-4']]

    return ((default_seqs[index][0] if cdr3_type == "aa"
             else default_seqs[index][1]),
            default_seqs[index][2],
            default_seqs[index][3])


@app.callback([Output('print-pgen', 'children'),
               Output('print-pgen', 'style'),
               Output('pgen-graph', 'figure'),
               Output('hidden-pgen-value', 'children')],
              [Input('submit-button-pgen', 'n_clicks')],
              [State('v_option', 'value'),
               State('j_option', 'value'),
               State('name-dropdown', 'value'),
               State('cdr3', 'value'),
               State('cdr3-type', 'value')])
def update_pgen(n_clicks, v_gene, j_gene, index, cdr3, type_cdr3):
    """
    Compute pgen with the user provided elements (cdr3 + cdr3-type / V / J).
    Print the value of pgen and modify the graph.
    """
    if n_clicks is None:  # default behaviour
        return None, {"display": "none"}, figure_pgen, None
    else:
        qm = MinimalSonia(qfiles[index],norms[index][1])
        cdr3 = str(cdr3).upper()
        if type_cdr3 == "nt":
            cdr3 = nt2aa(cdr3)
        current_pgen = compute_pgen(index, [cdr3, v_gene, j_gene])
        fig = generate_graph(data_tables[0].binning_pgen,
                             data_tables[0].pgen_gen.values,
                             data_tables[0].pgen_data.values,
                             vline=current_pgen,
                             title=r"$\log_{10}(P_{\mathrm{gen}})$",
                             showlegend=True)
        return (dcc.Markdown(
            f"- **P<sub>gen</sub>** of the sequence: **{float_to_html(current_pgen)}**",
            dangerously_allow_html=True),
                {"display": "inline"},
                fig,
                current_pgen)
                

@app.callback([Output('print-ppost', 'children'),
               Output('print-ppost', 'style'),
               Output('Q-graph', 'figure'),
               Output('ppost-graph', 'figure'),
               Output('hidden-ppost-value', 'children')],
              [Input('submit-button-ppost', 'n_clicks')],
              [State('v_option', 'value'),
               State('j_option', 'value'),
               State('name-dropdown', 'value'),
               State('cdr3', 'value'),
               State('cdr3-type', 'value')])
def update_ppost(n_clicks, v_gene, j_gene, seq_type, cdr3, type_cdr3):
    """
    Compute Q and ppost with the user provided elements (cdr3 + cdr3-type / V / J).
    Print the value of pgen and modify the graph.
    """
    if n_clicks is None:
        return None, {"display": "none"}, figure_Q, figure_ppost, None
    qm = MinimalSonia(qfiles[seq_type],norms[seq_type][1])
    cdr3 = str(cdr3).upper()
    if type_cdr3 == "nt":
        cdr3 = nt2aa(cdr3)
    list_seqs = []

    # contrary to Olga, Sonia cannot "sum over all Vs/Js", so the sum
    # must be done manually (if some v/j are not provided)
    genV, genJ=return_genes(int(seq_type)) #  return all v/j genes
    if v_gene is None and j_gene is None:
        list_seqs = [(str(cdr3).upper(), v, j) for v in genV for j in genJ]
    elif v_gene is None:
        list_seqs = [(str(cdr3).upper(), v, j_gene) for v in genV]
    elif j_gene is None:
        list_seqs = [(str(cdr3).upper(), v_gene, j) for j in genJ]
    else:
        list_seqs = [(str(cdr3).upper(), v_gene, j_gene)]

    pgens = compute_pgen_for_ppost(seq_type, list_seqs)
    qs = qm.compute_sel_factor(list_seqs)
    pposts = pgens * qs
    current_ppost = pposts.sum()

    figppost = generate_graph(data_tables[0].binning_ppost,
                              data_tables[0].ppost_gen.values,
                              data_tables[0].ppost_data.values,
                              vline=current_ppost,
                              title=r"$\log_{10}(P_{\mathrm{post}})$")

    figQ = generate_graph(data_tables[0].binning_Q,
                          data_tables[0].Q_gen.values,
                          data_tables[0].Q_data.values,
                          vline=np.mean(qs[0]),
                          title=r"$\log_{10}(Q)$")
    
    return (dcc.Markdown(
        f"- **P<sub>post</sub>** of the sequence: **{float_to_html(current_ppost)}**",
        dangerously_allow_html=True),
            {"display": "inline"},
            figQ,
            figppost,
            current_ppost)


@app.callback(Output('probability-to-observe', 'expression'),
              [Input('n-productive','value'),
               Input('print-ppost', 'children'),
               Input('submit-button-ppost', 'n_clicks'),
               Input('hidden-ppost-value', 'children')])
def number_found(n_prod, ppost_value, nc, ppost):
    """
    Given the user input (number of total nucleotide sequences, modify
    the element 'probability-to-observe' to give the ppost probability to
    find a sequence. 
    """
    n=float(n_prod)
    # not clear on what to use, ppost, really ? 
    if ppost is not None:
        final_value=(1-np.exp(-n*ppost))
        if final_value > 1e-50:
            return (r'P_{\mathrm{obs}} = 1 - \left( 1 - P_{\mathrm{post}}\right)^N \approx 1 - e^{-N P_{\mathrm{post}}} = ' + float_to_tex(final_value))
    return ''


@app.callback(Output('v_option', 'options'),
              [Input('name-dropdown', 'value')])
def update_V_dropdown(name):
    genV, genJ=return_genes(int(name))
    return [{'label': i, 'value': i} for i in genV]

@app.callback(Output('j_option', 'options'),
              [Input('name-dropdown', 'value')])
def update_J_dropdown(name):
    genV, genJ=return_genes(int(name))
    return [{'label': i, 'value': i} for i in genJ]


@app.callback([Output('dl-link', 'href'),
               Output('dl-link', 'style'),
               Output('generate-results', 'children')],
              [Input('generate-button', 'n_clicks')],
              [State('n-generation','value'),
               State('name-dropdown-copy', 'value'),
               State('sampling-origin', 'value'),
               State('seed', 'value')])
def generate_seqs(n_clicks, nb_gen_seq, cdr3_type, sampling_origin, seed):
    """
    Generate `nb_gen_seq` of type given by the dropdown menu. 
    Generate a dl link and display the sequences on the right
    """
    if n_clicks is None or nb_gen_seq is None:
        return '', {'display': 'none'}, html.Div()
    nb_gen_seq = int(nb_gen_seq)
    seqs = sample_olga(nb_gen_seq, cdr3_type, sampling_origin, seed)
    df = pd.DataFrame(seqs, columns=['nucleotide', 'amino_acid', 'v_gene', 'j_gene'])
    csv_string = df.to_csv(index=False, encoding='utf-8')
    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(csv_string)
    id_to_name = {'nucleotide': 'Nucleotide CDR3',
                  'amino_acid': 'Amino-Acid CDR3',
                  'v_gene': 'V', 'j_gene': 'J'}
    table_display = html.Div([html.H5("Generated sequences"),
                              dash_table.DataTable(
                                  data=df.to_dict('records'),
                                  columns=[{'name':id_to_name[i], 'id':i} for i in df.columns],
                                  style_cell={'textAlign': 'center'},
                                  style_table={
                                      'maxWidth': '100%',
                                      'overflowX': 'scroll',
                                      'maxHeight': '10cm',
                                      'overflowY': 'scroll'
                                      },
                                  style_data_conditional=[
                                      {
                                          'if': {'row_index': 'odd'},
                                          'backgroundColor': 'rgb(248, 248, 248)'
                                      }],
                                  style_as_list_view=True,
                                  style_header={
                                      'backgroundColor': 'white',
                                      'fontWeight': 'bold'
                                  },
                                  css=[
                                      {
                                          'selector': 'table',
                                          'rule': 'width: 100%;'
                                      }
                                  ]),
                              html.Hr()],
                             className='table-container')

    return csv_string, {'display': 'inline'}, table_display

@app.callback([Output('dl-A', 'href'),
               Output('dl-A', 'children')],
              [Input('upload-dataset-A', 'contents'),
               Input('radio-full', 'value')],
              [State('upload-dataset-A', 'filename')])
def upload_sample_A(contentsA, radio_full, filenameA):
    """
    When a sample is uploaded, hide its content in a div on the user side.
    """
    if contentsA is None:
        if radio_full == 'trb':
            return csv_string_sampleA, 'example A.csv'
        elif radio_full == 'full':
            return csv_string_sampleA_full, 'example A (full).csv'
        else:
            return '', ''
    return contentsA, filenameA


@app.callback([Output('dl-B', 'href'),
               Output('dl-B', 'children')],
              [Input('upload-dataset-B', 'contents'),
               Input('radio-full', 'value')],
              [State('upload-dataset-B', 'filename')])
def upload_sample_B(contentsB, radio_full, filenameB):
    """
    When a sample is uploaded, hide its content in a div on the user side.
    """
    if contentsB is None:
        if radio_full == 'trb':
            return csv_string_sampleB, 'example B.csv'
        elif radio_full == 'full':
            return csv_string_sampleB_full, 'example B (full).csv'
        else:
            return '', ''
    return contentsB, filenameB




@app.callback([Output('immprint-answer', 'children'),
               Output('immprint-results', 'children')],
              [Input('immprint-submit-button', 'n_clicks')],
              [State('dl-A', 'href'),
               State('dl-B', 'href'),
               State('radio-full', 'value'),
               State('radio-estimator', 'value')])
def immprit_web(nclicks, contentsA, contentsB, radio_full, radio_estimator):
    """ Main logic of the immprint webapp, load the two dataframes, 
        apply immprint and display the results 
    """

    if contentsA is None or contentsB is None:
        return '', ''

    _, contentsA = contentsA.split(',')
    decoded = base64.b64decode(contentsA)
    try:
        dfA = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=None, engine='python')
    except Exception as e:
        print(e)
        return html.Div(['There was an error processing the sample A file.']), ''

    _, contentsB = contentsB.split(',')
    decoded = base64.b64decode(contentsB)
    try:
        dfB = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=None, engine='python')
    except Exception as e:
        print(e)
        return html.Div(['There was an error processing the sample B file.']), ''

    if radio_full not in ['full', 'trb']:
        return html.Div(['There was an error.']), ''
    full = True if radio_full == 'full' else False
    
    dfA = immprint.rename_cols(dfA)
    dfB = immprint.rename_cols(dfB)
    
    S, I, shared_sequences, parms = immprint.immprint(dfA, dfB,
                                                      full=full,
                                                      use_counts=True,
                                                      onlyS=(radio_estimator == 's'),
                                                      γ=12)

    under_threshold = (S < parms['rS'] if I is None else I < parms['rI'])
    ambiguous = (under_threshold and parms['pv1'] > 1e-4) or (not under_threshold and parms['pv2'] > 1e-4)
    immprint_says = ""
    if ambiguous:
        immprint_says = ("Ambiguous case: maybe not enough sequences.\n"
                         f"- Probability the two samples come from the same patient <= {float_to_html(parms['pv1'])}\n"
                         f"- Probability the two samples come from different patients <= {float_to_html(parms['pv2'])}\n")
    else:
        if under_threshold: # different patients
            immprint_says = (f"**The samples come from two different patients.**\n"
                      f"P-value lower than {float_to_html(max(1e-100, parms['pv1']))}")
        if not under_threshold: # same patient
            immprint_says = (f"**The samples come from the same patient.**\n"
                             f"P-value lower than {float_to_html(max(1e-100, parms['pv2']))}")


    figure = None
    if I is None:
        xS1, yS1 = immprint.S_graph(parms['µS1'])
        xS2, yS2 = immprint.S_graph(parms['µS2'])
        figure = generate_graph_immprint(xS1, xS2, yS1, yS2, parms["rS"], S, title="S")
    else:
        xI1, yI1 = immprint.I_graph(parms["µS1"], parms["µ_logpgen"], parms["σ_logpgen"])
        xI2, yI2 = immprint.I_graph(parms["µS2"],
                                    parms["µ_logpgen_shared"],
                                    parms["σ_logpgen_shared"])
        figure = generate_graph_immprint(xI1, xI2, yI1, yI2, parms["rI"], I, title="I")
        
    answer = [html.P(f"Number of shared sequences S: {S}")]
    if I is not None:
        answer += [html.P(f"Immprint estimate I: {I:.2f}")]

    answer += [dcc.Markdown(immprint_says, dangerously_allow_html=True)]

    plots = [html.H6('Probability distributions')]
    plots += [dcc.Graph(id='plot-immprint',
                         config={'staticPlot': True},
                         figure=figure,
                         className='graph')]
    
    plots += [dash_table.DataTable(
        data=[{'Shared CDR3': s} for s in shared_sequences],
        columns=[{'name': 'Shared CDR3', 'id': 'Shared CDR3'}],
        style_table={
            'maxWidth': '100%',
            'overflowX': 'scroll',
            'maxHeight': '6cm',
            'overflowY': 'scroll'
        }),
        html.Hr()]

    return answer, plots


# app layout
# 4 tabs, each of them containing two windows left/right,
# one for interacting, the other for displaying results

with open("assets/about.md") as f:
    content = dcc.Markdown(f.read(), dangerously_allow_html=True)
    tab_about = dcc.Tab(label='About',
                        value='what-is',
                        className='top-tab',
                        children=[html.Div(content, className='interaction-div'),
                                  html.Div(
                                      html.Img(src='assets/explanation.png'),
                                      className='results-div')])

tab_olga = dcc.Tab(label='Evaluate',
                   className='top-tab',
                   children=[html.Div(layout_evaluate_tab,
                                      className='interaction-div'),
                             html.Div(pgen_results,
                                      className='results-div')])

tab_generate = dcc.Tab(label='Generate',
                       className='top-tab',
                       children=[html.Div(layout_generate_tab,
                                        className='interaction-div'),
                                 html.Div(generate_results,
                                          className='results-div')])

tab_immprit = dcc.Tab(label='ImmPrint',
                      className='top-tab',
                      children=[html.Div(layout_immprit_tab,
                                         className='interaction-div'),
                                html.Div(immprint_results,
                                         className='results-div')])

tabs = dcc.Tabs(value='what-is', className='tabs',
                children = [tab_about, tab_olga, tab_generate, tab_immprit],
                colors={'primary':'gold',"background": "cornsilk"})

app.layout = html.Div([hidden_layout, tabs], id='main')



if __name__ == '__main__':
    # This is used when running locally
    app.run_server(host='127.0.0.1', port=8080, debug=True)
