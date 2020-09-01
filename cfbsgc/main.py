# -*- coding: utf-8 -*-

from functools import lru_cache

from os import listdir
from os.path import dirname, join

import numpy as np
import pandas as pd
import pickle
import sqlite3

from bokeh.io import curdoc
from bokeh.layouts import row, column, widgetbox, layout, gridplot
from bokeh.models import ColumnDataSource, CustomJS, HoverTool, Div, Circle
from bokeh.palettes import Spectral6
from bokeh.transform import linear_cmap
from bokeh.models.widgets import Select, CheckboxGroup
from bokeh.plotting import figure

from io import StringIO
import base64

from sklearn import preprocessing
min_max_scaler = preprocessing.MinMaxScaler()

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# creat dictionaries for plotting selections
with open(join(dirname(__file__), 'axis_map.pkl'), 'rb') as handle:
    axis_map = pickle.load(handle)
with open(join(dirname(__file__), 'reverse_axis_map.pkl'), 'rb') as handle:
    reverse_axis_map = pickle.load(handle)
with open(join(dirname(__file__), 'library_map.pkl'), 'rb') as handle:
    library_map = pickle.load(handle)

# open database connection
conn = sqlite3.connect(join(dirname(__file__), 'commercial_fragment_libraries.sqlite3'))
cursor = conn.cursor()
fields = cursor.fetchall()
column_names = [x[1] for x in fields]

pca_map = {"Autocorr2D" : "CalcAUTOCORR2D-pca",
            "Autocorr3D" : "CalcAUTOCORR3D-pca",
            "GETAWAY" : "CalcGETAWAY-pca",
            "MORSE" : "CalcMORSE-pca",
            "RDF" : "CalcRDF-pca",
            "WHIM" : "CalcWHIM-pca",
            "Atom Pair Fingerprint" : "GetHashedAtomPairFingerprintAsBitVect-pca",
            "Topological Torsion Fingerprint" : "GetHashedTopologicalTorsionFingerprintAsBitVect-pca",
            "MACCS Keys Fingerprint" : "GetMACCSKeysFingerprint-pca",
            "Layered Fingerprint" : "LayeredFingerprint-pca",
            "Pattern Fingerprint" : "PatternFingerprint-pca",
            "RDKit Fingerprint" : "RDKFingerprint-pca",
            "Mol2Vec Features" : "mol2vec-pca"}


pca = "mol2vec-pca"
lib = ""

# Create selection controls
color_by = Select(title="Color By Values", options=sorted(list(axis_map.keys())), value='MolWt')
alpha_by = Select(title="Alpha by Values", options=sorted(list(axis_map.keys())), value='MolLogP')
size_by = Select(title="Size by Values", options=sorted(list(axis_map.keys())), value='NumRotatableBonds')
pca_select = Select(title="Principle Component Analysis", options=sorted(list(pca_map.keys())), value='RDKit Fingerprint')

# Create library selection check boxes
library_names = list(library_map.keys())
ordered_library_names = list(np.sort(library_names))
library = CheckboxGroup(labels=ordered_library_names, active=[])

# create a column data source for the plots to share
compounds = pd.DataFrame(columns=['ID', axis_map[color_by.value], axis_map[alpha_by.value], axis_map[size_by.value], 'png', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5'])

source = ColumnDataSource(data=dict(PC1=[], PC2=[], PC3=[], PC4=[], PC5=[], color=[], alpha=[], size=[], png=[]))

# Create scatter plot
TOOLTIPS=[
        ("size_real", "@size_real"),
        ("alpha_real", "@alpha_real"),
        ("color_real", "@color_real"),
        ("png", "@png")
]

hover = HoverTool(tooltips="""
        <div>
                <div>
                        <img
                                src="@png" height="150" alt="@png" width="150"
                                style="float: left; margin: 0px 15px 15px 0px;"
                                border="0"
                        ></img>
                </div>
                <div>color = @color_real</div>
                <div>alpha = @alpha_real</div>
                <div>size = @size_real</div>
        </div>
        """
)

TOOLS = ["pan, box_select, wheel_zoom, box_zoom, reset, save, help", hover]

selected1= Circle(fill_alpha=1, fill_color="grey", line_color=None)
nonselected1 = Circle(fill_alpha='alpha', fill_color='lightblue', line_color=None)

selected2 = Circle(fill_alpha=1, fill_color="grey", line_color=None)
nonselected2 = Circle(fill_alpha=0.2, fill_color="lightblue", line_color=None)

selected3 = Circle(fill_alpha=1, fill_color="grey", line_color=None)
nonselected3 = Circle(fill_alpha=0.2, fill_color="lightblue", line_color=None)

selected4 = Circle(fill_alpha=1, fill_color="grey", line_color=None)
nonselected4 = Circle(fill_alpha=0.2, fill_color="lightblue", line_color=None)

selected5 = Circle(fill_alpha=1, fill_color="grey", line_color=None)
nonselected5 = Circle(fill_alpha=0.2, fill_color="lightblue", line_color=None)

selected6 = Circle(fill_alpha=1, fill_color="grey", line_color=None)
nonselected6 = Circle(fill_alpha=0.2, fill_color="lightblue", line_color=None)

selected7 = Circle(fill_alpha=1, fill_color="grey", line_color=None)
nonselected7 = Circle(fill_alpha=0.2, fill_color="lightblue", line_color=None)

selected8 = Circle(fill_alpha=1, fill_color="grey", line_color=None)
nonselected8 = Circle(fill_alpha=0.2, fill_color="lightblue", line_color=None)

selected9 = Circle(fill_alpha=1, fill_color="grey", line_color=None)
nonselected9 = Circle(fill_alpha=0.2, fill_color="lightblue", line_color=None)

selected10 = Circle(fill_alpha=1, fill_color="grey", line_color=None)
nonselected10 = Circle(fill_alpha=0.2, fill_color="lightblue", line_color=None)

# create a new plot
s1 = figure(plot_width=250, plot_height=250, title='PC1 vs PC2', tools=TOOLS, toolbar_location='left')
s1.circle('PC1', 'PC2', size='size', alpha='alpha', color='color', source=source, name="circles1")

s1renderer = s1.select(name="circles1")
s1renderer.selection_glyph = selected1
s1renderer.nonselection_glyph = nonselected1

# create a new plot and share both ranges
s2 = figure(plot_width=250, plot_height=250, x_range=s1.x_range, title='PC1 vs PC3', tools=TOOLS)
s2.circle('PC1', 'PC3', size='size', alpha='alpha', color='color', source=source, name="circles2")

s2renderer = s2.select(name="circles2")
s2renderer.selection_glyph = selected2
s2renderer.nonselection_glyph = nonselected2

# create a new plot and share both ranges
s3 = figure(plot_width=250, plot_height=250, x_range=s1.x_range, title='PC1 vs PC4', tools=TOOLS)
s3.circle('PC1', 'PC4', size='size', alpha='alpha', color='color', source=source, name="circles3")

s3renderer = s3.select(name="circles3")
s3renderer.selection_glyph = selected3
s3renderer.nonselection_glyph = nonselected3

# create a new plot and share both ranges
s4 = figure(plot_width=250, plot_height=250, x_range=s1.x_range, title='PC1 vs PC5', tools=TOOLS)
s4.circle('PC1', 'PC5', size='size', alpha='alpha', color='color', source=source, name="circles4")

s4renderer = s4.select(name="circles4")
s4renderer.selection_glyph = selected4
s4renderer.nonselection_glyph = nonselected4

# create a new plot and share both ranges
s5 = figure(plot_width=250, plot_height=250, x_range=s1.y_range, y_range=s2.y_range, title='PC2 vs PC3', tools=TOOLS)
s5.circle('PC2', 'PC3', size='size', alpha='alpha', color='color', source=source, name="circles5")

s5renderer = s5.select(name="circles5")
s5renderer.selection_glyph = selected5
s5renderer.nonselection_glyph = nonselected5

# create a new plot and share both ranges
s6 = figure(plot_width=250, plot_height=250, x_range=s1.y_range, y_range=s3.y_range, title='PC2 vs PC4', tools=TOOLS)
s6.circle('PC2', 'PC4', size='size', alpha='alpha', color='color', source=source, name="circles6")

s6renderer = s6.select(name="circles6")
s6renderer.selection_glyph = selected6
s6renderer.nonselection_glyph = nonselected6

# create a new plot and share both ranges
s7 = figure(plot_width=250, plot_height=250, x_range=s1.y_range, y_range=s4.y_range, title='PC2 vs PC5', tools=TOOLS)
s7.circle('PC2', 'PC5', size='size', alpha='alpha', color='color', source=source, name="circles7")

s7renderer = s7.select(name="circles7")
s7renderer.selection_glyph = selected7
s7renderer.nonselection_glyph = nonselected7

# create a new plot and share both ranges
s8 = figure(plot_width=250, plot_height=250, x_range=s2.y_range, y_range=s3.y_range, title='PC3 vs PC4', tools=TOOLS)
s8.circle('PC3', 'PC4', size='size', alpha='alpha', color='color', source=source, name="circle8")

s8renderer = s8.select(name="circle8")
s8renderer.selection_glyph = selected8
s8renderer.nonselection_glyph = nonselected8

# create a new plot and share both ranges
s9 = figure(plot_width=250, plot_height=250, x_range=s2.y_range, y_range=s4.y_range, title='PC3 vs PC5', tools=TOOLS)
s9.circle('PC3', 'PC5', size='size', alpha='alpha', color='color', source=source, name="circles9")

s9renderer = s9.select(name="circles9")
s9renderer.selection_glyph = selected9
s9renderer.nonselection_glyph = nonselected9

# create a new plot and share both ranges
s10 = figure(plot_width=250, plot_height=250, x_range=s3.y_range, y_range=s4.y_range, title='PC4 vs PC5', tools=TOOLS)
s10.circle('PC4', 'PC5', size='size', alpha='alpha', color='color', source=source, name="circles10")

s10renderer = s10.select(name="circles10")
s10renderer.selection_glyph = selected10
s10renderer.nonselection_glyph = nonselected10

global p
p = gridplot([[s1, s2, s3, s4],
             [s5, s6, s7],
             [s8, s9, s10]])

# bokeh plottomg methods
def select_library(attr, old, new):
    global compounds
    for i in range(len(library_names)):
        ix = library_names.index((ordered_library_names[i]))
        substring = '-'.join(list(library_map.values())[ix].split('-')[:-1])+'-'
        if i not in library.active:
            compounds = compounds[compounds.ID.str.contains(substring) == False]
        if i in library.active:
            strings = ['-'.join(j.split('-')[:-1])+"-" for j in compounds.ID]
            if substring not in strings:
                resultDf = pd.read_sql_query('SELECT ID,zinc FROM data WHERE ID LIKE "'+substring+'%"', conn)
                png = [str]*len(resultDf)
                for z, ix in zip(resultDf.zinc, range(len(resultDf.zinc))):
                    png[ix] = "http://zinc15.docking.org/substances/"+z+".png"
                resultDf['png'] = png
                resultDf = resultDf[['ID', 'png']]
                tmpDf = pd.read_sql_query('SELECT ID,'+axis_map[color_by.value]+','+axis_map[alpha_by.value]+','+axis_map[size_by.value]+' FROM data WHERE ID LIKE"'+substring+'%"', conn)
                resultDf = pd.merge(resultDf, tmpDf, on=['ID'], sort=True)
                tmpDf = pd.read_pickle(join(dirname(__file__), pca_map[pca_select.value]+".pkl"))
                resultDf = pd.merge(resultDf, tmpDf, on=['ID'], sort=True)
                compounds = compounds.append(resultDf, sort=True)
    compounds = compounds.dropna()
    update()

def select_color(attr, old, new):
    global compounds
    att = color_by.value
    if att not in list(compounds):
        resultDf = pd.DataFrame(columns=['ID', axis_map[att]])
        for i in range(len(library_names)):
            ix = library_names.index(ordered_library_names[i])
            substring = '-'.join(list(library_map.values())[ix].split('-')[:-1])+'-'
            selected = axis_map[att]
            if i not in library.active:
                compounds = compounds[compounds.ID.str.contains(substring) == False]
            if i in library.active:
                strings = ['-'.join(j.split('-')[:-1])+"-" for j in compounds.ID]
                tmpDf = pd.read_sql_query('SELECT ID,"'+selected+'" FROM data WHERE ID LIKE "'+substring+'%"', conn)
                resultDf = resultDf.append(tmpDf)
        actives = [axis_map[color_by.value], axis_map[alpha_by.value], axis_map[size_by.value], 'ID', 'png', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
        for l in list(compounds):
            if l not in actives:
                compounds = compounds.drop(l, 1)
        compounds = pd.merge(compounds, resultDf, on=['ID'], sort=True)
        compounds = compounds.dropna()
    update()

def select_alpha(attr, old, new):
    global compounds
    att = alpha_by.value
    if att not in list(compounds):
        resultDf = pd.DataFrame(columns=['ID', axis_map[att]])
        for i in range(len(library_names)):
            ix = library_names.index(ordered_library_names[i])
            substring = '-'.join(list(library_map.values())[ix].split('-')[:-1])+'-'
            selected = axis_map[att]
            if i not in library.active:
                compounds = compounds[compounds.ID.str.contains(substring) == False]
            if i in library.active:
                strings = ['-'.join(j.split('-')[:-1])+"-" for j in compounds.ID]
                tmpDf = pd.read_sql_query('SELECT ID,"'+selected+'" FROM data WHERE ID LIKE "'+substring+'%"', conn)
                resultDf = resultDf.append(tmpDf)
        actives = [axis_map[color_by.value], axis_map[alpha_by.value], axis_map[size_by.value], 'ID', 'png', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
        for l in list(compounds):
            if l not in actives:
                compounds = compounds.drop(l, 1)
        compounds = pd.merge(compounds, resultDf, on=['ID'], sort=True)
        compounds = compounds.dropna()
    update()

def select_size(attr, old, new):
    global compounds
    att = size_by.value
    if att not in list(compounds):
        resultDf = pd.DataFrame(columns=['ID', axis_map[att]])
        for i in range(len(library_names)):
            ix = library_names.index(ordered_library_names[i])
            substring = '-'.join(list(library_map.values())[ix].split('-')[:-1])+'-'
            selected = axis_map[att]
            if i not in library.active:
                compounds = compounds[compounds.ID.str.contains(substring) == False]
            if i in library.active:
                strings = ['-'.join(j.split('-')[:-1])+"-" for j in compounds.ID]
                tmpDf = pd.read_sql_query('SELECT ID,"'+selected+'" FROM data WHERE ID LIKE "'+substring+'%"', conn)
                resultDf = resultDf.append(tmpDf)
        actives = [axis_map[color_by.value], axis_map[alpha_by.value], axis_map[size_by.value], 'ID', 'png', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
        for l in list(compounds):
            if l not in actives:
                compounds = compounds.drop(l, 1)
        compounds = pd.merge(compounds, resultDf, on=['ID'], sort=True)
        compounds = compounds.dropna()
    update()

def select_pca(attr, old, new):
    global compounds
    for p in ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']:
        compounds = compounds.drop(p,1)
    resultDf = pd.read_pickle(join(dirname(__file__), pca_map[pca_select.value]+".pkl"))
    compounds = pd.merge(compounds, resultDf, on=['ID'], sort=True)
    update()

def update():
    global compounds
    if len(compounds) > 0:
        cnorm = min_max_scaler.fit_transform(np.array(compounds[axis_map[color_by.value]]).reshape(-1, 1))
        cnorm = np.array(cnorm).T[0]
        colors = [
            "#%02x%02x%02x" % (int(r), int(g), int(b)) for r, g, b, _ in 255*plt.cm.coolwarm(cnorm)
        ]
        source.data = dict(
            size=10*min_max_scaler.fit_transform(np.array(compounds[axis_map[size_by.value]]).reshape(-1, 1))+1,
            size_real=compounds[axis_map[size_by.value]],
            alpha=(min_max_scaler.fit_transform(np.array(compounds[axis_map[alpha_by.value]]).reshape(-1, 1))+0.1)/1.1,
            alpha_real=compounds[axis_map[alpha_by.value]],
            color=colors,
            color_real=compounds[axis_map[color_by.value]],
            png=compounds["png"],
            PC1=compounds["PC1"],
            PC2=compounds["PC2"],
            PC3=compounds["PC3"],
            PC4=compounds["PC4"],
            PC5=compounds["PC5"]
        )

library.on_change('active', select_library)
color_by.on_change('value', select_color)
alpha_by.on_change('value', select_alpha)
size_by.on_change('value', select_size)
pca_select.on_change('value', select_pca)

# page layout
desc = Div(text=open(join(dirname(__file__), "description.html")).read(), width=1000)
footer = Div(text=open(join(dirname(__file__), "footer.html"), encoding="utf8").read(), width=1200)

sizing_mode = 'fixed'

library_box = widgetbox(pca_select, color_by, alpha_by, size_by, library, width=500)
main_row = row(library_box, column(p, footer))

l = layout([
    [desc],
    [main_row]
], sizing_mode=sizing_mode)

update()

curdoc().add_root(l)
curdoc().title = 'cfbsgC'
