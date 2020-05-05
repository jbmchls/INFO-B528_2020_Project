# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:19:19 2020
@author: jlbro

Jessie Browne
SOIC - IUPUI
MS Bioinformatics
INFO-B528 Janga

"""

import pandas as pd
import numpy as np
import operator as op
from plotly.offline import plot
import plotly.graph_objects as go
from statistics import mean, median, mode, stdev

################### DATA PROCESSING #####################
def flt(value):
    _fc=[]
    for h in value: 
        _fc.append(float(h))
    return _fc

def get_ft(sample): 
    _ft = []
    for value in sample: 
        _ft.append(op.itemgetter(0)(value))
    return _ft

# import abundance data for PRISM cohort
microb = pd.read_csv('microb_abund.csv', header=1, index_col=0)

# separate datasets based on [f-calprotectin] (fc) level
# and diagnosis (dg): Crohn's (CD), ulcerative colitis (UC)
fc = microb.loc['Fecal.Calprotectin']
fc = flt(fc)
dg = microb.loc['Diagnosis']
mg = microb.loc['SRA_metagenome_name']
ft = list(microb.columns)
sample_list = list(zip(ft, mg, dg, fc))

uc_hifc = []
uc_lofc = []
cd_hifc = []
cd_lofc = []
control_hifc = []
control_lofc = []
nofc = []
for i,j,k,l in sample_list: 
    if k == 'UC': 
        if l > 99.9:
            uc_hifc.append([i,j,l])
        if l < 100: 
            uc_lofc.append([i,j,l])
    if k == 'CD':
        if l > 99.9: 
            cd_hifc.append([i,j,l])
        if l < 100: 
            cd_lofc.append([i,j,l])
    if k == 'Control':
        if l > 99.9: 
            control_hifc.append([i,j,l])
        if l < 100: 
            control_lofc.append([i,j,l])
    if np.isnan(l): 
        nofc.append([i,j,k])
            
#verify all 220 samples are accounted for
divided = [uc_hifc, uc_lofc, cd_hifc, cd_lofc, control_hifc, control_lofc, nofc]
counts = []
for each in divided: 
    counts.append(len(each))
sum(counts)    

# statistics for [f-c] in UC subgroups
hi_uc_fc = []
for value in uc_hifc: 
    hi_uc_fc.append(op.itemgetter(2)(value))
lo_uc_fc = []
for value in uc_lofc: 
    lo_uc_fc.append(op.itemgetter(2)(value))
hi_cd_fc = []
for value in cd_hifc: 
    hi_cd_fc.append(op.itemgetter(2)(value))
lo_cd_fc = []
for value in cd_lofc: 
    lo_cd_fc.append(op.itemgetter(2)(value))
lo_control_fc = []
for value in control_lofc: 
    lo_control_fc.append(op.itemgetter(2)(value))
hi_control_fc = []
for value in control_hifc: 
    hi_control_fc.append(op.itemgetter(2)(value))
all_control_fc = lo_control_fc + hi_control_fc
all_uc_fc = hi_uc_fc + lo_uc_fc
all_cd = hi_cd_fc + lo_cd_fc
all_hi = hi_uc_fc + hi_cd_fc + hi_control_fc
all_lo = lo_uc_fc + lo_cd_fc + lo_control_fc
every = all_hi + all_lo

def stats(fc): 
    f.write('\nmean: ' + str(round(mean(fc),2)) + '\n')
    f.write('median: ' + str(median(fc)) + '\n')
    f.write('std dev: ' + str(round(stdev(fc),2)) + '\n')
    f.write('min: ' + str(min(fc)) + '\n')
    f.write('max: ' + str(max(fc)) + '\n')
    f.write('n=' + str(len(fc)) + '\n\n')

f = open('fc_statistics.txt', 'w')
f.write('Statistics for f-calprotectin concentration in each group: \n\n')
f.close()

f = open('fc_statistics.txt', 'a')
f.write('High [f-c] CD: ')
stats(hi_cd_fc)
f.write('Low [f-c] CD: ')
stats(lo_cd_fc)
f.write('High [f-c] UC: ')
stats(hi_uc_fc)
f.write('Low [f-c] UC: ')
stats(lo_uc_fc)
f.write('High [f-c] Control: ')
stats(hi_control_fc)
f.write('Low [f-c] Control: ')
stats(lo_control_fc)
f.write('All CD: ')
stats(all_cd)
f.write('All UC: ')
stats(all_uc_fc)
f.write('All Control: ')
stats(all_control_fc)
f.write('All High: ')
stats(all_hi)
f.write('All Low: ')
stats(all_lo)
f.write('All: ')
stats(every)
f.close()

################# ANALYSIS #######################
# Microbioal species abundance in UC and Control
#create a list of # Feature / Sample for subgroups
uc_ft_hi = get_ft(uc_hifc)
uc_ft_lo = get_ft(uc_lofc)

#slice original tables using feature list
uc_hi_microb = microb[uc_ft_hi]
uc_lo_microb = microb[uc_ft_lo]

# Heatmap for UC grouped by [fc]
conc_fc_hi_uc = hi_uc_fc
conc_ft_hiuc = list(zip(uc_ft_hi,conc_fc_hi_uc))
conc_ft_hiuc = sorted(conc_ft_hiuc, key=op.itemgetter(1), reverse=True)
ft_hiuc = get_ft(conc_ft_hiuc)
hm_hiuc = uc_hi_microb[ft_hiuc]
hm_hiuc = hm_hiuc[8:209].astype(float)
hm_hiuc['total_RA'] = hm_hiuc.sum(axis=1)
hm_hiuc = hm_hiuc.sort_values('total_RA', ascending=False)
hm_z_hiuc = hm_hiuc.iloc[0:30,:26]
hm_z = hm_z_hiuc.to_numpy()
hm_y = hm_z_hiuc.index.to_list()
hm_x = sorted(conc_fc_hi_uc )

fig1 = go.Figure()

fig1 = go.Figure(
    data=go.Heatmap(
        z=hm_z,
        x=hm_x,
        y=hm_y,
        transpose=False,
        colorscale='Solar'),
    layout=go.Layout(
        xaxis=dict(type='category'),
        yaxis=dict(autorange='reversed'),
        xaxis_title="[f-c] micrograms/g"
        )
)

plot(fig1,show_link = True, filename = 'heatmap.html')

# narrow in on microbes uniquely enriched in UC
# sum RA across samples, rank, and divide by n
uc_hi_microb = uc_hi_microb[8:209].astype(float)
uc_hi_microb['total_RA'] = uc_hi_microb.sum(axis = 1)
uc_hi_microb = uc_hi_microb.sort_values('total_RA', ascending=False)
uc_hi_microb['avg_RA'] = uc_hi_microb['total_RA']/26

uc_lo_microb = uc_lo_microb[8:209].astype(float)
uc_lo_microb['total_RA'] = uc_lo_microb.sum(axis = 1)
uc_lo_microb = uc_lo_microb.sort_values('total_RA', ascending=False)
uc_lo_microb['avg_RA'] = uc_lo_microb['total_RA']/20

control_ft_hi = get_ft(control_hifc)
control_ft_lo = get_ft(control_lofc) 
control_hi_microb = microb[control_ft_hi]
control_lo_microb = microb[control_ft_lo]

control_hi_microb = control_hi_microb[8:209].astype(float)
control_hi_microb['total_RA'] = control_hi_microb.sum(axis = 1)
control_hi_microb = control_hi_microb.sort_values('total_RA', ascending=False)
control_hi_microb['avg_RA'] = control_hi_microb['total_RA']/6

control_lo_microb = control_lo_microb[8:209].astype(float)
control_lo_microb['total_RA'] = control_lo_microb.sum(axis = 1)
control_lo_microb = control_lo_microb.sort_values('total_RA', ascending=False)
control_lo_microb['avg_RA'] = control_lo_microb['total_RA']/36

#print to csv files for future analysis
uc_hi_microb.to_csv('PRISM_uc_hifc_microb_abund.csv')
uc_lo_microb.to_csv('PRISM_uc_lofc_microb_abund.csv')
control_hi_microb.to_csv('PRISM_control_hifc_microb_abund.csv')
control_lo_microb.to_csv('PRISM_control_lofc_microb_abund.csv')

# print the top RA microbes to file 
hi_uc_30 = uc_hi_microb.iloc[0:30,[27]]
lo_uc_30 = uc_lo_microb.iloc[0:30,[21]]
hi_uc_30.to_csv('PRISM_top_30_RA_microbes_inflammed_uc.txt', sep='\t')
lo_uc_30.to_csv('PRISM_top_30_RA_microbes_lofc_uc.txt', sep='\t')

# find microbes that are more abundant in UC than Control
def get_RA(ft, zipped):
    RA = []
    for i in ft: 
        for x,y in zipped:  
           if i == x:
                RA.append(y) 
    return RA

mc1 = hi_uc_30.index.to_list()
mc2 = lo_uc_30.index.to_list()       
xuc = list(dict.fromkeys(mc1+mc2))

mc = uc_hi_microb.index.to_list()
ra = list(uc_hi_microb['avg_RA'])
mc_ra = list(zip(mc,ra))
y_hiuc = get_RA(xuc, mc_ra)

mc = control_lo_microb.index.to_list()
ra = list(control_lo_microb['avg_RA'])
mc_ra = list(zip(mc,ra))
y_locontrol = get_RA(xuc, mc_ra)

control_vs_uc = list(zip(xuc,y_locontrol,y_hiuc))
suspects = []
for i,j,k in control_vs_uc: 
    if j < k: 
        suspects.append([i,k,j])

f = open('PRISM_UC_suspects.txt', 'w')
f.write('These 18 species were found to be more abundant in UC than in controls \n\n'
        'Species' + '\t' + 'UC (high [f-c])' + '\t' + 'Control (low [f-c])' + '\n')
f.close()
f = open('PRISM_UC_suspects.txt', 'a')
for i,j,k in suspects: 
        f.write(i + '\t' + str(round(j,5)) + '\t' + str(round(k,5)) + '\n')
f.close()

x_suspect = get_ft(suspects)
y_suspect_uc = []
y_suspect_control = []
for value in suspects: 
    y_suspect_uc.append(op.itemgetter(1)(value))
    y_suspect_control.append(op.itemgetter(2)(value))
    
fig2 = go.Figure()
fig2.add_trace(go.Bar(
    x=x_suspect,
    y=y_suspect_control,
    name='Control (low [f-c])',
    marker_color='cadetblue'
))
fig2.add_trace(go.Bar(
    x=x_suspect,
    y=y_suspect_uc,
    name='UC (high [f-c])',
    marker_color='indianred'
))
fig2.update_layout(
    xaxis_tickfont_size=14,
    yaxis=dict(
        title='Average Relative Abundance',
        titlefont_size=16,
        tickfont_size=14,
    ),
    legend=dict(
        x=0.75,
        y=1.0,
        bgcolor='rgba(255, 255, 255, 0)',
        bordercolor='rgba(255, 255, 255, 0)'
    ))
plot(fig2,show_link = True, filename = 'avg_RA_plot.html')

ratios = np.asarray(y_suspect_uc) / np.asarray(y_suspect_control)
species_ratios = list(zip(x_suspect, ratios.tolist()))
f = open('PRISM_UC_vs_control_ratios.txt', 'w')
f.write('Ratio of high [f-c] UC species RA by low [f-c] control \n\n'
        'Species' + '\t' + 'Ratio' + '\n')
f.close()
f = open('PRISM_UC_vs_control_ratios.txt', 'a')
for i,j in species_ratios: 
        f.write(i + '\t' + str(round(j,5)) + '\n')
f.close()

colors = ['maroon',] * 18
colors[2] = 'goldenrod'
colors[7] = 'goldenrod'
colors[11] = 'goldenrod'
colors[14] = 'goldenrod'

fig3 = go.Figure(go.Bar(
            x=ratios,
            y=x_suspect,
            orientation='h',
            marker_color=colors
))
fig3.update_layout(
    xaxis=dict(
        title='RA(UC) / RA(Control)',
        tick0=0,
        dtick=5)
    )

plot(fig3,show_link = True, filename = 'ratio_plot.html')





















