#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
from collections import defaultdict
from scipy import stats
import time
import csv
import xlsxwriter
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def create_excel_table(fname):
    path =  os.getcwd()
    print(path)
    excelName = os.path.join(path, fname)

    os.chdir(path)

    # Couple each file with its creation date
    datesdict = {}
    for fl in os.listdir():
        if fl.split('.')[-1] == 'csv':
            with open(fl) as csv_file:
                csv_reader = csv.reader(csv_file)
                rows = list(csv_reader)
                datesdict[fl]= rows[6][1]

    # sort file names by their date of creation value
    from datetime import datetime
    nombres = []
    for n, d in sorted(datesdict.items(), key=lambda x: datetime.strptime(x[1], '%m/%d/%Y %H:%M:%S')):
        nombres.append(n)

    # Sort file names by date and append the desired rows to daticosList
    nombresList = []
    daticosList = []
    rowsChulas = np.array([68,72,70,125,126,127,141,150,159,167]) - 1
    for fname in nombres:
        # OPEN CSV
        if fname.split('.')[-1] == 'csv':
            # append file names to "nombresList"
            day = fname.split('.')[0].replace('(000)', '')
            nombresList.append(day)
            with open(fname) as csv_file:
                singleCSV = []
                csv_reader = csv.reader(csv_file, delimiter=',')
                line_count = 0
                for fname in csv_reader:
                    if line_count in rowsChulas:
                        singleCSV.append(fname)
                    line_count += 1
                daticosList.append(singleCSV)

    data = defaultdict(dict)

    for day, ds in zip(nombresList, daticosList):
        for row in ds:
            key = row[0]
            data[day][key] = row[1:]

    for day in nombresList:
        keys = list(data[day].keys())
        vals = list(data[day].values())
        idx = np.argsort(data[day]['Treatment'])

    ### append baseline values for normalization to norms dictionary
    norms = dict()
    for name, quantity in data['bl'].items():
        if name not in ['Well', 'Well Coloring', 'Treatment']:
            norms[name] = []
            for val in quantity:
                try:
                    norms[name].append(float(val))
                except ValueError:
                    norms[name].append(None)

    ### create excel workbook
    with xlsxwriter.Workbook(excelName) as workbook:
        styles = dict()
        data_0 = data[nombresList[0]]
        sheets_names = []
        for key in data_0.keys():
            if key not in ['Well Coloring', 'Well', 'Treatment']:
                sheets_names.append(key)
        for name in sheets_names:
            workbook.add_worksheet(name)
        for name in sheets_names:
            workbook.add_worksheet(name='n{}'.format(name))
        sheets = [sh for sh in workbook.worksheets() if sh.name in sheets_names]
        norm_sheets = [sh for sh in workbook.worksheets() if sh.name not in sheets_names]
        for i, day_data in enumerate(data.values()):
            for color, treatment in zip(day_data['Well Coloring'], day_data['Treatment']):
                if not treatment in styles:
                    styles[treatment] = workbook.add_format({'bg_color': color})
            for sh in sheets:
                sh.write(i + 2, 0, nombresList[i])
                for j in range(len(day_data['Well'])):
                    sh.write(0, j + 1, day_data['Well'][idx[j]])
                    treatment = day_data['Treatment'][idx[j]]
                    sh.write(1, j + 1, treatment)
                    raw = day_data[sh.name][idx[j]]
                    try:
                        val = float(raw)
                        if val == 0:
                            val = None
                    except ValueError as e:
                        #print(str(e).replace(':', ' in') + " {}, {}, {}".format(nombresList[i],
                                                            #sh.name,
                                                            #day_data['Well'][idx[j]]))
                        val = None
                    sh.write(i + 2, j + 1, val, styles[treatment])
            for sh in norm_sheets:
                sh.write(i + 2, 0, nombresList[i])
                for j in range(len(day_data['Well'])):
                    sh.write(0, j + 1, day_data['Well'][idx[j]])
                    treatment = day_data['Treatment'][idx[j]]
                    sh.write(1, j + 1, treatment)
                    raw = day_data[sh.name[1:]][idx[j]]
                    try:
                        val = float(raw) / norms[sh.name[1:]][idx[j]]
                        if val == 0:
                            val = None
                    except ValueError:
                        val = None
                    except TypeError:
                        val = None
                    sh.write(i + 2, j + 1, val, styles[treatment])

                    
def plot_metrics(fname):
    ### create pandas dataframe
    df = pd.read_excel(fname, index_col=0, header=1, sheet_name=None)
    print(df.keys())

    ### create metrics dictionary from pandas dataframe for plotting data
    metrics = dict()

    for metric_name, df_metric in df.items():
        metrics[metric_name] = dict(N=[], C=[], R=[])
        for day in df_metric.index:
            for trt in metrics[metric_name]:
                vals_over_trt = np.array([df_metric.loc[day][key] for key in df_metric.keys() if key.startswith(trt)])
                with np.errstate(invalid='ignore'):
                    metrics[metric_name][trt].append(np.mean(vals_over_trt, where=~np.isnan(vals_over_trt)))

    for name, mean_metr in metrics.items():
        if name.startswith('n'):
            continue
        devel = []
        day_of_trt = []
        for val in df[name].index:
            if val == 'bl':
                break
            devel.append(val)
        for val in df[name].index[len(devel):]:
            if val.startswith('div'):
                break
            day_of_trt.append(val)
        rest = df[name].index.values[len(devel) + len(day_of_trt):]
        fig, axs = plt.subplots(3, 1, figsize=(10, 15), dpi=100, gridspec_kw=dict(hspace=0.3))
        for ind, ax in zip([devel, day_of_trt, rest], axs):
            for trt, row in mean_metr.items():
                x = range(len(ind))
                ax.plot(x, np.array(row)[x], label=trt)
                ax.set_ylabel(name)
                ax.set_xticks(x)
                ax.set_xticklabels(ind,
                                   fontdict=dict(size=10, rotation=45))
                ax.legend(loc='best')
        axs[0].set_title('Days Before Treatment')
        axs[1].set_title('Overnight Treatment')
        axs[2].set_title('Days After Treatment')

        
if __name__ == '__main__':
    fname = 'Advanced metrics tables.xlsx'
    create_excel_table(fname)
    plot_metrics(fname)

