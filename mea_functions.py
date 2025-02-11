# Modules for averageMetricsandDates
import numpy as np
import math
from datetime import datetime
import inspect
from collections import OrderedDict

# Modules for averageAmps
# import math

# Modules for statistics_chart
from copy import copy
from copy import deepcopy
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import statsmodels.api as sm
from bioinfokit.analys import stat
from statsmodels.formula.api import ols
import pingouin as pg
from statsmodels.graphics.factorplots import interaction_plot
# from brokenaxes import brokenaxes
import matplotlib.patheffects as pe
from matplotlib import container
import tkinter as tk
import tkinter.filedialog as fd

def averageMetricsAndDates(indexes,dataList,cleandatesdict,datesforaveragemerging,substring):
    # First block calculates average metric values for whole duration of argument 'substring' recordings (use for normalization)
    # Second block averages metrics for each hour (6x10-minute time points) of pretreatment, bl and ac
    # Third block creates a remapped dates dictionary of the affected time points
    if len(indexes) == 1:
        dataListWholeAvg = dataList[indexes[0]]
        dataListperHour = [dataList[indexes[0]]]
        defdates = OrderedDict()
        defdates[list(cleandatesdict.keys())[indexes[0]]] = list(cleandatesdict.values())[indexes[0]]
        print(f'{substring} has only one time point: no dataList entries were averaged')
    elif len(indexes) > 1:
        dataListWholeAvg = []
        dataListperHour = [[] for _ in range(math.ceil(len(indexes)/6))]
        datestoaverage = []
        metric_for_average = []
    # First block
        treatment_info_length = 4 if len(dataList[0][4]) <= 25 else 7
        for i in range(treatment_info_length, len(dataList[0])):
            if any(substring in e.strip() for e in cleandatesdict.keys()):
                # for all time points
                for timepoint, metric in enumerate(dataList[min(indexes):max(indexes)+1]):
                    temp = [float(_ if _ else 0) for _ in metric[i][1:]]
                    arr = np.array(temp, dtype=np.double)
                    metric_for_average.append(arr)
                dataListWholeAvg.append(np.mean(metric_for_average, axis=0).tolist())
                dataListWholeAvg[-1].insert(0, metric[i][0])
                metric_for_average = []
        # Second block 
                # Average the metrics in the input 'indexes' argument time points by segments of one hour
                for timepoint, metric in enumerate(dataList[min(indexes):max(indexes)+1]):
                    temp = [float(_ if _ else 0) for _ in metric[i][1:]]
                    arr = np.array(temp, dtype=np.double)
                    metric_for_average.append(arr)
                    if len(indexes) > 6 and (timepoint+1) >= 6 and ((timepoint+1) % 6 == 0 or (timepoint+1) == len(indexes) % 6):
                        dataListperHour[math.ceil((timepoint+1)/6-1)].append(np.mean(metric_for_average, axis=0).tolist())
                        dataListperHour[math.ceil((timepoint+1)/6-1)][-1].insert(0, metric[i][0])
                        metric_for_average = []        
                    elif len(indexes) < 6 and (timepoint+1) == len(indexes):
                        dataListperHour[0].append(np.mean(metric_for_average, axis=0).tolist())
                        dataListperHour[0][-1].insert(0, metric[i][0])
                        metric_for_average = []
        # Add the first four elements (well, coloring, treatment and concentration) from the first instance of 'substring' in dataList to the newly created 'averaged data' instances
        for infolist in reversed(dataList[min(indexes)][:treatment_info_length]):
            dataListWholeAvg.insert(0, infolist)
        for dat in dataListperHour:
            for infolist in reversed(dataList[min(indexes)][:treatment_info_length]):
                dat.insert(0, infolist)
        # Third block
        for n, ind in enumerate(indexes):
            datestoaverage.append(list(datesforaveragemerging.values())[ind])
        segments = [datestoaverage[i:i + 6] for i in range(0, len(datestoaverage), 6)]
        avgTimes = []
        for dates in segments:
            avgTimes.append(datetime.strftime(datetime.fromtimestamp(sum(map(datetime.timestamp,dates))/len(dates)),'%Y-%m-%d %H:%M:%S'))
        newdates = list(cleandatesdict.keys())[min(indexes):min(indexes)+len(dataListperHour)]
        defdates = OrderedDict()
        for i, _ in enumerate(newdates):
            defdates[_] = avgTimes[i]
        print(f'{substring} has {len(indexes)} time points: dataList entries were averaged')
    else:
        dataListWholeAvg = []
        dataListperHour = []
        defdates = {}
        def retrieve_name(var):
            callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
            return [var_name for var_name, var_val in callers_local_vars if var_val is var]
        print(f'Input argument {retrieve_name(indexes)} for {substring} is empty: cannot average time points')
    return dataListWholeAvg, dataListperHour, defdates

def data_and_date_reconstruction(dataList, messydatesdict, cleandatesdict, datesnottobemergeddict,
                                 dataListWholeAvgbl, dataListperHourbl, defdatesbl,
                                 dataListWholeAvgpre, dataListperHourpre, defdatespre,
                                 dataListWholeAvgac, dataListperHourac, defdatesac):
    
    if any((dataListWholeAvgbl, dataListperHourbl, defdatesbl,
       dataListWholeAvgpre, dataListperHourpre, defdatespre,
       dataListWholeAvgac, dataListperHourac, defdatesac)):
        
        # If necessary, remake name|date dictionary using the new averaged time points and
        # Merge the averaged results in the adequate order within the remaining non-averaged data
        deftimepointsdict = defdatesbl|defdatespre|defdatesac
        defcleandateslist = [None for _ in range(len(cleandatesdict))]

        specialtimepointsdataList = dataListperHourbl+dataListperHourpre+dataListperHourac
        averageddataList = [None for _ in range(len(dataList))]

        for k, v in deftimepointsdict.items():
            defcleandateslist[list(cleandatesdict.keys()).index(k)] = (k, deftimepointsdict[k])
            averageddataList[list(cleandatesdict.keys()).index(k)] = specialtimepointsdataList[list(deftimepointsdict.keys()).index(k)]
        for k, v in datesnottobemergeddict.items():
            defcleandateslist[list(cleandatesdict.keys()).index(k)] = (k, datesnottobemergeddict[k])
            averageddataList[list(cleandatesdict.keys()).index(k)] = dataList[list(cleandatesdict.keys()).index(k)]

        defcleandateslist = list(filter(None, defcleandateslist))
        defcleandatesdict = OrderedDict(defcleandateslist)
        averageddataList = list(filter(None, averageddataList))

        # Replace the original ones with the refurbished ones
        dataList = averageddataList
        cleandatesdict = defcleandatesdict

        print(f'Original dataList and cleandatesdict were replaced.\nThere are: {len(defcleandatesdict)} definitive time points, while there were {len(messydatesdict)} original time points\n', *defcleandatesdict,
              sep="\n")
        return dataList, cleandatesdict
    else:
        print('Neither bl, pre, nor ac had multiple time points')
        return dataList, cleandatesdict

def averageAmps(amp_day_uV, replicates_list, substring):
    # First block: average all bl, pre and/or ac treatment amplitudes
    amp_day_uV_merged = {}
    mergeddict = {w: [] for w in replicates_list}
    mydict = [k for k in amp_day_uV if substring.lower() in k.lower()]

    for el in mydict:
        for k, v in amp_day_uV[el].items():
            mergeddict[k].extend(v)
    try:
        amp_day_uV_merged[mydict[0]] = mergeddict
        print(f'All {substring} amplitudes merged for each well')
    except IndexError:
        print(f'No {substring} time points found')
    
    # Second block: piecewise average of every 60 minutes of recording
    finalmydict = mydict[:len([None for _ in range(math.ceil(len(mydict)/6))])]
    amp_day_uV_merged_per_hour = {k: {w:[] for w in replicates_list} for k in finalmydict}

    for i, el in enumerate(mydict):
        for k, v in amp_day_uV[el].items():
            amp_day_uV_merged_per_hour[finalmydict[math.ceil((i+1)/6-1)]][k].extend(v)       
    if amp_day_uV_merged_per_hour:
        print(f'All {substring} amplitudes merged per hour of recording')
    else:
        print(f'No {substring} time points found: could not perform per-hour calculations')
    return amp_day_uV_merged, amp_day_uV_merged_per_hour

# Function to use ideal x tick labels with ideal separation in time instead of the real ones retrieved from the csv files
def ideal_x_axis(df):
    def_labels = {}
    ac_loc = df['Number_of_Spikes'].index.get_loc('ac')
    bl_loc = df['Number_of_Spikes'].index.get_loc('bl')
    hours = -ac_loc * 24
    for dd in df['Number_of_Spikes'].index:
        if dd.startswith('div'):
            hours += 24
        elif dd.startswith('bl'):
            hours = -1
        elif dd.lower().startswith('pre'):
            hours = -0.5
        elif dd.startswith('ac'):
            hours = 0
        elif dd.startswith('on') or dd.startswith('aft'):
            hours += 1
        def_labels[dd] = hours
    def_labels.update(bl='Baseline -0.5', ac='Treatment 0')
    return def_labels
    
def statistics_chart(df):
    # Create a composition with all the charts using add_gridspec and add_subplot
    local_df = deepcopy(df) # Make a deepcopy of df to avoid modifying the original df

    # Ignore these warnings:
    warnings.filterwarnings(action='default',
                            message='(.*divide by zero.*|.*Mean of empty slice.*|.*Degrees of freedom.*|.*invalid value.*|.*Epsilon values.*)')

    # Create parameter data and SEM averages for each treatment directly from pandas dataframe content
    two_way_rm_anovas = {}
    ancovas = {}
    two_way_unb_anova = {}
    # three_way_anova = {}
    tukeys_test = {}
    levenes_test = {}
    parameters = dict()
    sem = dict()

    # with interactive plotting (refreshes itself and shows individual plots after every iteration of the loop)
    with plt.ion():
        for metric_name, df_metric in tqdm_notebook(local_df.items(), desc = 'Creating DataFrames for rm_ANOVA on ' + str(time.ctime(time.time()))):
        # Extract numbers and calculate Average and SEM for each treatment group and time point before the local_df DataFrame gets altered
            print('Calculating statistics for: ' + metric_name)
            parameters[metric_name] = {k: [] for k in res}
            sem[metric_name] = {k: [] for k in res}
            for day in df_metric.index:
                for trt in parameters[metric_name]:
                    vals_over_trt = np.array([df_metric.loc[day][key] for key in df_metric.keys() if key.startswith(trt)])
                    with np.errstate(invalid='ignore'):
                        parameters[metric_name][trt].append(np.nanmean(vals_over_trt))
                        sem[metric_name][trt].append(stats.sem(vals_over_trt, axis=0, ddof=1, nan_policy='omit'))

            # Create grid
            slide = plt.figure(constrained_layout=False, figsize=(25,35), dpi=100)
            #         gs1 = slide.add_gridspec(nrows=1, ncols=1, left=0, right=1, wspace=0.1, hspace=0.1) # For the time vs average parameter with SEM plot
            gs2 = slide.add_gridspec(nrows=6, ncols=3, left=0, right=1, top=1, bottom=0,
                                     height_ratios=[1,0.5,0.5,1,1,1], width_ratios=[1,1,1], wspace=0.2, hspace=0.3) # For the statistical test and data assessment plots
            #         slide_ax0 = slide.add_subplot(gs1[:, :])
            slide_ax1 = slide.add_subplot(gs2[0, :])
            slide_ax2 = slide.add_subplot(gs2[1, :])
            slide_ax3 = slide.add_subplot(gs2[2, :])
            slide_ax4 = slide.add_subplot(gs2[3, :])
            slide_ax5 = slide.add_subplot(gs2[4, 0])
            slide_ax5_1 = slide.add_subplot(gs2[4, 1:])
            slide_ax6 = slide.add_subplot(gs2[5, 0])
            slide_ax7 = slide.add_subplot(gs2[5, 1])
            slide_ax8 = slide.add_subplot(gs2[5, 2])

            allaxes = slide.get_axes()

            # Create a long-format table for statistical calculations
            cc = pd.DataFrame(local_df[metric_name])
            cc.columns.name = 'Treatment'
            cc.index.name = 'Time'
            cc = cc.rename({c: c + '.0' for c in local_df[metric_name].columns if c in res}, axis=1, inplace=False)
            cc.columns = metric_name + '_' + cc.columns
            cc['Time'] = cc.index
        #     cc = cc.drop(cc.columns[[cc.columns.get_loc(c) for c in cc.columns if 'Inactive' in c or 'Dead' in c or 'Untreated' in c]], axis=1, inplace=False)

            # Still some limitations in treatment name recognition through regex:
            # logic needs to be improved, as it fails to recognize some treatment names

            cc = pd.wide_to_long(cc, stubnames = metric_name, i = 'Time', j = 'Treatment.Replicate', sep='_',
                                 suffix='({}.\\d{{1,2}})'.format('.\\d{1,2}|'.join(list(res.keys())).replace('/', '\\'+'/').replace('+','\\'+'+').replace('%', '\\'+'%').replace('(', '\\'+'(').replace(')', '\\'+')')))
        #     cc.to_csv(path_or_buf=r'C:\Users\lab\Desktop\DataFrame.csv', encoding='UTF-8-sig') # export active DataFrame to a .csv file on the Desktop
            cc = cc.reset_index(inplace=False)
            cc[['Treatment','Replicate']] = cc['Treatment.Replicate'].str.split(".", n=1, expand=True)
            del cc['Treatment.Replicate']
        #     cc.columns = cc.columns.str.replace(" ", "_")
        #     cc.columns.str.replace("[_(μV)]", "",regex='False')
        #     cc = cc.replace(np.nan, 0, regex=True) # replace all NaN values with 0 in the long-format DataFrame
            # End of the DataFrame formatting

            # Perform statistical tests for each metric/worksheet within the excel

        #     cc_full = cc.replace(np.nan, 0, regex=True)
            try:
                two_way_rm_anovas[metric_name] = pg.rm_anova(cc, dv=metric_name, within=['Time', 'Treatment'], subject='Replicate',
                                                             detailed=True, correction=True, effsize='n2').round(6)
                ancovas[metric_name] = pg.rm_anova(cc, dv=metric_name, within=['Time', 'Treatment'], subject='Replicate',
                                                             detailed=True, correction=True, effsize='n2').round(6)
                two_way_unb_anova[metric_name] = pg.anova(cc, dv=metric_name, between=['Time', 'Treatment'],
                                                          detailed=True, effsize="n2").round(6)
    #             three_way_anova[metric_name] = pg.anova(cc, dv=metric_name, between=['Time', 'Treatment', 'Replicate'],
    #                                                     detailed=True, ss_type=3, effsize="n2").round(6) # We do not usually have the appropiate data for a three-way ANOVA
                tukeys_test[metric_name] = pg.pairwise_tukey(cc, dv=metric_name, between='Treatment', effsize='eta-square')
            except KeyError as e:
                print('Cannot perform two-way repeated measures Anova on {}: {} only contains NaN values'.format(metric_name.replace(' ', '_'), e))
                continue

            strip_plot = sns.stripplot(x='Treatment', y=metric_name, hue='Replicate', data=cc, ax=slide_ax1, size=7, edgecolor='black', linewidth=0.5, dodge=True) # Replicate means well
            box_plot = sns.boxplot(x='Treatment', y=metric_name, hue='Treatment', data=cc, ax=slide_ax1, dodge=False, palette = [v for v in res.values()])
        #     box_plot.legend(bbox_to_anchor=(1.03, 1), borderaxespad=0., labelspacing = 1, fontsize=17)
            slide_ax1.set_xlabel('Treatment', loc='center', labelpad=20, fontsize=18)
            slide_ax1.set_title(metric_name, pad=20, fontsize=20)
            slide_ax1.tick_params(axis='both', which='major', labelsize=17)
            slide_ax1.set_ylabel(metric_name, labelpad=25, fontsize=20)

            # Add repeated measures two-way ANOVA results table to the graph
            two_way_rm_anovas_text = []
            for row in range(len(two_way_rm_anovas[metric_name])):
                two_way_rm_anovas_text.append(two_way_rm_anovas[metric_name].iloc[row])
            table1 = slide_ax2.table(cellText=two_way_rm_anovas_text, colLabels=two_way_rm_anovas[metric_name].columns, loc='bottom', bbox=[0, -0.1, 1, 1])
            slide_ax2.set_title('Two-way repeated measures Anova test within time and treatment groups for {}'.format(metric_name),
                               fontsize=18, pad=0)
            slide_ax2.axis('off')
    #         table1.auto_set_font_size(False)
    #         table1.set_fontsize(12)

            # Add unbalanced two-way ANOVA table to the graph
            two_way_unb_anova_text = []
            for row in range(len(two_way_unb_anova[metric_name])):
                two_way_unb_anova_text.append(two_way_unb_anova[metric_name].iloc[row])
            table2 = slide_ax3.table(cellText=two_way_unb_anova_text, colLabels=two_way_unb_anova[metric_name].columns, loc='bottom', bbox=[0, 0, 1, 1])
            slide_ax3.set_title('Two-way unbalanced samples Anova test within time and treatment groups for {}'.format(metric_name),
                               fontsize=18, pad=20)
            slide_ax3.axis('off')
    #         table2.auto_set_font_size(False)
    #         table2.set_fontsize(11)

               # Add Tukey's test table to the graph
            tukeys_text = []
            for row in range(len(tukeys_test[metric_name])):
                tukeys_text.append(tukeys_test[metric_name].iloc[row])
            table3 = slide_ax4.table(cellText=tukeys_text, colLabels=tukeys_test[metric_name].columns, loc='bottom', bbox=[0, 0, 1, 1])
    #         table3.auto_set_font_size(False)
            slide_ax4.set_title('Pairwise Tukey\'s post-hoc test (table) within treatment groups for {}'.format(metric_name),
                               fontsize=18, pad=20)
            slide_ax4.axis('off')
    #         table3.set_fontsize(11)

            interaction_plot(x=cc['Replicate'], trace=cc['Treatment'], response=cc[metric_name],
                             ax=slide_ax5, plottype='b', colors=[v for v in res.values()])
            slide_ax5.set_title('Interaction plot treatments ~ replicates {}'.format(metric_name),
                               fontsize=18, pad=20)

    #         interaction_plot(x=cc['Time'], trace=cc['Treatment'], response=cc[metric_name],
    #                          ax=slide_ax5_1, plottype='b', colors=[v for v in res.values()])
            sns.lineplot(x='Time', y=metric_name, data=cc, hue='Treatment', palette=[v for v in res.values()],
                         estimator='mean', ax=slide_ax5_1, err_style='bars', n_boot=1000, sort=True, ci=95)
            slide_ax5_1.set_title('Interaction plot treatments ~ time {} 95% CI bars'.format(metric_name),
                                 fontsize=18, pad=20)
            plt.setp(slide_ax5_1.xaxis.get_majorticklabels(), rotation=-25, ha='left', rotation_mode='anchor')
    #         field = "Time"
    #         time_order = [list(df[metric_name].index)]
    #         cc.set_index(field).loc[time_order].plot

            residuals = stat()
            residuals.anova_stat(df=cc, res_var='Q(\'{}\')'.format(metric_name),
                                 anova_model='{}~C(Treatment)+C(Time)+C(Treatment):C(Time)'.format('Q(\'{}\')'.format(metric_name)))
            pg.qqplot(residuals.anova_std_residuals, dist='norm', confidence=False, ax = slide_ax6)
            pg.qqplot(residuals.anova_std_residuals, dist='expon', ax=slide_ax7)

            # histogram
            slide_ax8 = plt.hist(residuals.anova_model_out.resid, bins='auto', histtype='bar', ec='k', zorder=1) 
            plt.xlabel('Residuals')
            plt.ylabel('Frequency')
            plt.title('Histogram of residuals\' distribution')

            # Shapiro-Wilkins test for normal distribution of residuals 
            # (Null hypothesis: data were drawn from a normal distribution)
            w, pvalue = stats.shapiro(residuals.anova_model_out.resid)
            textstr1 = '\n'.join((r'$\bf{Shapiro-Wilkins test}$',
                                  r'$\mathrm{w}=%.30f$' % (w, ),
                                  r'$\mathrm{p-value}=%.30f$' % (pvalue, )))
            # these are matplotlib.patch.Patch properties for the text box
            hist_props = dict(boxstyle='round', facecolor='wheat', alpha=.5)
            # add text box with the data to the plot
            plt.text(0.555, 0.99, textstr1, transform=allaxes[8].transAxes, fontsize=9,
                    verticalalignment='top', bbox=hist_props)

            # Levene's test for homogeneity of variances
            residuals.levene(df=cc, res_var=metric_name, xfac_var=['Treatment', 'Time'])   
            lev_text = []
            levene_cell_colors = [["wheat","wheat"],["wheat","wheat"],[ "wheat","wheat"]]
            for row in range(len(residuals.levene_summary)):
                lev_text.append(residuals.levene_summary.iloc[row])
            table4 = plt.table(cellText=lev_text, colLabels=residuals.levene_summary.columns,
                              transform=allaxes[8].transAxes, loc='left', edges='closed',
                              cellColours=levene_cell_colors, colColours=["gainsboro","gainsboro"],
                               bbox=[0,0.75,0.25,0.25], zorder=2, alpha=.5)
            # Make table4 semitransparent
            for cell in table4._cells:
                table4._cells[cell].set_alpha(.5)
            plt.text(0.01, 1.035, 'Levene\'s test summary', transform=allaxes[8].transAxes,
                    verticalalignment='top', fontsize=9, weight='bold')

            # Save each figure within the loop
            if os.path.isdir(os.path.join(output_dir,'Descriptive statistics and plots')) == False:
                os.mkdir(os.path.join(output_dir,'Descriptive statistics and plots'))
            slide.savefig(os.path.join(output_dir,'Descriptive statistics and plots\\') + metric_name, edgecolor='none', transparent=False)

    #         plt.show(close=None, block=None) # Inline plotting consumes many resources and hugely increases running time
    #         plt.close(slide)

    plt.close("all") # Close all figures and do not display any of them inline after running the loop