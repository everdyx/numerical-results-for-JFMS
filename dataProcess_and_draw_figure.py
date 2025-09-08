# -*- coding: utf-8 -*-
import pandas as pd
import os
import re
import numpy as np
import matplotlib.pyplot as plt


def plot_boxplot(series_list, title,
                 xlabel=None, ylabel=None, figsize=(10, 6),
                 linewidth=2.5, tick_labels=None, y_lim=None):

    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor('none')  
    ax.patch.set_facecolor('none')  


    if tick_labels is None:
        tick_labels = [s.name for s in series_list]


    bp = ax.boxplot(series_list,
                    patch_artist=True,  
                    tick_labels=tick_labels,
                    widths=0.6)  


    for box in bp['boxes']:
        box.set_facecolor('#87CEFA')
        box.set_edgecolor('#4682B4')
        box.set_linewidth(linewidth)


    for median in bp['medians']:
        median.set_color('#FF6347')
        median.set_linewidth(linewidth + 0.5)


    for whisker in bp['whiskers']:
        whisker.set_color('#4682B4')
        whisker.set_linestyle('-')
        whisker.set_linewidth(linewidth)

    for cap in bp['caps']:
        cap.set_color('#4682B4')
        cap.set_linewidth(linewidth)


    for flier in bp['fliers']:
        flier.set_color('#4682B4')
        flier.set_marker('o')
        flier.set_markersize(6)
        flier.set_markeredgewidth(linewidth - 0.5)


    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

  
    ax.tick_params(axis='x')
    ax.tick_params(axis='y')

    if y_lim is not None:
        ax.set_ylim(y_lim)  


    ax.yaxis.grid(True, linestyle='--', alpha=0.7, linewidth=linewidth - 1)

    plt.tight_layout(pad=0.1, h_pad=1, w_pad=0.1)
    return fig, ax


def get_file_names(folder_path):

    file_names = []
    try:

        if not os.path.exists(folder_path):
            raise FileNotFoundError(f"wrong：folder '{folder_path}' does not exist")


        if not os.path.isdir(folder_path):
            raise NotADirectoryError(f"wrong：'{folder_path}' is not a directory")


        items = os.listdir(folder_path)


        for item in items:
            item_path = os.path.join(folder_path, item)
            if os.path.isfile(item_path):
                file_names.append(item)

    except PermissionError:
        print(f"wrong：no permission '{folder_path}'")
    except Exception as e:
        print(f"wrong：unknown error - {e}")

    return file_names

#%%
df_exact={}
df_small1={}
df_small2={}
df_large={}
df_meta={}
scen = ["P2M2V2 CT3.0", "P2M2V2 CT8.0", "P2M2V2 CT13.0", "P3M2V2 CT3.0", "P3M2V2 CT8.0", "P3M2V2 CT13.0", "P10M2V2 CT3.0", "P10M2V2 CT8.0", "P10M2V2 CT13.0", "P20M2V2 CT3.0", "P20M2V2 CT8.0", "P20M2V2 CT13.0", "P10M8V2 CT3.0", "P10M8V8 CT3.0", "P10M8V14 CT3.0", "P10M8V20 CT3.0", "P20M8V2 CT3.0", "P20M8V8 CT3.0", "P20M8V14 CT3.0", "P20M8V20 CT3.0", ]

algo=["SD and CP", "JFMS\\small and medium","JFMS\\large", "LNS"]
current_dir = os.getcwd()
for i in range(len(algo)):
    folder_path=os.path.join(current_dir, algo[i])
    data_files=get_file_names(folder_path)
    pattern = r'^P(?P<P>\d+)M(?P<M>\d+)V(?P<V>\d+) CT(?P<CT>\d+\.\d+)(?=_|\.|$)'
    extracted_data = []
    for file in data_files:
        match = re.match(pattern, file)
        if match:
            p = int(match.group('P'))
            m = int(match.group('M'))
            v = int(match.group('V'))
            ct = float(match.group('CT'))
            prefix = f'P{p}M{m}V{v} CT{ct}'
            extracted_data.append((p, m, v, ct, prefix, file))


    sorted_data = sorted(extracted_data, key=lambda x: (x[0], x[1], x[2], x[3]))

    for data in sorted_data:
        path=os.path.join(folder_path, data[5])
        scenario=data[4]

        print(scenario)
        print(algo[i])
        df = pd.read_csv(path, sep=',')
        if algo[i]== "SD and CP":
            if len(df.columns)>10:
                df = df.drop(df.columns[-1], axis=1)
            rename_dict = {
                'SDgap': 'SD_Gap',
                'E': 'CP',
                'Egap': 'CP_Gap',
                'time_E':'time_CP',
                'prodSeq':'CP_prodSeq',
                'logiSeq':'CP_logiSeq'
            }
            df = df.rename(columns=rename_dict)
            df = df.drop('Gap-SD-E', axis=1)
            df = df.iloc[np.repeat(np.arange(len(df)), 5)].reset_index(drop=True)
            df_exact[scenario]=df
        elif algo[i]== "JFMS\\small and medium":
            if scenario in ["P2M2V2 CT3.0","P2M2V2 CT8.0","P2M2V2 CT13.0","P3M2V2 CT3.0","P3M2V2 CT8.0","P3M2V2 CT13.0"]:
                if len(df.columns)>20:
                    df = df.drop(df.columns[-1], axis=1)
                rename_dict = {
                    'UB':'JFMS',
                    'time_H':'time_JFMS',
                    'FCFF':'RULE',
                    'prodSeq':'JFMS_prodSeq',
                    'logiSeq':'JFMS_logiSeq',
                    'UBNoLS':'UBNL',
                    'LB':'LBL',
                    'LBNoLift':'LBNL',
                }
                df = df.rename(columns=rename_dict)
                df = df.drop(['SD','E','time_SD','time_E','improve','G(UB-OPT)','G(UB-LB)','LBImprov','UBImprov','G(FCFS-UB)'], axis=1)
                df['UBNL'] = df['UBNL'].where(df['UBNL'] <= df['RULE'], df['RULE'])
                df_small1[scenario] = df
            else:
                if len(df.columns)>14:
                    df = df.drop(df.columns[-1], axis=1)
                rename_dict = {
                    'UB': 'JFMS',
                    'time_H': 'time_JFMS',
                    'FCFF': 'RULE',
                    'prodSeq': 'JFMS_prodSeq',
                    'logiSeq': 'JFMS_logiSeq',
                    'UBNoLS': 'UBNL',
                    'LB': 'LBL',
                    'LBNoLift': 'LBNL',
                }
                df = df.rename(columns=rename_dict)
                df = df.drop([ 'G(UB-LB)', 'LBImprov', 'UBImprov','G(FCFS-UB)'], axis=1)
                df['UBNL'] = df['UBNL'].where(df['UBNL'] <= df['RULE'], df['RULE'])
                df_small2[scenario]=df
        elif algo[i]== "JFMS\\large":
            if len(df.columns)>14:
                df = df.drop(df.columns[-1], axis=1)
            rename_dict = {
                'UB': 'JFMS',
                'time': 'time_JFMS',
                'FCFS': 'RULE',
                'prodSeq': 'JFMS_prodSeq',
                'logiSeq': 'JFMS_logiSeq',
                'UBNoLS': 'UBNL',
                'LB': 'LBL',
                'LBNoLift': 'LBNL',
            }
            df = df.rename(columns=rename_dict)
            df = df.drop(['G(UB-LB)', 'LBImprov', 'UBImprov', 'G(FCFS-UB)'],axis=1)
            df['UBNL'] = df['UBNL'].where(df['UBNL'] <= df['RULE'], df['RULE'])
            df_large[scenario]=df
        elif algo[i] == "LNS":

            if len(df.columns)>6:
                df = df.drop(df.columns[-1], axis=1)
            rename_dict = {
                'time': 'time_LNS',
                'z_FCFS': 'RULE',
            }
            df = df.rename(columns=rename_dict)
            df = df.drop(['Gap-FCFS-LNS'], axis=1)
            df_meta[scenario] = df
        else:
            pass

#%%
df={}
for i in scen:
    if i in [ "P2M2V2 CT3.0","P2M2V2 CT8.0","P2M2V2 CT13.0","P3M2V2 CT3.0","P3M2V2 CT8.0","P3M2V2 CT13.0"]:
        df[i]=pd.concat([df_exact[i], df_small1[i],df_meta[i]], axis=1)
    elif i in ["P10M2V2 CT3.0","P10M2V2 CT8.0","P10M2V2 CT13.0","P20M2V2 CT3.0","P20M2V2 CT8.0","P20M2V2 CT13.0"]:
        df[i]=pd.concat([df_exact[i], df_small2[i],df_meta[i]], axis=1)
    else:
        df[i]=pd.concat([df_large[i],df_meta[i]], axis=1)
    df[i] = df[i].loc[:, ~df[i].columns.duplicated()]

all_columns=set(df[scen[0]].columns).union(set(df[scen[6]].columns)).union(set(df[scen[12]].columns))

for i in scen:
    for j in all_columns:
        if j not in df[i].columns:
            df[i][j]=np.nan
    df[i]['Best/Opt']=df[i][['SD','CP','JFMS','LNS','RULE']].min(axis=1)
    df[i]['G^{SD-CP}_{CP}'] = 100 * (df[i]['SD'] - df[i]['CP']) / df[i]['CP']
    df[i]['G^{JFMS-CP}_{CP}'] = 100 * (df[i]['JFMS'] - df[i]['CP']) / df[i]['CP']
    df[i]['G^{JFMS-LBL}_{JFMS}']=100*(df[i]['JFMS']-df[i]['LBL'])/df[i]['JFMS']
    df[i]['G^{RULE-JFMS}_{JFMS}'] = 100 * (df[i]['RULE'] - df[i]['JFMS']) / df[i]['JFMS']
    df[i]['G^{LBL-LBNL}_{LBNL}'] = 100 * (df[i]['LBL'] - df[i]['LBNL']) / df[i]['LBNL']
    df[i]['G^{UBNL-JFMS}_{JFMS}'] = 100 * (df[i]['UBNL'] - df[i]['JFMS']) / df[i]['JFMS']
    df[i]['G^{LNS-JFMS}_{JFMS}'] = 100 * (df[i]['LNS'] - df[i]['JFMS']) / df[i]['JFMS']

    df[i]['G^{SD-Best/Opt}_{SD}'] = 100 * (df[i]['SD'] - df[i]['Best/Opt']) / df[i]['SD']
    df[i]['G^{CP-Best/Opt}_{CP}'] = 100 * (df[i]['CP'] - df[i]['Best/Opt']) / df[i]['CP']
    df[i]['G^{JFMS-Best/Opt}_{JFMS}'] = 100 * (df[i]['JFMS'] - df[i]['Best/Opt']) / df[i]['JFMS']
    df[i]['G^{LNS-Best/Opt}_{LNS}'] = 100 * (df[i]['LNS'] - df[i]['Best/Opt']) / df[i]['LNS']
    df[i]['G^{RULE-Best/Opt}_{RULE}'] = 100 * (df[i]['RULE'] - df[i]['Best/Opt']) / df[i]['RULE']


    df[i]['LNS'] = df[i][['LNS', 'RULE']].min(axis=1)
    df[i]['UBNL'] = df[i][['UBNL', 'RULE']].min(axis=1)
    df[i] = df[i].loc[:, ~df[i].columns.duplicated()]

df_agg={}
for i in scen:
    df_agg[i]=df[i].groupby('instance').agg(
        SD=('SD', 'mean'),
        SD_Gap=('SD_Gap', 'mean'),
        CP=('CP', 'mean'),
        CP_Gap=('CP_Gap', 'mean'),
        time_SD_mean=('time_SD', 'mean'),
        time_CP_mean=('time_CP', 'mean'),
        LBL_best=('LBL', 'max'),
        LBL_mean=('LBL', 'mean'),
        JFMS_best=('JFMS', 'min'),
        JFMS_mean=('JFMS', 'mean'),
        LBNL_best=('LBNL', 'max'),
        LBNL_mean=('LBNL', 'mean'),
        UBNL_best=('UBNL', 'min'),
        UBNL_mean=('UBNL', 'mean'),
        RULE_best=('RULE', 'min'),
        RULE_mean=('RULE', 'mean'),
        LNS_best=('LNS', 'min'),
        LNS_mean=('LNS', 'mean'),
        time_JFMS_mean=('time_JFMS', 'mean'),
        time_LNS_mean=('time_LNS', 'mean'),
    ).reset_index()
    df_agg[i]['Best/Opt'] = df_agg[i][['SD', 'CP', 'JFMS_best', 'LNS_best', 'RULE_best']].min(axis=1)
    df_agg[i]['G^{SD-CP}_{SD}'] = 100 * (df_agg[i]['SD'] - df_agg[i]['CP']) / df_agg[i]['SD']
    df_agg[i]['G^{JFMS-CP}_{JFMS}'] = 100 * (df_agg[i]['JFMS_mean'] - df_agg[i]['CP']) / df_agg[i]['JFMS_mean']
    df_agg[i]['G^{JFMS-CP}_{JFMS}_best'] = 100 * (df_agg[i]['JFMS_best'] - df_agg[i]['CP']) / df_agg[i]['JFMS_best']
    df_agg[i]['G^{JFMS-LBL}_{JFMS}'] = 100 * (df_agg[i]['JFMS_mean'] - df_agg[i]['LBL_mean']) / df_agg[i]['JFMS_mean']
    df_agg[i]['G^{JFMS-LBL}_{JFMS}_best'] = 100 * (df_agg[i]['JFMS_best'] - df_agg[i]['LBL_best']) / df_agg[i]['JFMS_best']
    df_agg[i]['G^{JFMS-RULE}_{JFMS}'] = 100 * (df_agg[i]['JFMS_mean'] - df_agg[i]['RULE_mean']) / df_agg[i]['JFMS_mean']
    df_agg[i]['G^{LBNL-LBL}_{LBNL}'] = 100 * (df_agg[i]['LBNL_mean'] - df_agg[i]['LBL_mean']) / df_agg[i]['LBNL_mean']
    df_agg[i]['G^{JFMS-UBNL}_{JFMS}'] = 100 * (df_agg[i]['JFMS_mean'] - df_agg[i]['UBNL_mean']) / df_agg[i]['JFMS_mean']
    df_agg[i]['G^{JFMS-LNS}_{JFMS}'] = 100 * (df_agg[i]['JFMS_mean'] - df_agg[i]['LNS_mean']) / df_agg[i]['JFMS_mean']

    df_agg[i]['G^{SD-Opt}_{SD}'] = 100 * (df_agg[i]['SD'] - df_agg[i]['Best/Opt']) / df_agg[i]['SD']
    df_agg[i]['G^{CP-Opt}_{CP}'] = 100 * (df_agg[i]['CP'] - df_agg[i]['Best/Opt']) / df_agg[i]['CP']
    df_agg[i]['G^{JFMS-Opt}_{JFMS}'] = 100 * (df_agg[i]['JFMS_mean'] - df_agg[i]['Best/Opt']) / df_agg[i]['JFMS_mean']
    df_agg[i]['G^{LNS-Opt}_{LNS}'] = 100 * (df_agg[i]['LNS_mean'] - df_agg[i]['Best/Opt']) / df_agg[i]['LNS_mean']
    df_agg[i]['G^{RULE-Opt}_{RULE}'] = 100 * (df_agg[i]['RULE_mean'] - df_agg[i]['Best/Opt']) / df_agg[i]['RULE_mean']

df_all=pd.DataFrame(columns=['scen_id','num_SD','num_CP',
                             'SD_Gap','CP_Gap','num_JFMS',
                             'num_LNS','num_RULE','SD_Time',
                             'CP_Time','JFMS_Time','LNS_Time',
                             'G^{SD-CP}_{SD}','G^{JFMS-CP}_{JFMS}',
'G^{JFMS-CP}_{JFMS}_best','G^{JFMS-LBL}_{JFMS}','G^{JFMS-LBL}_{JFMS}_best',
'G^{JFMS-RULE}_{JFMS}','G^{LBNL-LBL}_{LBNL}','G^{JFMS-UBNL}_{JFMS}',
'G^{JFMS-LNS}_{JFMS}','G^{SD-Opt}_{SD}','G^{CP-Opt}_{CP}','G^{JFMS-Opt}_{JFMS}','G^{LNS-Opt}_{LNS}','G^{RULE-Opt}_{RULE}'])
err=0.3
for i in range(len(scen)):
    num_SD = (abs(df[scen[i]]['SD'] - df[scen[i]]['Best/Opt']) <= err).sum()/5
    num_CP = (abs(df[scen[i]]['CP'] - df[scen[i]]['Best/Opt']) <= err).sum()/5
    SD_Gap = df[scen[i]]['SD_Gap'].mean()*100
    CP_Gap = df[scen[i]]['CP_Gap'].mean()*100
    # num_JFMS=(abs(df[scen[i]]['JFMS'] - df[scen[i]]['Best/Opt'])<=err).sum()/5
    # num_LNS =(abs(df[scen[i]]['LNS'] - df[scen[i]]['Best/Opt']) <= err).sum()/5
    # num_RULE =(abs(df[scen[i]]['RULE'] - df[scen[i]]['Best/Opt']) <= err).sum()/5
    num_JFMS=(abs(df_agg[scen[i]]['JFMS_best'] - df_agg[scen[i]]['Best/Opt'])<=err).sum()
    num_LNS =(abs(df_agg[scen[i]]['LNS_best'] - df_agg[scen[i]]['Best/Opt']) <= err).sum()
    num_RULE =(abs(df_agg[scen[i]]['RULE_best'] - df_agg[scen[i]]['Best/Opt']) <= err).sum()
    SD_time=df[scen[i]]['time_SD'].mean()
    CP_time=df[scen[i]]['time_CP'].mean()
    JFMS_time=df[scen[i]]['time_JFMS'].mean()
    LNS_time=df[scen[i]]['time_LNS'].mean()
    a1=df_agg[scen[i]]['G^{SD-CP}_{SD}'].mean()
    a2=df_agg[scen[i]]['G^{JFMS-CP}_{JFMS}'].mean()
    a3=df_agg[scen[i]]['G^{JFMS-CP}_{JFMS}_best'].mean()
    a4=df_agg[scen[i]]['G^{JFMS-LBL}_{JFMS}'].mean()
    a5=df_agg[scen[i]]['G^{JFMS-LBL}_{JFMS}_best'].mean()
    a6=df_agg[scen[i]]['G^{JFMS-RULE}_{JFMS}'].mean()
    a7=df_agg[scen[i]]['G^{LBNL-LBL}_{LBNL}'].mean()
    a8=df_agg[scen[i]]['G^{JFMS-UBNL}_{JFMS}'].mean()
    a9=df_agg[scen[i]]['G^{JFMS-LNS}_{JFMS}'].mean()

    a10=df_agg[scen[i]]['G^{SD-Opt}_{SD}'].mean()
    a11 = df_agg[scen[i]]['G^{CP-Opt}_{CP}'].mean()
    a12 = df_agg[scen[i]]['G^{JFMS-Opt}_{JFMS}'].mean()
    a13 = df_agg[scen[i]]['G^{LNS-Opt}_{LNS}'].mean()
    a14 = df_agg[scen[i]]['G^{RULE-Opt}_{RULE}'].mean()

    df_all.loc[i]=[i,num_SD,num_CP,SD_Gap,CP_Gap,num_JFMS,num_LNS,num_RULE,SD_time,CP_time,JFMS_time,LNS_time,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14]

df_all['scen_id']=df_all['scen_id'].astype(int)
df_all['num_SD']=df_all['num_SD'].astype(int)
df_all['num_CP']=df_all['num_CP'].astype(int)
df_all['SD_Gap']=df_all['SD_Gap'].astype(float).round(2)
df_all['CP_Gap']=df_all['CP_Gap'].astype(float).round(2)
df_all['num_JFMS']=df_all['num_JFMS'].astype(int)
df_all['num_LNS']=df_all['num_LNS'].astype(int)
df_all['num_RULE']=df_all['num_RULE'].astype(int)
df_all['SD_Time']=df_all['SD_Time'].astype(float).round(2)
df_all['CP_Time']=df_all['CP_Time'].astype(float).round(2)
df_all['JFMS_Time']=df_all['JFMS_Time'].astype(float).round(2)
df_all['LNS_Time']=df_all['LNS_Time'].astype(float).round(2)
df_all['G^{SD-CP}_{SD}']=df_all['G^{SD-CP}_{SD}'].astype(float).round(2)
df_all['G^{JFMS-CP}_{JFMS}']=df_all['G^{JFMS-CP}_{JFMS}'].astype(float).round(2)
df_all['G^{JFMS-CP}_{JFMS}_best']=df_all['G^{JFMS-CP}_{JFMS}_best'].astype(float).round(2)
df_all['G^{JFMS-LBL}_{JFMS}']=df_all['G^{JFMS-LBL}_{JFMS}'].astype(float).round(2)
df_all['G^{JFMS-LBL}_{JFMS}_best']=df_all['G^{JFMS-LBL}_{JFMS}_best'].astype(float).round(2)
df_all['G^{JFMS-RULE}_{JFMS}']=df_all['G^{JFMS-RULE}_{JFMS}'].astype(float).round(2)
df_all['G^{LBNL-LBL}_{LBNL}']=df_all['G^{LBNL-LBL}_{LBNL}'].astype(float).round(2)
df_all['G^{JFMS-UBNL}_{JFMS}']=df_all['G^{JFMS-UBNL}_{JFMS}'].astype(float).round(2)
df_all['G^{JFMS-LNS}_{JFMS}']=df_all['G^{JFMS-LNS}_{JFMS}'].astype(float).round(2)

df_all['G^{SD-Opt}_{SD}']=df_all['G^{SD-Opt}_{SD}'].astype(float).round(2)
df_all['G^{CP-Opt}_{CP}']=df_all['G^{CP-Opt}_{CP}'].astype(float).round(2)
df_all['G^{JFMS-Opt}_{JFMS}']=df_all['G^{JFMS-Opt}_{JFMS}'].astype(float).round(2)
df_all['G^{LNS-Opt}_{LNS}']=df_all['G^{LNS-Opt}_{LNS}'].astype(float).round(2)
df_all['G^{RULE-Opt}_{RULE}']=df_all['G^{RULE-Opt}_{RULE}'].astype(float).round(2)

#%%
fontsize=36
plt.rcParams.update({

    "text.usetex": True,

    "font.family": "serif",
    "font.serif": ["Times New Roman", "Computer Modern Roman"],
    "text.latex.preamble": r"\usepackage{amsmath}",

    "font.size": fontsize,  
    "axes.labelsize": fontsize,  
    "axes.titlesize": fontsize,  
    "xtick.labelsize": fontsize,  
    "ytick.labelsize": fontsize,  
    "legend.fontsize": fontsize,  
    "figure.titlesize": fontsize  

})

series_list=[]
for i in range(12):
    series_list.append(df_agg[scen[i]]['G^{SD-CP}_{SD}'])
fig, ax=plot_boxplot(series_list, title=None, xlabel='Scenario', ylabel=r'$G^{\text{SD}}_{\text{CP}}$', figsize=(10, 6), linewidth=2.5, tick_labels=range(1,13),y_lim=[-50,50])
plt.savefig('G-SD-CP')
plt.show()


series_list=[]
for i in range(12):
    series_list.append(df_agg[scen[i]]['G^{JFMS-CP}_{JFMS}'])
fig, ax=plot_boxplot(series_list, title=None, xlabel='Scenario', ylabel=r'mean $G^{\text{JFMS}}_{\text{CP}}$', figsize=(10, 6), linewidth=2.5, tick_labels=range(1,13),y_lim=[-150,15])
plt.savefig('G-JFMS-CP')
plt.show()



series_list=[]
for i in range(12):
    series_list.append(df_agg[scen[i]]['G^{JFMS-LBL}_{JFMS}'])
fig, ax=plot_boxplot(series_list, title=None, xlabel='Scenario', ylabel=r'mean $G^{\text{JFMS}}_{\text{LBL}}$', figsize=(10, 6), linewidth=2.5, tick_labels=range(1,13),y_lim=[0,15])
plt.savefig('G-JFMS-LBL')
plt.show()



series_list=[]
for i in range(20):
    series_list.append(df_agg[scen[i]]['G^{JFMS-RULE}_{JFMS}'])
fig, ax=plot_boxplot(series_list, title=None, xlabel='Scenario', ylabel=r'mean $G^{\text{JFMS}}_{\text{RULE}}$', figsize=(16, 6), linewidth=2.5, tick_labels=range(1,21),y_lim=[-40,2])
plt.savefig('G-JFMS-RULE')
plt.show()


series_list=[]
for i in range(12):
    series_list.append(df_agg[scen[i]]['G^{LBNL-LBL}_{LBNL}'])
fig, ax=plot_boxplot(series_list, title=None, xlabel='Scenario', ylabel=r'mean $G^{\text{LBNL}}_{\text{LBL}}$', figsize=(10, 6), linewidth=2.5, tick_labels=range(1,13),y_lim=[-40,2])
plt.savefig('G-LBNL-LBL')
plt.show()

series_list=[]
for i in range(20):
    series_list.append(df_agg[scen[i]]['G^{JFMS-UBNL}_{JFMS}'])
fig, ax=plot_boxplot(series_list, title=None, xlabel='Scenario', ylabel=r'mean $G^{\text{JFMS}}_{\text{UBNL}}$', figsize=(16, 6), linewidth=2.5, tick_labels=range(1,21),y_lim=[-30,2])
plt.savefig('G-JFMS-UBNL')
plt.show()


series_list=[]
for i in range(20):
    series_list.append(df_agg[scen[i]]['G^{JFMS-LNS}_{JFMS}'])
fig, ax=plot_boxplot(series_list, title=None, xlabel='Scenario', ylabel=r'mean $G^{\text{JFMS}}_{\text{LNS}}$', figsize=(16, 6), linewidth=2.5, tick_labels=range(1,21),y_lim=[-20,5])
plt.savefig('G-JFMS-LNS')
plt.show()


df_all['scen_id']=df_all['scen_id']+1
df_all.to_excel(
    'results for table 4.xlsx',
    sheet_name='Sheet1',  
    index=False,  
    engine='openpyxl'  
)

#%%
N1=[]
N2=[]
N3=[]
N4=[]
for i in range(12):
    N1.append(int(( (df_exact[scen[i]]['time_CP']<900)).sum()/5))
    N2.append(int((df_exact[scen[i]]['SD'] >df_exact[scen[i]]['CP']+0.1).sum()/5))
    N3.append(int((df_exact[scen[i]]['SD'] < df_exact[scen[i]]['CP']-0.1).sum()/5))
    N4.append(int((df_exact[scen[i]]['SD'] == df_exact[scen[i]]['CP']).sum()/5))


fig, ax = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor('none')  
ax.patch.set_facecolor('none')  


x = range(1, len(N1) + 1)


line1, = ax.plot(x, N2, label=r'$G^{SD}_{CP}>0$', color='blue', marker='o',
                 linestyle='-', linewidth=4, markersize=12)
line3, = ax.plot(x, N3, label=r'$G^{SD}_{CP}<0$', color='red', marker='^',
                 linestyle='-', linewidth=4, markersize=12)


for i, (xi, n2, n3) in enumerate(zip(x, N2, N3)):

    ax.text(xi, n2 + 2,  
            f'{n2}',  
            fontsize=fontsize - 6,  
            ha='center')  
    if n3==21:
        ax.text(xi, n3 - 10.0,  
                f'{n3}',
                fontsize=fontsize - 6,
                ha='center')
    else:
        ax.text(xi, n3 + 2.0,  
                f'{n3}',
                fontsize=fontsize - 6,
                ha='center')


ax.set_xlabel('Scenario', fontsize=fontsize)
ax.set_ylabel('Number of instances', fontsize=fontsize)


ax.set_xticks(x)
ax.set_yticks([0, 25, 50, 75, 100])


ax.grid(True, linestyle='--', alpha=0.7)


ax.legend(fontsize=fontsize-6)


y_min, y_max = ax.get_ylim()
ax.set_ylim(y_min, y_max * 1.1)  

plt.tight_layout(pad=0.1, h_pad=1, w_pad=0.1)
plt.savefig('num ins')
plt.show()

#%%

fontsize = 36  
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Computer Modern Roman"],
    "text.latex.preamble": r"\usepackage{amsmath}",
    "font.size": fontsize,
    "axes.labelsize": fontsize,
    "axes.titlesize": fontsize,
    "xtick.labelsize": fontsize,
    "ytick.labelsize": fontsize,
    "legend.fontsize": fontsize ,  
    "figure.titlesize": fontsize
})


n_scenarios = 20  
lns_series_list = []  
rule_series_list = []  


for i in range(n_scenarios):
    scenario = scen[i]
    lns_series = df_agg[scenario]['G^{JFMS-LNS}_{JFMS}']
    lns_series.name = f'Scenario {i + 1} (LNS)'
    lns_series_list.append(lns_series)


    rule_series = df_agg[scenario]['G^{JFMS-RULE}_{JFMS}']
    rule_series.name = f'Scenario {i + 1} (RULE)'
    rule_series_list.append(rule_series)


combined_series = []
for l_ser, r_ser in zip(lns_series_list, rule_series_list):  
    combined_series.append(l_ser)
    combined_series.append(r_ser)


fig, ax = plt.subplots(figsize=(20, 6))  
fig.patch.set_facecolor('none')  
ax.patch.set_facecolor('none')  


positions = []
for i in range(n_scenarios):
    positions.append(i * 2.2 + 1)  
    positions.append(i * 2.2 + 1.8)  


bp = ax.boxplot(
    combined_series,
    positions=positions,
    patch_artist=True,
    widths=0.6,  
    showfliers=True
)


color_rule = {'box': '#FFA07A', 'edge': '#DC143C', 'median': '#006400'}  
color_lns = {'box': '#87CEFA', 'edge': '#4682B4', 'median': '#FF6347'}  


for i, (box, median, whisker, cap, flier) in enumerate(zip(
        bp['boxes'], bp['medians'], bp['whiskers'], bp['caps'], bp['fliers']
)):
    if i % 2 == 0:  
        box.set_facecolor(color_lns['box'])
        box.set_edgecolor(color_lns['edge'])
        box.set_linewidth(2.5)
        median.set_color(color_lns['median'])
        median.set_linewidth(3.0)
        whisker.set_color(color_lns['edge'])
        whisker.set_linestyle('-')
        whisker.set_linewidth(2.5)
        cap.set_color(color_lns['edge'])
        cap.set_linewidth(2.5)
        flier.set_color(color_lns['edge'])
        flier.set_marker('s')  
        flier.set_markersize(6)
        flier.set_markeredgewidth(2.0)
    else:  #
        box.set_facecolor(color_rule['box'])
        box.set_edgecolor(color_rule['edge'])
        box.set_linewidth(2.5)
        median.set_color(color_rule['median'])
        median.set_linewidth(3.0)
        whisker.set_color(color_rule['edge'])
        whisker.set_linestyle('-')
        whisker.set_linewidth(2.5)
        cap.set_color(color_rule['edge'])
        cap.set_linewidth(2.5)
        flier.set_color(color_rule['edge'])
        flier.set_marker('o')  
        flier.set_markersize(6)
        flier.set_markeredgewidth(2.0)


ax.set_ylim(-45, 3)
ax.set_yticks([0, -10, -20, -30, -40])
ax.set_ylabel(r'Gap (\%)', fontsize=fontsize)  


x_tick_positions = [i * 2.2 + 1.4 for i in range(n_scenarios)]  
x_tick_labels = [f'{i + 1}' for i in range(n_scenarios)]
ax.set_xticks(x_tick_positions)
ax.set_xticklabels(x_tick_labels, rotation=0)  
ax.set_xlabel('Scenario', fontsize=fontsize)


ax.yaxis.grid(True, linestyle='--', alpha=0.7, linewidth=1.5)


from matplotlib.patches import Patch

legend_elements = [
    Patch(facecolor=color_lns['box'], edgecolor=color_lns['edge'], label=r'mean $G^{\text{JFMS}}_{\text{LNS}}$'),
    Patch(facecolor=color_rule['box'], edgecolor=color_rule['edge'], label=r'$G^{\text{JFMS}}_{\text{RULE}}$')
]

ax.legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(1.0, 0.0), frameon=True)


plt.tight_layout(pad=0.1, h_pad=1, w_pad=0.1)  

plt.savefig('G-JFMS-LNS-RULE_Combined', dpi=300, bbox_inches='tight')
plt.show()

#%%
current_dir = os.getcwd()
folder_path=os.path.join(current_dir, "prop value")
df_out = pd.DataFrame(columns=['without q and b', 'with q', 'with b', 'n without q and b', 'n with q', 'n with b'])
# csv_files = glob(os.path.join(folder_path, '*.csv'))
scen=["P2M2V2 CT3.0", "P2M2V2 CT8.0", "P2M2V2 CT13.0"]
df={}
i=0
for scenario in scen:
    df[scenario]=pd.read_csv(os.path.join(folder_path, scenario+".csv"))
    a=(df[scenario]['time_CP'].mean())
    b=(df[scenario]['time_CP_q'].mean())
    c=(df[scenario]['time_CP_b'].mean())

    d=(df[scenario].query('time_CP <= 900').shape[0])
    e=(df[scenario].query('time_CP_q <= 900').shape[0])
    f=(df[scenario].query('time_CP_b <= 900').shape[0])
    df_out.loc[i]=[a,b,c,d,e,f]
    i=i+1

df_out.to_excel(
    'results for table 5.xlsx',
    sheet_name='Sheet1',
    index=False,
    engine='openpyxl'
)

#%%
x = [0,50,100,150,200]
y=[1656.1033198389416,1756.1033333158566,1867.4366531795608,2020.6866666181402,2227.4366531658707]

plt.figure(figsize=(10, 6))

plt.plot(x, y, color='blue', linestyle='-', linewidth=2, marker='o', markersize=5)

plt.xlabel('Distance between component loading point\n and charging station', fontsize=fontsize)
plt.ylabel('Objective value', fontsize=fontsize)

plt.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout(pad=0.1, h_pad=1, w_pad=0.1)
plt.savefig('lemma value',transparent=True)
plt.show()



