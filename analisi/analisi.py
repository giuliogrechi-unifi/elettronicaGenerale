#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedFormatter
import numpy as np
import scipy.fftpack
import scipy.signal
import tabulate
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#matplotlib.rc('font', **{'family' : "serif"})
#params= {'text.latex.preamble' : [r'\usepackage{libertine}']}
#plt.rcParams.update(params)
#matplotlib.rc('font',**{'family':'serif','serif':['Libertine']})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.style.use(['science'])
fig_size = 18/2.54, 7/2.54 #dnl inl alle varie frequenze
fig_size2 = 18/2.54, 14/2.54 #andamenti dnl inl
fig_size3 = 18/2.54, 7/2.54 #ricostruzione onde
directory = '.'
clock_period = (1/66.6)*1E-6 #in secondi
clock_period_err = (1E-4)*clock_period #100 ppm frequency stability
clock_period_arduino = (1/16.0)*1E-6 #in secondi
clock_period_arduino_err = (1E-4)*clock_period_arduino #100 ppm frequency stability
num_samples_arduino = 3500
bits_ADC = 10
V_range = 5000.0 #in mV
V_range = 5172.0 #in mV misurato con voltmetro, tensione usb (alimentazione arduino)
V_range_err = 28 #in mV
INL0_1 = (12.1/V_range)*1023 - 0.5#parametro da definire!
INL0_1_err = ((12.1*0.01 + 0.2)/V_range)*1023
INL_0_mV =INL0_1*V_range/1023
INL0 = 0

result_total = []
result_fit = []
result_bsl = []

#cerca nella directory path_to_dir tutti i file che hanno la giusta formattazione:
#il nome del file contiene già i parametri impostati durante la misura e la funzione
#ne estrai i valori numerici. Ritorna una lista [ [nomefile,freq,presc],[. . . ],...]
def find_mis(path_to_dir):
    fn  = os.listdir(path_to_dir)
    mis = []
    for file in fn:
        if 'dati_seriale' in file:
            param = file.strip('.txt').strip('dati_seriale-f').split('-p')
            m = [file,int(param[0]),int(param[1])]
            mis.append(m)
    mis = sorted(mis,key=getkey)
    #print mis
    #sys.exit()
    return mis

def getkey(item):
    return item[1]

#cerca nella lista wave_list il valore value e restituisce una lista 
#i cui elementi sono gli indici nei quali sono stati trovati
def search_code(wave_list,value):
    x_points = []
    for i, data in enumerate(wave_list):
        if data == value:
            x_points.append(i)
    return x_points

#legge il file e resituisce una lista di interi
def read_data_file(data_file):
    with open(data_file,'r') as f:
        data_list_int = [int(x) for x in f.readlines()]
    #with open(data_file,'r') as f:
    #    data_list_int_2 = [int(x) for x in f.readlines()]
    return data_list_int #+ data_list_int_2

#dati gli edges dell'istogramma, calcola il punto centrale per disegnare il 
#grafico
def list_central_points(bin_counts,bins_edges):
    bin_central = []
    for i,value in enumerate(bin_counts):
        central = (bins_edges[i] + bins_edges[i+1])/2.0
        bin_central.append(central)
    return bin_central

#dati una lista di interi che contiene l'onda e la frequenza di campionamento
#disegna nel subplot num_subplot la trasformata di fourier del segnale
#ritorna null
def plotFFT(data_list_int, samp_freq, num_subplot):
    plt.subplot(num_subplot)
    N = len(data_list_int)
    T = 1.0 / samp_freq
    x = range(len(data_list_int))
    y = data_list_int
    yf = scipy.fftpack.fft(y)
    xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
    plt.yscale('log')
    plt.grid(True)
    yfft = 2.0/N * np.abs(yf[:N//2])
    list_yfft = [x for x in yfft]
    fft_mx_index = list_yfft.index(max(yfft[1:]))
    plt.plot(xf, yfft, linewidth=0.5, markersize=2)
    plt.axvline(x=xf[fft_mx_index],color='r',alpha=0.50)
    print ("frequenza fondamentale FT: %s Hz"%(xf[fft_mx_index]))
    #print (yfft[10:])
    return


#data la lista di interi (che rappresenta il segnale) cerca tutti i codici compresi
#tra min e max, prendendo solo gli estremi trovati e calcolando la differenza,
#nel caso siano su un fronte dello stesso segno. Nel caso non trovi valori, restituisce
#comunque il primo dei dati che ha trovato: il risultato deve essere controllato 
#volta per volta dal plot
def codes_search_difference(data_list_int,min_srch,max_srch):
    diff = []
    codice = []
    maxmin_code = []
    for d in range(min_srch,max_srch):
        x_values_found = search_code(data_list_int,d)
        #print (x_values_found)
        if len(x_values_found) > 1:
            x1 = x_values_found[0]
            i = 1
            x2 = x_values_found[-i]
            sign1 = data_list_int[x1-5] - data_list_int[x1]
            sign2 = data_list_int[x2-5] - data_list_int[x2]
            if sign1 == 0 or sign2 == 0:
                pass
            else:
                sign1 = sign1/abs(sign1)
                sign2 = sign2/abs(sign2)
                while sign1 != sign2:
                    i = i + 1
                    x2 = x_values_found[-i]
                    sign1 = data_list_int[x1-5] - data_list_int[x1]
                    sign2 = data_list_int[x2-5] - data_list_int[x2]
                    if sign1 == 0 or sign2 == 0:
                        pass
                    else:
                        sign1 = sign1/abs(sign1)
                        sign2 = sign2/abs(sign2)
            codes = [x1,x2]
            difference = x2 - x1
            maxmin_code.append(codes)
            codice.append(d)
            diff.append(difference)
    return maxmin_code, codice, diff

def plot_wave2(data_list_int,max_list,codes,frq_set,presc_st):
    matplotlib.pyplot.close("all") #cambiare se necessario!!
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=fig_size3)
    t = range(len(data_list_int))
    axes.plot(t, data_list_int,
                    linewidth=0.5, markersize=1)
    if True: #len(max_list) > 1:
        for xvalue_zero in max_list:
            axes.axvline(x=xvalue_zero,color='g',linewidth=0.5,alpha=0.50)
        axes.axvline(x=codes[0],color='r',linewidth=0.5,alpha=0.50)
        axes.axvline(x=codes[1],color='r',linewidth=0.5,alpha=0.50)
    axes.set_xlabel(r'numero campione')
    axes.set_ylabel(r'codice ADC')
    fig.suptitle(r'Frequenza: $%.3f$ Hz - Prescaler: $%s$'%(frq_set,presc_st))
    if frq_set > 260.0 :
        axes.set_xlim([0,1000])
        if frq_set > 800.0 :
            axes.set_xlim([0,250])
            if frq_set > 2000.0 :
                axes.set_xlim([0,150])

    fig.savefig('img/_-_triang%.0f-p%s.pgf'%(frq_set,presc_st), dpi=300)
    return

def sampl_freq(data_file_to_read, num_samples_ar,count_period_clk_counts,frq_set,prescaler_set):
    period_wave = count_period_clk_counts * clock_period * 4095 * 2
    period_wave_err = count_period_clk_counts * clock_period_err * 4095 * 2
    data_list_int_all = read_data_file(data_file_to_read)
    data_list_int = data_list_int_all[:num_samples_ar]
    maxmin_code, codice, diff = codes_search_difference(data_list_int,10,1015)
    diff_max = max(diff)
    dif_mx_index = diff.index(diff_max)
    codes = maxmin_code[dif_mx_index]
    data_ls_restr = data_list_int[codes[0]:codes[1]]
    max_list, _ = scipy.signal.find_peaks(data_ls_restr,distance=3,height=950)
    max_list = [x+codes[0] for x in max_list]
    if (len(max_list)) == 0:
        print("errore nella misura della freq di campionamento")
        result_total.append(0)
    else:
        samp_freq = diff_max/(len(max_list)*period_wave)
        samp_freq_err = samp_freq*np.sqrt((2/diff_max)**2+(period_wave_err/period_wave)**2)
        print ("frequenza campionamento misurata: %.4f+-%.4fkHz"%(samp_freq/1000,samp_freq_err/1000))
        result_total.append(samp_freq/1000)
        result_total.append(samp_freq_err/1000)
        #plotFFT(data_list_int, samp_freq, 212)
    plot_wave2(data_list_int,max_list,maxmin_code[dif_mx_index],frq_set,prescaler_set)
    plot_graph2(data_list_int_all,frq_set,prescaler_set)

    #fig.savefig(file_to_save, dpi=fig.dpi)
    #plt.show()

def count_period_clk_counts(frequency):
    periodo_sec = 1.0/frequency
    periodo_clk_counts = periodo_sec/clock_period
    clock_per_count = int(periodo_clk_counts/(2*4095))
    count_clock = 1/(frequency * (4095*2*clock_period))
    periodo = clock_period*clock_per_count*1E3*4095.0*2.0 #in mssampl_freq
    periodo_err = clock_period_err*clock_per_count*1E3*4095.0*2.0 #in ms
    frequenza = 1000/periodo
    frequenza_err = frequenza*(periodo_err/periodo)
    print ("colpi clock ad ogni incremento: %s"%clock_per_count)
    print ("periodo: %.4f+-%.4f ms"%(periodo,periodo_err))
    print ("frequenza impostata: %.3f+-%.3f Hz"%(frequenza,frequenza_err))
    result_total.append(frequenza)
    result_total.append(frequenza_err)
    return clock_per_count,frequenza

def calc_DNL(list_bin_norm,n_tot):
    P_D = 1.0/len(list_bin_norm)
    DNL = []
    DNL_sig = []
    for Counts_D in list_bin_norm:
        DNL_D = (Counts_D - P_D)*V_range
        #DNL_sig_D = abs(DNL_D)/np.sqrt(Counts_D*n_tot) #d'alessandro
        #DNL_sig_D = V_range*np.sqrt(Counts_D/n_tot)   #nostro
        #DNL_sig_D = (V_range/1023)*(abs((DNL_D*1023/V_range) + 1)/np.sqrt(Counts_D*n_tot)) #ricciarini
        DNL_sig_D = (V_range)*np.sqrt(Counts_D/n_tot) #leandro
        #print (DNL_D)
        DNL.append(DNL_D) #dnl in millivolt
        DNL_sig.append(DNL_sig_D)
    return DNL,DNL_sig

def calc_INL(list_dnl,sig_dnl):
    DNL = list_dnl
    INL = []
    INL_sig = []
    for i_d,dnl in enumerate(DNL):
        INL_D = 0.0 #DNL[0]
        INL_var_D = 0.0
        for j in range(i_d + 1):
            INL_D = INL_D + DNL[j]
            INL_var_D = INL_var_D + sig_dnl[j]**2
        INL_D = INL_D - DNL[i_d]/2.0 - DNL[0]/2.0 + INL0
        #non faccio la somma in quadratura, perché sto togliendo dei pezzi a qualcosa
        #che esiste già.. altrimenti lo conto due volte (basta fare il conto)
        INL_var_D = INL_var_D - (sig_dnl[i_d]**2)/2 - (sig_dnl[0]**2)/2
        INL.append(INL_D)
        if INL_var_D < 0:
            print ('***************************************INL_var_D %s'%INL_var_D)
        INL_sig.append(np.sqrt(INL_var_D))
    return INL,INL_sig

def func(x,m,q):
    return m*x + q

def calc_AD(list_dnl,sig_dnl):
    DNL = list_dnl
    AD = []
    AD_sig = []
    for i_d,dnl in enumerate(DNL):
        AD_D = 0.0 #DNL[0]
        AD_var_D = 0.0
        for j in range(i_d + 1):
            AD_D = AD_D + DNL[j] + 1
            AD_var_D = AD_var_D + sig_dnl[j]**2
        AD_D = AD_D - DNL[i_d]/2.0 - DNL[0]/2.0 + INL0_1
        #non faccio la somma in quadratura, perché sto togliendo dei pezzi a qualcosa
        #che esiste già.. altrimenti lo conto due volte (basta fare il conto)
        AD_var_D = AD_var_D - (sig_dnl[i_d]**2)/2 - (sig_dnl[0]**2)/2 + INL0_1_err**2
        AD.append(AD_D)
        if AD_var_D < 0:
            print ('***************************************AD_var_D %s'%AD_var_D)
        AD_sig.append(np.sqrt(AD_var_D))
    return AD,AD_sig

def graphAD(ad, ad_sig,frq_set,presc_st):
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(18/2.54, 18/2.54))
    x = range(1,1023)
    y = [value*V_range/1023 for value in ad]
    y_err = [value*V_range/1023 for value in ad_sig]
    y_err[0] = min(y_err[1:])
    with open('data_AD/retta_cal_data_f%.0f_p%s.txt'%(frq_set,presc_st), 'w') as f:
        for i,data in enumerate(y):
            line_str = "%s %s %s\n"%(x[i],data,y_err[i])
            f.write(line_str)
    plt.errorbar(x, y, yerr=y_err, fmt=',',elinewidth=0.3,capsize=0)
    plt.suptitle('Retta di calibrazione - frequenza: %.2f - prescaler: %s'
                                                           %(frq_set,presc_st))
    plt.xlabel('codice')
    plt.ylabel('Tensione [mV]')
    #plt.show()
    x_a = np.array(x) 
    y_a = np.array(y)
    y_err_a = np.array(y_err)
    popt, pcov = curve_fit(func, x_a, y_a,sigma=y_err_a, p0=[5.2, INL_0_mV])
    perr = np.sqrt(np.diag(pcov))
    param_mis = [frq_set,presc_st]
    result_fit.append([param_mis,popt,perr])
    plt.plot(x_a,func(x_a,*popt),linewidth=0.3)
    residui = [abs(y[d] - (popt[0]*(d+1) + popt[1])) for d in range(0,1022)]
    residui_err= [np.sqrt(y_err[d]**2 + ((perr[0]**2)*(d+1) + perr[1]**2)) for d in range(0,1022)]
    INL_sl_volt = max(residui)
    INL_sl_volt_err = residui_err[residui.index(INL_sl_volt)]
    INL_sl_lsb = INL_sl_volt*1023/V_range
    INL_sl_lsb_err = INL_sl_volt_err*1023/V_range
    #print('INL con metodo best straight line: %.4f pm %.4f'%(INL_sl_lsb,INL_sl_lsb_err))
    result_bsl.append([frq_set,presc_st,INL_sl_lsb,INL_sl_lsb_err])
    axins = zoomed_inset_axes(axes, 10, loc=4) # zoom-factor: 2.5, location: upper-left
    axins.errorbar(x, y, yerr=y_err, fmt=',',elinewidth=0.7,capsize=0.5)
    axins.plot(x_a,func(x_a,*popt),linewidth=0.7, color='g')
    axins.plot(x_a,[d/0.197796 for d in x_a],linewidth=0.7,color='r')
    mark_inset(axes, axins, loc1=1, loc2=2, fc="none", ec="0.5")
    x1, x2, y1, y2 = 980, 1035, 4980, 5220 # specify the limits
    axins.set_xlim(x1, x2) # apply the x-limits
    axins.set_ylim(y1, y2) # apply the y-limits
    plt.yticks(visible=False)
    plt.xticks(visible=False)

    axins_2 = zoomed_inset_axes(axes, 10, loc=2) # zoom-factor: 2.5, location: upper-left
    axins_2.errorbar(x, y, yerr=y_err, fmt=',',elinewidth=0.7,capsize=0.5)
    axins_2.plot(x_a,func(x_a,*popt),linewidth=0.7, color='g')
    axins_2.plot(x_a,[d/0.197796 for d in x_a],linewidth=0.7,color='r')
    mark_inset(axes, axins_2, loc1=3, loc2=4, fc="none", ec="0.5")
    x1, x2, y1, y2 = -10, 30, -20, 150 # specify the limits
    axins_2.set_xlim(x1, x2) # apply the x-limits
    axins_2.set_ylim(y1, y2) # apply the y-limits
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    fig.savefig('img_AD/_-_retta_cal_%.0f-p%s.pgf'%(frq_set,presc_st), dpi=300)
    fig.savefig('img_AD/_-_retta_cal_%.0f-p%s.png'%(frq_set,presc_st), dpi=300)
    #plt.show()

def plot_graph2(int_list,frq_set,presc_st):
    matplotlib.pyplot.close("all") #cambiare se necessario!!
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=fig_size)
    plt.subplots_adjust(wspace=0.02,bottom=0.16,top=0.89,right=0.88,left=0.10)
    bins_list = list(range(1,1024)) #da 1 a 1023 compreso
    for i,x in enumerate(bins_list): #da 0.5 a 1022.5 compresi
        bins_list[i] = x - 0.5
    bin_values, bins_edges, null= axes[0].hist(int_list,
                                               bins=bins_list, 
                                               density=False,
                                               alpha=0.80)
    #num_bins = len(bin_values)
    N_tot =  np.sum(bin_values) #numero totale di eventi
    bin_values_orig = bin_values
    bin_values = [x/N_tot for x in bin_values]
    #print ('prima: %s'%np.sum(bin_values_orig))
    #print ('dopo: %s'%np.sum(bin_values))
    #N_tot =  np.sum(bin_values) #numero totale di eventi
    DNL,DNL_sig = calc_DNL(bin_values,N_tot) #in mV
    DNL_LSB = [x/(V_range/1023) for x in DNL] #trasformazione in LSBsampl_freq
    DNL_sig_LSB = [x/(V_range/1023) for x in DNL_sig] #trasformazione in LSB
    INL,INL_sig = calc_INL(DNL_LSB,DNL_sig_LSB)            #gia' in LSB (credo... ripensarci)
    AD,AD_sig = calc_AD(DNL_LSB,DNL_sig_LSB)
    DNL_max = max([abs(x) for x in DNL_LSB])
    try:
        DNL_max_index = DNL_LSB.index(DNL_max)
    except ValueError:
        DNL_max_index = DNL_LSB.index(-DNL_max)
    DNL_max_sig = DNL_sig_LSB[DNL_max_index]
    #print ('errore DNL: %s percento'%((100*DNL_max_sig)/DNL_max))

    INL_max = max([abs(x) for x in INL])
    try:
        INL_max_index = INL.index(INL_max)
    except ValueError:
        INL_max_index = INL.index(-INL_max)
    INL_max_sig = INL_sig[INL_max_index]
    print ("DNL: %.2f+-%.2f lsb"%(DNL_max, DNL_max_sig))
    print ("INL: %.2f+-%.2f lsb"%(INL_max, INL_max_sig))
    result_total.append(DNL_max)
    result_total.append(DNL_max_sig)
    result_total.append(INL_max)
    result_total.append(INL_max_sig)
    
    bin_central = list_central_points(bin_values,bins_edges)
    axes[0].clear()
    axes[0].plot(bin_central,DNL_LSB,
                    linewidth=0.2, markersize=0,
                    alpha=0.80)
    axes[1].plot(bin_central,INL,
                    linewidth=0.2, markersize=0,
                    alpha=0.8)
    if frq_set > 11:
        axes[0].set_ylim([-0.75,0.75])
        axes[1].set_ylim([-2.0,2.0])
        if frq_set < 6:
            axes[1].set_ylim([-9,5])
        elif frq_set < 11:
            axes[1].set_ylim([-1,6])
        elif frq_set > 800:
            axes[0].set_ylim([-2,2])
            axes[1].set_ylim([-2,6])
            if frq_set > 2000:
                axes[1].set_ylim([-2,50])
    else:
        axes[0].set_ylim([-1.7,1.7])
        axes[1].set_ylim([-15,15])
    axes[1].set_xlabel(r'Codice\ ADC')
    axes[0].set_xlabel(r'Codice\ ADC')
    axes[1].set_ylabel(r'INL\ [$lsb$]')
    axes[0].set_ylabel(r'DNL\ [$lsb$]')
    ticks0 = axes[0].get_yticks()
    ticks1 = axes[1].get_yticks()
    ll0 = ['%.2f' % a for a in ticks0]
    ll1 = ['%.2f' % a for a in ticks1]
    axes[0].yaxis.set_major_formatter(FixedFormatter(ll0))
    axes[1].yaxis.set_major_formatter(FixedFormatter(ll1))
    #axes[0].set_yticklabels(ll0)
    #axes[1].set_yticklabels(ll1)
    fig.suptitle(r'Frequenza: $%.3f$ Hz - Prescaler: $%s$'%(frq_set,presc_st))
    #axes[1][0].set_title('Probabilita')
    #axes[0].set_title('DNL')
    #axes[1].set_title('INL')
    #axes[1].yaxis.tick_right()
    axes[1].tick_params(axis='y',labelleft=False, labelright=True)
    axes[1].yaxis.set_label_position("right")
    #matplotlib.rcParams['ytick.left'] = True
    #matplotlib.rcParams['ytick.right'] = True
    fig.savefig('img/_-_dnlINL_f%.0f-p%s.pgf'%(frq_set,presc_st), dpi=300)
    #plt.show()
    matplotlib.pyplot.close("all") #cambiare se necessario!!
    graphAD(AD,AD_sig,frq_set,presc_st)

def plot_graph(int_list,frq_set):
    fig, axes = plt.subplots(2,2)
    bins_list = list(range(1,1024)) #da 1 a 1023 compreso
    for i,x in enumerate(bins_list): #da 0.5 a 1022.5 compresi
        bins_list[i] = x - 0.5
    bin_values, bins_edges, null= axes[1][0].hist(int_list,
                                               bins=bins_list, 
                                               density=False,
                                               alpha=0.80)
    #num_bins = len(bin_values)
    N_tot =  np.sum(bin_values) #numero totale di eventi
    bin_values_orig = bin_values
    bin_values = [x/N_tot for x in bin_values]
    #print ('prima: %s'%np.sum(bin_values_orig))
    #print ('dopo: %s'%np.sum(bin_values))
    #N_tot =  np.sum(bin_values) #numero totale di eventi
    DNL,DNL_sig = calc_DNL(bin_values,N_tot) #in mV
    DNL_LSB = [x/(V_range/1023) for x in DNL] #trasformazione in LSB
    DNL_sig_LSB = [x/(V_range/1023) for x in DNL_sig] #trasformazione in LSB
    INL,INL_sig = calc_INL(DNL_LSB,DNL_sig_LSB)            #gia' in LSB (credo... ripensarci)
    DNL_max = max([abs(x) for x in DNL_LSB])
    try:
        DNL_max_index = DNL_LSB.index(DNL_max)
    except ValueError:
        DNL_max_index = DNL_LSB.index(-DNL_max)
    DNL_max_sig = DNL_sig_LSB[DNL_max_index]
    #print ('errore DNL: %s percento'%((100*DNL_max_sig)/DNL_max))

    INL_max = max([abs(x) for x in INL])
    try:
        INL_max_index = INL.index(INL_max)
    except ValueError:
        INL_max_index = INL.index(-INL_max)
    INL_max_sig = INL_sig[INL_max_index]
    print ("DNL: %.2f+-%.2f lsb"%(DNL_max, DNL_max_sig))
    print ("INL: %.2f+-%.2f lsb"%(INL_max, INL_max_sig))
    result_total.append(DNL_max)
    result_total.append(DNL_max_sig)
    result_total.append(INL_max)
    result_total.append(INL_max_sig)
    matplotlib.pyplot.close("all") #cambiare se necessario!!
    return
    bin_central = list_central_points(bin_values,bins_edges)

    axes[0][0].plot(list(range(len(int_list))), int_list,
                    linewidth=1, markersize=2,
                    alpha=0.80)
    axes[0][1].plot(bin_central,DNL_LSB,
                    linewidth=1, markersize=2,
                    alpha=0.80)
    axes[1][1].plot(bin_central,INL,
                    linewidth=1, markersize=2,
                    alpha=0.80)
    axes[1][0].set_ylabel('Probabilita')
    axes[1][0].set_title('istogramma normalizzato')
    axes[1][0].set_xlabel('Codice ADC')
    axes[1][1].set_xlabel('Codice ADC')
    axes[1][1].set_ylabel('LSB')
    axes[0][1].set_ylabel('LSB')
    axes[0][0].set_ylabel('codice adc')
    #axes[0][1].set_ylabel('DNL')
    #axes[1][0].set_title('Probabilita')
    axes[0][1].set_title('DNL')
    axes[1][1].set_title('INL')
    #fig.savefig(file_to_save, dpi=fig.dpi)

    #plt.show()

def format_result(list_res):
    list = []
    #print('numero dati nella lista %s'%len(list_res))
    for m in range(0,len(list_res),12):
        #print m
        data = []
        for i in range(12):
            #print m+i
            data.append(list_res[m+i])
        list.append(data)
    return list

if __name__ == "__main__":
    try:
        misure = find_mis(directory)
        numero_misure = len(misure)
        #print ('numero misure: %s'%len(misure))
        for i,m in enumerate(misure):
            print('')
            print (m[0])
            print('freq: %s - prescaler: %s'%(m[1],m[2]))
            data_fil = m[0]
            prescaler = m[2]
            frq = m[1]
            result_total.append(m[0])
            result_total.append(m[1])
            result_total.append(m[2])
            clk_cnt,freq_impost = count_period_clk_counts(int(frq))
            result_total.append(clk_cnt)
            sampl_freq(data_fil,num_samples_arduino,clk_cnt,freq_impost,prescaler)
            #print result_total
            #choice = raw_input('prosegui (y/n)')
            #if choice != 'n':
            #    pass
            #else:
            #    sys.exit(0)
        #print result_total
        result = format_result(result_total)
        prescaler_list = [128,64,32]
        p128 = []
        p64 = []
        p32 = []
        freq_128 = []
        freq_64  = []
        freq_32 = []
        freq_128_err = []
        freq_64_err  = []
        freq_32_err = []
        val_128 = []
        val_64  = []
        val_32 = []
        val_128_err = []
        val_64_err  = []
        val_32_err = []

        for r in result:
            #print r
            if r[2] == prescaler_list[0]:
                freq_128.append(r[3])
                freq_128_err.append(r[4])
                val_128.append([r[8],r[10],r[6]])
                val_128_err.append([r[9],r[11],r[7]])
                p128.append(r)
            elif r[2] == prescaler_list[1]:
                p64.append(r)
                freq_64.append(r[3])
                freq_64_err.append(r[4])
                val_64.append([r[8],r[10],r[6]])
                val_64_err.append([r[9],r[11],r[7]])
            elif r[2] == prescaler_list[2]:
                p32.append(r)
                freq_32.append(r[3])
                freq_32_err.append(r[4])
                val_32.append([r[8],r[10],r[6]])
                val_32_err.append([r[9],r[11],r[7]])

        #print(freq_128[3:9])
        dnl_average_128 = np.average([x[0] for x in val_128[3:9]], weights=[1/x[0]**2 for x in val_128_err[3:9]]) #media pesata
        dnl_average_128_err = np.sqrt(1/sum([1/x[0]**2 for x in val_128_err[3:9]]))
        print('dnl medio 128 %.3f+-%.4f'%(dnl_average_128,dnl_average_128_err))
        inl_average_128 = np.average([x[1] for x in val_128[3:9]], weights=[1/x[1]**2 for x in val_128_err[3:9]]) #media pesata
        inl_average_128_err = np.sqrt(1/sum([1/x[1]**2 for x in val_128_err[3:9]]))
        print('inl medio 128 %.3f+-%.4f'%(inl_average_128,inl_average_128_err))

        dnl_average_64 = np.average([x[0] for x in val_64], weights=[1/x[0]**2 for x in val_64_err]) #media pesata
        dnl_average_64_err = np.sqrt(1/sum([1/x[0]**2 for x in val_64_err]))
        print('dnl medio 64 %.3f+-%.4f'%(dnl_average_64,dnl_average_64_err))
        inl_average_64 = np.average([x[1] for x in val_64], weights=[1/x[1]**2 for x in val_64_err]) #media pesata
        inl_average_64_err = np.sqrt(1/sum([1/x[1]**2 for x in val_64_err]))
        print('inl medio 64 %.3f+-%.4f'%(inl_average_64,inl_average_64_err))

        dnl_average_32 = np.average([x[0] for x in val_32], weights=[1/x[0]**2 for x in val_32_err]) #media pesata
        dnl_average_32_err = np.sqrt(1/sum([1/x[0]**2 for x in val_32_err]))
        print('dnl medio 32 %.3f+-%.4f'%(dnl_average_32,dnl_average_32_err))
        inl_average_32 = np.average([x[1] for x in val_32], weights=[1/x[1]**2 for x in val_32_err]) #media pesata
        inl_average_32_err = np.sqrt(1/sum([1/x[1]**2 for x in val_32_err]))
        print('inl medio 32 %.3f+-%.4f'%(inl_average_32,inl_average_32_err))

        #sys.exit()

        table = []
        print ('')
        print('prescaler 128')
        smp_frq = []
        smp_frq_err = []
        table.append(['Freq. onda [Hz]','Freq. campionamento [kHz]', 'DNL [lsb]','INL [lsb]'])
        table.append(['prescaler 128'])
        table.append(['\hline'])
        for i,v in enumerate(val_128):
            print('f: %.4f+-%.4f   \t sr: %.5f+-%.5fHz \t DNL: %.3f+-%.3f \t INL: %.3f+-%.3f'%(
                       freq_128[i], freq_128_err[i],
                       v[2], val_128_err[i][2],
                       v[0], val_128_err[i][0],
                       v[1], val_128_err[i][1]))
            data = [
                    r'$%.3f\pm%.3f$'%(freq_128[i], freq_128_err[i]),
                    r'$%.3f\pm%.3f$'%(v[2], val_128_err[i][2]),
                    r'$%.2f\pm%.2f$'%(v[0], val_128_err[i][0]),
                    r'$%.1f\pm%.1f$'%(v[1], val_128_err[i][1])
                    ]
            table.append(data)
            if v[2] > 7 and v[2] < 12:
                smp_frq.append(v[2])
                smp_frq_err.append(val_128_err[i][2])
        smp_average_128 = np.average(smp_frq, weights=[1/x**2 for x in smp_frq_err]) #media pesata
        smp_average_128_err = np.sqrt(1/sum([1/x**2 for x in smp_frq_err]))
        print('frequenza di campionamento media:   %.4f +- %.4f'%(smp_average_128,smp_average_128_err))
        smp_teo_128 = 1/(1000*clock_period_arduino*128*13)
        smp_teo_128_err = smp_teo_128*(clock_period_arduino_err/clock_period_arduino)
        print('frequenza di campionamento teorica: %.4f +- %.4f'%(smp_teo_128,smp_teo_128_err))

        print('prescaler 64')
        smp_frq = []
        smp_frq_err = []
        table.append([''])
        table.append(['\hline'])
        table.append(['prescaler 64'])
        table.append(['\hline'])
        for i,v in enumerate(val_64):
            print('f: %.4f+-%.4f   \t sr: %.5f+-%.5fHz \t DNL: %.3f+-%.3f \t INL: %.3f+-%.3f'%(
                        freq_64[i], freq_64_err[i], 
                        v[2], val_64_err[i][2],
                        v[0], val_64_err[i][0],
                        v[1], val_64_err[i][1]))
            data = [
                    r'$%.3f\pm%.3f$'%(freq_64[i], freq_64_err[i]),
                    r'$%.3f\pm%.3f$'%(v[2], val_64_err[i][2]),
                    r'$%.2f\pm%.2f$'%(v[0], val_64_err[i][0]),
                    r'$%.1f\pm%.1f$'%(v[1], val_64_err[i][1])
                    ]
            table.append(data)
            if v[2] > 17 and v[2] < 21:
                smp_frq.append(v[2])
                smp_frq_err.append(val_64_err[i][2])
        smp_average_64 = np.average(smp_frq, weights=[1/x**2 for x in smp_frq_err]) #media pesata
        smp_average_64_err = np.sqrt(1/sum([1/x**2 for x in smp_frq_err]))
        print('frequenza di campionamento media: %.4f +- %.4f'%(smp_average_64,smp_average_64_err))
        smp_teo_64 = 1/(1000*clock_period_arduino*64*13)
        smp_teo_64_err = smp_teo_64*(clock_period_arduino_err/clock_period_arduino)
        print('frequenza di campionamento teorica: %.4f +- %.4f'%(smp_teo_64,smp_teo_64_err))
        
        print('prescaler 32')
        smp_frq = []
        smp_frq_err = []
        table.append([''])
        table.append(['\hline'])
        table.append(['prescaler 32'])
        table.append(['\hline'])
        for i,v in enumerate(val_32):
            print('f: %.4f+-%.4f   \t sr: %.5f+-%.5fHz \t DNL: %.3f+-%.3f \t INL: %.3f+-%.3f'%(
                       freq_32[i], freq_32_err[i],
                        v[2], val_32_err[i][2],
                        v[0], val_32_err[i][0],
                        v[1], val_32_err[i][1]))
            data = [
                    r'$%.3f\pm%.3f$'%(freq_32[i], freq_32_err[i]),
                    r'$%.3f\pm%.3f$'%(v[2], val_32_err[i][2]),
                    r'$%.2f\pm%.2f$'%(v[0], val_32_err[i][0]),
                    r'$%.1f\pm%.1f$'%(v[1], val_32_err[i][1])
                    ]
            table.append(data)
            if v[2] > 35 and v[2] < 41:
                smp_frq.append(v[2])
                smp_frq_err.append(val_32_err[i][2])
        smp_average_32 = np.average(smp_frq, weights=[1/x**2 for x in smp_frq_err]) #media pesata
        smp_average_32_err = np.sqrt(1/sum([1/x**2 for x in smp_frq_err]))
        print('frequenza di campionamento media: %.4f +- %.4f'%(smp_average_32,smp_average_32_err))
        smp_teo_32 = 1/(1000*clock_period_arduino*32*13)
        smp_teo_32_err = smp_teo_32*(clock_period_arduino_err/clock_period_arduino)
        print('frequenza di campionamento teorica: %.3f +- %.3f'%(smp_teo_32,smp_teo_32_err))

        #print( tabulate.tabulate(table, tablefmt="latex_raw"))

        matplotlib.pyplot.close("all")
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=fig_size2,sharex=True)
        plt.subplots_adjust(hspace=0.07,bottom=0.11,top=0.94,right=0.94,left=0.10)
        #fig = plt.figure(figsize=fig_size2)
        #plt.subplot(211)
        axes[0].errorbar(freq_128, [dnl[0] for dnl in val_128], yerr=[err[0] for err in val_128_err],
                     fmt='_',elinewidth=0.5,capsize=1.5,label=r'$%.2f$ kHz'%smp_average_128)
        axes[0].errorbar(freq_64, [dnl[0] for dnl in val_64],  yerr=[err[0] for err in val_64_err],
                     fmt='_',elinewidth=0.5,capsize=1.5,label=r'$%.2f$ kHz'%smp_average_64)
        axes[0].errorbar(freq_32, [dnl[0] for dnl in val_32],  yerr=[err[0] for err in val_32_err],
                     fmt='_',elinewidth=0.5,capsize=1.5,label=r'$%.2f$ kHz'%smp_average_32)
        #axes = plt.gca()
        axes[0].axhline(y=0.5,color='r',linewidth=0.5)
        #axes.set_title('DNL')
        axes[0].set_ylabel(r'DNL [$lsb$]')
        #axes[0].set_xlabel(r'$frequenza onda [Hz]$')
        yticks = [0.0,0.5,1.0,1.5,2.0,2.5,3.0]
        axes[0].set_yticks(yticks)
        x_ticks = [10,50,100,150,254,406,625,813]
        #axes[0].set_xticks([0.0 + i for i in range(0,1000,100)]+[50])
        #axes[0].set_xticks(x_ticks)
        axes[0].set_ylim([0,2.5])
        axes[0].set_xlim([-40,950])
        axes[0].legend(loc=4,framealpha=1,fancybox=False,frameon=True,
                       markerscale=0.7,title=r'freq. campionamento')
        #plt.subplot(212)
        """
        axes[1].errorbar(freq_128, [inl[1] for inl in val_128], yerr=[err[1] for err in val_128_err],
                     fmt='x',elinewidth=0.5,capsize=1,label=r'$%.2f$ kHz'%smp_average_128)
        axes[1].errorbar(freq_64, [inl[1] for inl in val_64],  yerr=[err[1] for err in val_64_err],
                     fmt='x',elinewidth=0.5,capsize=1,label=r'$%.2f$ kHz'%smp_average_64)
        axes[1].errorbar(freq_32, [inl[1] for inl in val_32],  yerr=[err[1] for err in val_32_err],
                     fmt='x',elinewidth=0.5,capsize=1,label=r'$%.2f$ kHz'%smp_average_32)
        """
        axes[1].errorbar(freq_128, [inl[1] for inl in val_128],
                     fmt='x',label=r'$%.2f$ kHz'%smp_average_128)
        axes[1].errorbar(freq_64, [inl[1] for inl in val_64],
                     fmt='4',label=r'$%.2f$ kHz'%smp_average_64)
        axes[1].errorbar(freq_32, [inl[1] for inl in val_32],
                     fmt='+',label=r'$%.2f$ kHz'%smp_average_32)
        #axes = plt.gca()
        axes[1].axhline(y=1.25,color='r',linewidth=0.5)
        #axes.set_title('INL')
        axes[1].set_ylabel(r'INL [$lsb$]')
        axes[1].set_xlabel(r'frequenza onda [$Hz$]')
        yticks = [0.0,1.0,2.0,3.0,4.0]
        axes[1].set_yticks(yticks)
        ticks0 = axes[0].get_yticks()
        ticks1 = axes[1].get_yticks()
        ll0 = ['%2.2f' % a for a in ticks0]
        ll1 = ['%2.2f' % a for a in ticks1]
        axes[0].yaxis.set_major_formatter(FixedFormatter(ll0))
        axes[1].yaxis.set_major_formatter(FixedFormatter(ll1))
        #axes[1].set_xticks([0.0 + i for i in range(0,1000,100)]+[50])
        #axes[1].set_xticks(x_ticks)
        axes[1].set_xlim([-40,950])
        axes[1].set_ylim([0,3.5])
        axes[1].legend(loc=4,framealpha=1,markerscale=0.7,fancybox=False,frameon=True,
                       title=r'freq. campionamento')
        #plt.subplot(313)
        #plt.errorbar(freq_128, [frq[2] for frq in val_128], yerr=[err[2] for err in val_128_err],
        #             fmt='v',elinewidth=1,capsize=0.5,label='samp_freq p128')
        #plt.errorbar(freq_64, [frq[2] for frq in val_64],   yerr=[err[2] for err in val_64_err],
        #             fmt='v',elinewidth=1,capsize=0.5,label='samp_freq p64')
        #plt.errorbar(freq_32, [frq[2] for frq in val_32],   yerr=[err[2] for err in val_32_err],
        #             fmt='v',elinewidth=1,capsize=0.5,label='samp_freq p32')
        #axes = plt.gca()
        #axes.set_xlim([0,750])
        #plt.legend()
        fig.savefig('img/_-_dnl_inl-0-850.pgf', dpi=300)
        #plt.show()
        m128 = []
        m128_err = []
        q128 = []
        q128_err = []
        m64 = []
        m64_err = []
        q64 = []
        q64_err = []
        m32 = []
        m32_err = []
        q32 = []
        q32_err = []
        for data in result_fit:
            if data[0][1] == 128:
                        m128.append(data[1][0])
                        m128_err.append(data[2][0])
                        q128.append(data[1][1])
                        q128_err.append(data[2][1])
            if data[0][1] == 64:
                        m64.append(data[1][0])
                        m64_err.append(data[2][0])
                        q64.append(data[1][1])
                        q64_err.append(data[2][1])
            if data[0][1] == 32:
                        m32.append(data[1][0])
                        m32_err.append(data[2][0])
                        q32.append(data[1][1])
                        q32_err.append(data[2][1])
            #print ('frequenza: %.2f - prescaler: %s '%(data[0][0],data[0][1]))
            #print ('m: %s \t pm \t %s'%(data[1][0],data[2][0]))
            #print ('q: %s \t pm \t %s'%(data[1][1],data[2][1]))

        m_fit_average128 = np.average(m128[3:9], weights=[1/x**2 for x in m128_err[3:9]]) #media pesata
        m_fit_average128_err = np.sqrt(1/sum([1/x**2 for x in m128_err]))
        q_fit_average128 = np.average(q128[3:9], weights=[1/x**2 for x in q128_err[3:9]]) #media pesata
        q_fit_average128_err = np.sqrt(1/sum([1/x**2 for x in q128_err[3:9]]))
        m_fit_average64 = np.average(m64, weights=[1/x**2 for x in m64_err]) #media pesata
        m_fit_average64_err = np.sqrt(1/sum([1/x**2 for x in m64_err]))
        q_fit_average64 = np.average(q64, weights=[1/x**2 for x in q64_err]) #media pesata
        q_fit_average64_err = np.sqrt(1/sum([1/x**2 for x in q64_err]))
        m_fit_average32 = np.average(m32, weights=[1/x**2 for x in m32_err]) #media pesata
        m_fit_average32_err = np.sqrt(1/sum([1/x**2 for x in m32_err]))
        q_fit_average32 = np.average(q32, weights=[1/x**2 for x in q32_err]) #media pesata
        q_fit_average32_err = np.sqrt(1/sum([1/x**2 for x in q32_err]))

        print ('128')
        print ('m: %.5f \t pm \t %.5f'%(m_fit_average128,m_fit_average128_err))
        print ('q: %.3f \t pm \t %.3f'%(q_fit_average128,q_fit_average128_err))
        V_1023_128 = 1023*m_fit_average128 + q_fit_average128
        V_1023_128_err = np.sqrt((1023*m_fit_average128_err)**2 + q_fit_average128_err**2)
        gain_errore_128 = (-V_range + V_1023_128)*1023/V_range
        gain_errore_128_err = V_1023_128_err*1023/V_range
        print ('V(1023) = %.2f pm %.2f'%(V_1023_128,V_1023_128_err))
        print ('gain error: %.3f pm %.3f'%(gain_errore_128,gain_errore_128_err))

        print ('64')
        print ('m: %.5f \t pm \t %.5f'%(m_fit_average64,m_fit_average64_err))
        print ('q: %.3f \t pm \t %.3f'%(q_fit_average64,q_fit_average64_err))
        V_1023_64 = 1023*m_fit_average64 + q_fit_average64
        V_1023_64_err = np.sqrt((1023*m_fit_average64_err)**2 + q_fit_average64_err**2)
        gain_errore_64 = (-V_range + V_1023_64)*1023/V_range
        gain_errore_64_err = V_1023_64_err*1023/V_range
        print ('V(1023) = %.2f pm %.2f'%(V_1023_64,V_1023_64_err))
        print ('gain error: %.3f pm %.3f'%(gain_errore_64,gain_errore_64_err))

        print ('32')
        print ('m: %.5f \t pm \t %.5f'%(m_fit_average32,m_fit_average32_err))
        print ('q: %.3f \t pm \t %.3f'%(q_fit_average32,q_fit_average32_err))
        V_1023_32 = 1023*m_fit_average32 + q_fit_average32
        V_1023_32_err = np.sqrt((1023*m_fit_average32_err)**2 + q_fit_average32_err**2)
        gain_errore_32 = (-V_range + V_1023_32)*1023/V_range
        gain_errore_32_err = V_1023_32_err*1023/V_range
        print ('V(1023) = %.2f pm %.2f'%(V_1023_32,V_1023_32_err))
        print ('gain error: %.3f pm %.3f'%(gain_errore_32,gain_errore_32_err))

        inl_bl_128 = []
        inl_bl_128_err = []
        print ('prescaler 128')
        for i,data in enumerate(result_bsl):
            if data[1] == 128:
                print ('frequenza %.2f'%(data[0]))
                print('INL con metodo best straight line: %.1f pm %.1f'%(data[2],data[3]))
                inl_bl_128.append(data[2])
                inl_bl_128_err.append(data[3])
        inl_bl_128_average = np.average(inl_bl_128[3:9], weights=[1/x**2 for x in inl_bl_128_err[3:9]]) #media pesata
        inl_bl_128_average_err = np.sqrt(1/sum([1/x**2 for x in inl_bl_128_err]))
        print('media pesata tra i valori %.2f pm %.2f'%(inl_bl_128_average,inl_bl_128_average_err))        

        inl_bl_64 = []
        inl_bl_64_err = []       
        print ('prescaler 64')
        for i,data in enumerate(result_bsl):
            if data[1] == 64:
                print ('frequenza %.2f'%(data[0]))
                print('INL con metodo best straight line: %.1f pm %.1f'%(data[2],data[3]))
                inl_bl_64.append(data[2])
                inl_bl_64_err.append(data[3])
        inl_bl_64_average = np.average(inl_bl_64[3:9], weights=[1/x**2 for x in inl_bl_64_err[3:9]]) #media pesata
        inl_bl_64_average_err = np.sqrt(1/sum([1/x**2 for x in inl_bl_64_err]))
        print('media pesata tra i valori %.2f pm %.2f'%(inl_bl_64_average,inl_bl_64_average_err))

        inl_bl_32 = []
        inl_bl_32_err = []   
        print ('prescaler 32')
        for i,data in enumerate(result_bsl):
            if data[1] == 32:
                print ('frequenza %.2f'%(data[0]))
                print('INL con metodo best straight line: %.1f pm %.1f '%(data[2],data[3]))
                inl_bl_32.append(data[2])
                inl_bl_32_err.append(data[3])
        inl_bl_32_average = np.average(inl_bl_32[3:9], weights=[1/x**2 for x in inl_bl_32_err[3:9]]) #media pesata
        inl_bl_32_average_err = np.sqrt(1/sum([1/x**2 for x in inl_bl_32_err]))                
        print('media pesata tra i valori %.2f pm %.2f'%(inl_bl_32_average,inl_bl_32_average_err))        



    except KeyboardInterrupt:
        print ('')
        sys.exit(0)
    sys.exit(0)

