import serial
import sys
import time
import glob
import serial.tools.list_ports
import numpy as np
import matplotlib.pyplot as plt
import datetime

"""
inizializzazione delle variabili globali
"""
clock_period = (1/66.6)*1E-6 #periodo del master clock - in secondi
shap = '11'                  #forma del segnale impostato '11' -> triangolare
ampl_scelta = 7              #ampiezza massima del segnale tra 0 7
BAUDRATE = 57600
num_samples_arduino = 3500
bits_ADC = 10
V_range = 5000.0 #in mV


"""
valori di frequenza sui cui iterare - in kHz
numero di campionamenti da acquisire (diviso 100000)
prescaler su cui iterare
"""
freq_list =            [50,100,150,250,400,600]
num_points_acq_list =  [2 ,2  ,1  ,1  ,1  ,1  ]
prescaler_list = [128,64,32]

"""
creazione delle liste nella forma adatta alla successiva analisi
"""
freq_list = freq_list[::-1] #inversione dell'ordine della lista frequenze
num_points_acq_list = num_points_acq_list[::-1]
num_points_acq_list =  [x*100000 for x in num_points_acq_list]

"""
generazione dei valori da inviare via seriale.
Richiesti in ingresso i parametri dell'onda, cioè forma ampiezza e frequenza
restituisce in uscita una lista con i 3 byte formattati come richiesto 
dal software implementato nella CPLD
"""
def generate_bytes(clock_period,shape,ampl,frequency):
    bin_ampl = format(ampl, '03b')
    byte1 = '0b00%s0%s'%(shape,bin_ampl)
    byte1_int = int(byte1,2)
    #frequency = float(input("frequenza onda (2 - 8000 Hz): "))
    periodo_sec = 1.0/frequency
    periodo_clk_counts = periodo_sec/clock_period
    clock_per_count = int(periodo_clk_counts/(2*4095))
    count_clock = 1/(frequency * (4095*2*clock_period))
    periodo_scelta = clock_per_count - 1 
    byte3_int = periodo_scelta >> 6
    byte2_int = periodo_scelta - (byte3_int << 6)
    #print (frequency)
    #print ("byte1: %s"%byte1_int)
    #print ("byte2: %s"%byte2_int)
    #print ("byte3: %s"%byte3_int)
    bytes = [chr(byte1_int),chr(byte2_int),chr(byte3_int)]
    periodo = clock_period*clock_per_count*1E3*4095.0*2.0 #in ms
    print ("colpi clock ad ogni incremento: %s"%clock_per_count)
    print ("periodo: %f ms"%periodo)
    print ("frequenza impostata: %f kHz"%(1/periodo))
    return bytes,clock_per_count

"""
trova la porta seriale al quale è connesso il dispositivo
richiede in ingresso la stringa che identifica il produttore
maker_str = Arduino oppure maker_str = FTDI
resistuisce una stringa contenente il dispositivo seriale usato
"""
def find_port(maker_str):
    ports = list(serial.tools.list_ports.comports())
    port = 0
    for p in ports:
        if maker_str in p.manufacturer:
            port = p.device
    if port != 0:
        return port
    else:
        print ("errore apertura seriale %s"%maker_str)
        sys.exit(1)

"""
apre il file con le misure, effettua un casting a interi 
e restituisce una lista (di interi)
"""
def read_data_file(data_file):
    with open(data_file,'r') as f:
        data_list_int = [int(x) for x in f.readlines()]
    return data_list_int

"""
calcola il DNL a partire da una lista di interi, che rappresenta
il numero di conteggi, normalizzati, di un istogramma.
non calcola l'errore associato
"""
def calc_DNL(list_bin_norm):
    P_D = 1.0/len(list_bin_norm)
    DNL = []
    for Counts_D in list_bin_norm:
        DNL_D = (Counts_D - P_D)*V_range
        #print (DNL_D)
        DNL.append(DNL_D) #dnl in millivolt
    return DNL

"""
calcola l'INL a partire da una lista di interi, che rappresenta
il numero di conteggi, normalizzati, di un istogramma.
non calcola l'errore associato
"""
INL0 = 0 #ponendo INL0 si utilizza il metodo fixed end point per il calcolo di INL
def calc_INL(list_dnl):
    DNL = list_dnl
    INL = []
    for i_d,dnl in enumerate(DNL):
        INL_D = 0.0
        for j in range(i_d):
            INL_D = INL_D + DNL[j]
        INL_D = INL_D - DNL[i_d]/2.0 - DNL[0]/2.0 + INL0
        INL.append(INL_D)
    return INL

"""
crea una lista contenente le coordinate dei punti centrali di ogni bin, da
utilizzare nei plot
non è strettamente necessario, se i bin sono correttamente creati
"""
def list_central_points(bin_counts,bins_edges):
    bin_central = []
    for i,value in enumerate(bin_counts):
        central = (bins_edges[i] + bins_edges[i+1])/2.0
        bin_central.append(central)
    return bin_central

"""
funzione che, a partire da una lista di interi, riempe un istogramma, calcolando 
poi DNL ed INL. 
L'istogramma viene creato, normalizzato, escludendo i valori estremi.
I grafici vengono poi disegnati e salvati su un file, per controllo.
"""
def data_analysis(int_list,file_to_save):
    fig, axes = plt.subplots(2,2)
    bins_list = list(range(1,1024)) #da 1 a 1023 compreso
    for i,x in enumerate(bins_list): #da -0.5 a 1022.5 compresi
        bins_list[i] = x - 0.5
    bin_values, bins_edges, null= axes[1][0].hist(int_list,
                                               bins=bins_list, 
                                               density=True,
                                               alpha=0.80)
    #num_bins = len(bin_values)
    #print (np.sum(bin_values) #per vedere se l'istogramma e' normalizzato)
    DNL = calc_DNL(bin_values) #in mV
    DNL_LSB = [x/(5000.0/1024) for x in DNL] #trasformazione in LSB
    INL = calc_INL(DNL_LSB)                  #gia' in LSB (credo... ripensarci)

    print (max([abs(x) for x in DNL_LSB]))
    print (max([abs(x) for x in INL]))
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
    axes[0][1].set_title('DNL')
    axes[1][1].set_title('INL')
    fig.savefig(file_to_save, dpi=fig.dpi)

"""
funzione che gestisce l'acquisizione dei dati dall'Arduino.
La funzione prima apre la comunicazione seriale, invia il valore di prescaler opportuno, 
poi rimane in ascolto di un numero di misure decise in precedenza.
Viene anche letto il valore della frequenza di campionamento e salvato a parte.
Infine viene chiusa la porta seriale, calcolata la frequenza di campionamento 
media e restituita una lista con le misure.
"""
def arduino_sampling(prescaler,num_acquisizioni):
    num_set_letture = int(num_acquisizioni/num_samples_arduino)
    ard = serial.Serial(find_port("Arduino"),BAUDRATE,timeout=5)
    time.sleep(10)
    toSend = "%s#"%prescaler #il programma arduino richiede # come carattere finale
    try:
        ard.write(toSend)
        time.sleep(1)
        line_int = []
        freq_acq = []
        print ("acquisizione %s punti"%num_acquisizioni)
        ard.reset_input_buffer()
        ard.reset_output_buffer()
        time.sleep(0.2)
        index_set_letture = 0
        for index in range(num_acquisizioni):
            str = ard.readline().strip()
            #print ("letto %s"%str)
            if "kHz" in str:
                #print (str)
                freq_acq.append(float(str.strip("kHz\n")))
                index_set_letture = index_set_letture + 1
                print ("%s di %s - freq acquisizione: %s"%(index_set_letture,num_set_letture,str))
            else:
                line_int.append(int(str))
        if len(freq_acq) != 0:
            freq_med = sum(freq_acq)/len(freq_acq)
            print ("frequenza media di acquisizione: %fkHz"%freq_med)
        else:
            print("non misurata frequenza di acquisizione arduino")
        print ("chiusura seriale")
        ard.close()
        return line_int
    except Exception as ex:
        print (ex)
        print ("errore - chiusura seriale")
        ard.close()
        return 1

"""
funzione che gestisce l'invio dei byte via seriale alla CPLD
richiede in ingresso i parametri dell'onda triangolare
"""
def set_triangular_parameters(shp,amp,frq):
    bytes_to_send,clk_per_cnt = generate_bytes(clock_period,shp,amp,frq)
    try:
        usbToSer = serial.Serial(find_port("FTDI"),BAUDRATE,timeout=5)
        time.sleep(0.1)
        #bytes_to_send[0] = bytes(bytes_to_send[0], encoding="ascii")
        usbToSer.write(bytes_to_send[0])
        time.sleep(0.1)
        usbToSer.write(bytes_to_send[1])
        time.sleep(0.1)
        usbToSer.write(bytes_to_send[2])
        time.sleep(0.5)
        usbToSer.close()
        print ("impostata onda triangolare di frequenza %sHz"%frq)
    except Exception as ex:
        print (ex)
        usbToSer.close()
        print ("errore - chiusura seriale")
    return clk_per_cnt

def set_dimension_plot(a,b):
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = a
    fig_size[1] = b
    plt.rcParams["figure.figsize"] = fig_size

"""
funzione che itera su tutti i valori di prescaler e di frequenze, impostando
i parametri dell'onda e dell'ADC. Salva ogni set di campionamenti in un file
da utilizzare per la successiva analisi dati. Il nome del file contiene i parametri
utilizzati per la misura.
"""
def start_acq():
    set_dimension_plot(16,12)
    for p_index, presc in enumerate(prescaler_list):
        for f_index, frequency in enumerate(freq_list):
            print ("")
            print(datetime.datetime.now().time())
            print ("%sHz, presc: %s"%(frequency,presc))
            data_file = 'dati_seriale-f%s-p%s.txt'%(frequency,presc)
            graph_file = 'plot-f%s-p%s.png'%(frequency,presc)
            print (data_file)
            clk_cnt = set_triangular_parameters(shap,ampl_scelta,frequency)
            line_read_int = 1
            while(line_read_int == 1):
                line_read_int = arduino_sampling(presc,num_points_acq_list[f_index])    
            if line_read_int != 1:
                print ("salvati %s punti"%num_points_acq_list[f_index])
                data_analysis(line_read_int,graph_file)
                with open(data_file, 'w') as f:
                    for lin in line_read_int:
                        line_str = "%s\n"%lin
                        f.write(line_str)
            else:
                print ("errore nel campionamento")
    return

def main():
    start_acq()

if __name__ == "__main__":
    main()
    sys.exit(0)
