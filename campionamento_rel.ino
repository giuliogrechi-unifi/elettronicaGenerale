/*
NUM_SAMPLES e' il numero di campioni consecutivi acquisiti dall'ADC. E' limitato
dalla quantita' di memoria ram a disposizione (8kbyte per ATMega2560);
raggiunto il numero preimpostato (campionando una volta ogni 13 colpi di clock ADC),
la conversione si arresta,i dati vengono inviati via seriale e poi ricomincia.
La trasmissione seriale dei dati impiega molto piu' del tempo impiegato per campionarli.
PRE_SCALER e' il coefficiente di divisione del clock di sistema (16MHz) per formare
il clock ADC. Il datasheet richiede una frequenza di clock ADC tra 50 e 200 KHz per avere
piena risoluzione. I valori possibili per il
prescaler sono 8 16 32 64 128 (implementati qui da 32 a 128).
Il sampling rate e' [ADC clock] / [prescaler] / [conversion clock cycles]
*/
#define NUM_SAMPLES 3500
int PRE_SCALER;
volatile int numSamples=0;
int i = 0;
long t, t0;
volatile short values[NUM_SAMPLES];
void setup(){
    pinMode(13,OUTPUT);               //per debug
    digitalWrite(13,0);
    Serial.begin(57600);
    while(Serial.available() == 0){}  //idle fino a quando non riceve un segnale
    while (Serial.available() > 0) {
        PRE_SCALER = Serial.parseInt();
        //il programma richiede che il numero di prescaler inviato venga
        //seguito dal carattere di controllo #
        //se e' vero il led sulla scheda si accende per segnalare la corretta
        //ricezione di un valore
        if (Serial.read() == '#') {
            if(PRE_SCALER == 128 or PRE_SCALER == 64 or PRE_SCALER == 32){
                digitalWrite(13,1);
            }
        }
    }
    DIDR0  |= 0xFF;         //disattiva input digitali sui pin analogici
    DIDR2  |= 0xFF;         //solo per atmega2560
    ADCSRA = 0;             // clear registro ADCSRA da impostazioni di default
    ADCSRB = 0;             // clear registro ADCSRB da impostazioni di default
    ADMUX |= (0 & 0x07);    // imposta tutti 0 nel registro ADMUX
    ADMUX |= (1 << REFS0);  // impostazione reference voltage alimentazione
    //ADMUX |= (1 << REFS1);  // impostazione reference voltage interna 1.1v

    switch(PRE_SCALER){
        case 128:
            ADCSRA |= (1 << ADPS2) | (1 << ADPS1) | (1 << ADPS0); // 128 prescaler per 9.61 KHz
            break;
        case 64:
            ADCSRA |= (1 << ADPS2) | (1 << ADPS1);                // 64 prescaler per 19.2 KHz
            break;
        case 32:
            ADCSRA |= (1 << ADPS2) | (1 << ADPS0);                // 32 prescaler per 38.5 KHz
            break;
    }
    ADCSRA |= (1 << ADATE); // enable auto trigger
    ADCSRA |= (1 << ADIE);  // enable interrupts a conversione ADC completata
    ADCSRA |= (1 << ADEN);  // enable ADC
    ADCSRB |= (0 & 0x07);   // free running mode, trigger automatico al completamento della conversione
    delay(1000);            // per garantire stabilità e per sincronizzare con il software lato pc
    ADCSRA |= (1 << ADSC);  // start ADC
    t0 = micros();          // per misurare il tempo trascorso
}

/*
funzione chiamata quando si verifica l'interrupt alla fine della conversione ADC
I valori convertiti si trovano all'interno dei due registri ad 8 bit:
ADCL, che contiene gli 8 bit meno significativi, e ADCH, che contiene i due
bit piu' significativi. Il risultato finale viene salvato in una variabile di tipo
short dopo aver costruito il numero a 16 bit utilizzando un'operazione di bit shift.
Infine viene incrementato un contatore del numero di eventi.
*/
ISR(ADC_vect){
    volatile short y = ADCL;  // legge 8 bit dall' ADC (8 meno significativi)
    volatile short x = ADCH;  // legge 8 bit dall' ADC (2 piu' significativi)
    y |= x << 8;
    values[numSamples] = y;
    numSamples++;
}

/*
Il loop contiene solo una condizione if, che si verifica quando e' stato campionato
un numero di segnali pari al massimo definito da NUM_SAMPLES.
All'interno della condizione if viene per prima cosa calcolato il tempo trascorso,
disabilitato l'ADC e la generazione degli interrupt, poi inviati via seriale i dati,
infine calcolata e inviata sempre via seriale la frequenza di campionamento (seguita
dalla stringa kHz).
Un ritardo variabile in modo casuale permette di rendere casuale il punto
da cui comincia il successivo campionamento.
*/
void loop(){
    if (numSamples>NUM_SAMPLES){
        t = micros()-t0;         // tempo trascorso
        ADCSRA &= ~(1 << ADATE); // disable auto trigger
        ADCSRA &= ~(1 << ADIE);  // disable interrupts a conversione ADC completata
        ADCSRA &= ~(1 << ADSC);  // stop ADC
        for (int j=0;j<NUM_SAMPLES;j++){
            Serial.println(values[j]);
        }
        Serial.print((float)1000*NUM_SAMPLES/t);
        Serial.println("kHz");
        delay(random(0,500));
        // restart
        ADCSRA |= (1 << ADATE); // enable auto trigger
        ADCSRA |= (1 << ADIE);  // enable interrupts a conversione ADC completata
        ADCSRA |= (1 << ADSC);  // start ADC
        t0 = micros();
        numSamples=0;
    }
}