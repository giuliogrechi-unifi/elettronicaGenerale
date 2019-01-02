library ieee;
use ieee.std_logic_1164.ALL;
use ieee.numeric_std.all;

entity project is
    port (
        --gli input sono il clock, la seriale ed
        --un reset fisico
        i_Clk : in std_logic;
        i_RX_Serial : in std_logic;
        i_reset : in std_logic;
        --gli output sono i 12 bit da inviare 
        --al DAC ed i 3 bit da inviare al
        --binary DAC. eventualmente da aggiungere
        --un debug sui led
        o_ToDAC : out std_logic_vector (11 downto 0);
        o_ToBIN : out std_logic_vector (2 downto 0);
        o_wr : out std_logic; --segnale di write
        o_square : out std_logic; -- onda quadra =1 quando il conteggio decresce
        --o_value_read : out std_logic_vector (11 downto 0);
        o_led : out std_logic;
        o_led_2 : out std_logic;
        o_led_3 : out std_logic
    );
end project;

architecture archit of project is
    component rs232VHDL is
        port (
            i_reset_n     : in  std_logic := '1';
            i_Clk       : in  std_logic;
            i_RX_Serial : in  std_logic;
            o_RX_DV     : out std_logic;
            o_RX_Byte   : out std_logic_vector(7 downto 0)
            );
    end component;
    component decode_serial_data is 
        port (
            i_Clk       : in  std_logic;
            i_RX_Byte   : in  std_logic_vector (7 downto 0);
            i_RX_DV     : in std_logic;
            o_D_Ready   : out std_logic;
            o_Contr_shap: out std_logic_vector (1 downto 0);
            o_Ampl      : out std_logic_vector (2 downto 0);
            o_RX_Byte1  : out  std_logic_vector (7 downto 0);
            o_RX_Byte2  : out  std_logic_vector (7 downto 0);
            o_RX_Byte3  : out  std_logic_vector (7 downto 0);
            o_Value  : out  std_logic_vector (11 downto 0);
            o_B1_Ready   : out std_logic;
            o_B2_Ready   : out std_logic
            );
    end component decode_serial_data;
    component user_countSM is
        port (
            i_count_max_1 : in integer range 4095 downto 0 := 66;
            i_reset_n_sinc : in std_logic := '1';
            i_Clk : in std_logic;
            i_tri : in std_logic_vector (1 downto 0) := "11";
            i_const_value : in integer range 4095 downto 0; --12 bit per l'eventuale valore costante
            o_wr : out std_logic;
            o_falling : out std_logic;
            count_1 : out std_logic_vector (6 downto 0); --fino a 128
            count_2 : out std_logic_vector (11 downto 0) --fino a 4096
        );
    end component user_countSM;
--da cambiare, questi valori vanno bene per la simulazione
    component blink_led is
        generic ( 
            g_period_blink : integer := 670; --numero colpi di clock in un microsecondo: 
            g_period_blink_fast : integer := 3;
            g_time_blink_fast : integer := 20
            --i_count_max_1 : integer := 100
        );
        port (
            i_reset_n : in std_logic := '1';
            i_Clk : in std_logic;
            i_data_ready : in std_logic;
            o_led_control : out std_logic
            --o_count : out integer range 670 downto 0;
            --o_count_DR : out integer range 67 downto 0
        );
    end component blink_led;
        --SEGNALI INTERNi
        --
        --
    signal int_Serial_Ready : std_logic;
    signal int_Serial_Byte_Read : std_logic_vector(7 downto 0);
    signal int_Data_Decoded_Ready : std_logic;
    signal int_Shape : std_logic_vector (1 downto 0);
    signal int_DAC_value : std_logic_vector(11 downto 0);
    signal output_reset : std_logic;
    signal int_led2 : std_logic;
    signal int_led3 : std_logic;

    begin
        outReset : process (i_Clk, i_reset_n, int_Data_Decoded_Ready)
        begin
            if rising_edge(i_Clk) then
                output_reset <= i_reset_n or (NOT int_Data_Decoded_Ready);
            end if;
        end process outReset;
        
        toggle_led : process (i_Clk, i_reset_n, int_Data_Decoded_Ready)
        begin
            if rising_edge(i_Clk) then
                o_led_2 <= not int_led2;
                o_led_3 <= not int_led3;
            end if;
        end process toggle_led;

        rs232 : rs232VHDL port map(
            i_reset_n     => i_reset_n,
            i_Clk       => i_Clk,
            i_RX_Serial => i_RX_Serial,
            o_RX_DV     => int_Serial_Ready,
            o_RX_Byte   => int_Serial_Byte_Read
            );

        decode : decode_serial_data port map(
            i_Clk       => i_Clk,
            i_RX_Byte   => int_Serial_Byte_Read,
            i_RX_DV     => int_Serial_Ready,
            o_D_Ready   => int_Data_Decoded_Ready,
            o_Contr_shap => int_Shape,
            o_Ampl      => o_ToBIN, --ampiezza -> binary resistor
            --o_RX_Byte1  : out  std_logic_vector (7 downto 0);
            --o_RX_Byte2  : out  std_logic_vector (7 downto 0);
            --o_RX_Byte3  : out  std_logic_vector (7 downto 0);
            o_Value  => int_DAC_value,
            o_B1_Ready => int_led2,
            o_B2_Ready => int_led3
        );

        output_generate : user_countSM port map(
            i_count_max_1 => to_integer(unsigned(int_DAC_value)),
            i_reset_n_sinc => output_reset,
            i_Clk => i_Clk,
            i_tri => int_Shape,
            i_const_value => to_integer(unsigned(int_DAC_value)),
            --outputs
            o_wr => o_wr,
            o_falling => o_square,
            --count_1 : out std_logic_vector (6 downto 0); --fino a 128
            count_2 => o_ToDAC
        );

        blink : blink_led
        generic map(
                g_period_blink => 67, --numero colpi di clock toggle led idle: 
                g_period_blink_fast => 30, -- numero colpi di clock durata singolo blink data ready
                g_time_blink_fast => 130 --numero colpi di clock durata totale data ready
        )
            port map (
                i_reset_n => i_reset_n,
                i_Clk => i_Clk,
                i_data_ready => int_Data_Decoded_Ready,
                o_led_control => o_led
                --o_count : out integer range 670 downto 0;
                --o_count_DR : out integer range 67 downto 0
            );
    
        --o_value_read <= int_DAC_value; --debig

    end archit;