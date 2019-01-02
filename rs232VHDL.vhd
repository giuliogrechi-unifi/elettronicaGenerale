-- g_CLKS_PER_BIT = (Frequenza di i_Clk)/(Frequenza of rs232)
-- Esempio: 10 MHz Clock, 115200 baud
-- (10000000)/(115200) = 87
--66MHz clock, 57600 baudrate
--  66000000/57600 = 1146

library ieee;
use ieee.std_logic_1164.ALL;
use ieee.numeric_std.all;

entity rs232VHDL is
     generic (
        g_CLKS_PER_BIT : integer := 1156     -- Cambiare se viene cambiato il baudrate
    );
    port (
        i_reset_n     : in  std_logic := '1';
        i_Clk       : in  std_logic;
        i_RX_Serial : in  std_logic;
        o_RX_DV     : out std_logic;
        o_RX_Byte   : out std_logic_vector(7 downto 0)
    );
end rs232VHDL;

architecture rtl of rs232VHDL is

    type t_SM_Main is (s_Idle, s_RX_Start_Bit, s_RX_Data_Bits,
                       s_RX_Stop_Bit, s_Cleanup);
    signal r_SM_Main : t_SM_Main := s_Idle;
    signal r_RX_Data   : std_logic := '0';
    signal r_Clk_Count : integer range 0 to g_CLKS_PER_BIT-1 := 0;
    signal r_Bit_Index : integer range 0 to 7 := 0;  -- 8 Bits totali
    signal r_RX_Byte   : std_logic_vector(7 downto 0) := (others => '0');
    signal r_RX_DV     : std_logic := '0';

begin
    -- campionamento della linea in ingresso
    sample : process (i_Clk)
    begin
        if rising_edge(i_Clk) then
            r_RX_Data <= not i_RX_Serial;
        end if;
    end process sample;

    -- macchina a stati per la lettura rs232
    rs232_RX : process (i_Clk, i_reset_n)
    begin
        if (i_reset_n = '0') then
                r_RX_DV     <= '0';
                r_Clk_Count <= 0;
                r_Bit_Index <= 0;
                r_RX_Byte <= "00000000";
                o_RX_Byte <= "00000000";

        elsif rising_edge(i_Clk) then
            case r_SM_Main is
                when s_Idle =>
                    r_RX_DV     <= '0';
                    r_Clk_Count <= 0;
                    r_Bit_Index <= 0;
                    if r_RX_Data = '0' then       -- Start bit
                        r_SM_Main <= s_RX_Start_Bit;
                    else
                        r_SM_Main <= s_Idle;
                    end if;

                -- verifica che alla metà dei cicli di clock previsti
                -- si abbia il solito valore sulla linea
                when s_RX_Start_Bit =>
                    if r_Clk_Count = (g_CLKS_PER_BIT-1)/2 then
                        if r_RX_Data = '0' then
                            r_Clk_Count <= 0;  -- reset del contatore a metà 
                                               --del campionamento del primo bit
                            r_SM_Main <= s_RX_Data_Bits;
                        else
                            r_SM_Main <= s_Idle;
                        end if;
                    else
                        r_Clk_Count <= r_Clk_Count + 1;
                        r_SM_Main   <= s_RX_Start_Bit;
                    end if;

                -- aspetta il numero di cicli di clock previsti per leggere 
                -- un nuovo bit prosegue quando sono stati letti 8 bit
                when s_RX_Data_Bits =>
                    if r_Clk_Count < g_CLKS_PER_BIT-1 then
                        r_Clk_Count <= r_Clk_Count + 1;
                        r_SM_Main   <= s_RX_Data_Bits;
                    else
                        r_Clk_Count            <= 0;
                        r_RX_Byte(r_Bit_Index) <= r_RX_Data;
                        if r_Bit_Index < 7 then
                            r_Bit_Index <= r_Bit_Index + 1;
                            r_SM_Main   <= s_RX_Data_Bits;
                        else
                            r_Bit_Index <= 0;
                            r_SM_Main   <= s_RX_Stop_Bit;
                        end if;
                    end if;

                -- Riceve lo Stop bit.  Stop bit = 1
                when s_RX_Stop_Bit =>
                    if r_Clk_Count < g_CLKS_PER_BIT-1 then
                        r_Clk_Count <= r_Clk_Count + 1;
                        r_SM_Main   <= s_RX_Stop_Bit;
                    else
                        r_RX_DV     <= '1';
                        r_Clk_Count <= 0;
                        r_SM_Main   <= s_Cleanup;
                        --così dovrebbe aggiornare l'uscita
                        --solo quando finisce un ciclo di lettura
                        --e non ogni volta che salva un bit
                        o_RX_Byte <= r_RX_Byte;
                    end if;

                 -- Rimane in stato di cleanup per un ciclo di clock (e resetta il dataReady)
                when s_Cleanup =>
                    r_SM_Main <= s_Idle;
                    r_RX_DV   <= '0';

                when others =>
                    r_SM_Main <= s_Idle;

                end case;
        end if;
    end process rs232_RX;

    ff_DV: process(i_Clk)
    begin
        if rising_edge(i_Clk) then
                o_RX_DV   <= r_RX_DV;
        end if;
    end process ff_DV;
--o_RX_Byte <= r_RX_Byte;

end rtl;