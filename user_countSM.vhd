library ieee;
use IEEE.std_logic_1164.all;
--use ieee.std_logic_arith.all; --deprecato
--use ieee.std_logic_unsigned.all; --deprecato
use ieee.numeric_std.all;

entity user_countSM is
    generic ( 
        g_count_max_2 : integer := 4095
        --g_count_max_2 : integer := 100
    );
    port (
        i_count_max_1 : in integer range 4095 downto 0 := 66; --66 per 1MHz
        i_reset_n_sinc : in std_logic := '1';
        i_Clk : in std_logic;
        i_tri : in std_logic_vector (1 downto 0) :="11";
        i_const_value : in integer range 4095 downto 0; --12 bit per l'eventuale valore costante
        o_wr : out std_logic;
        o_falling : out std_logic;
        count_1 : out std_logic_vector (6 downto 0); --fino a 128
        count_2 : out std_logic_vector (11 downto 0) --fino a 4096
    );
end entity user_countSM;

architecture archSM of user_countSM is
    signal int_count_1: integer range 4095 downto 0 := 0;
    signal int_count_2: integer range 4095 downto 0 := 0;
    signal int_count_max_1 : integer range 4095 downto 0;
    signal int_const_value : integer range 4095 downto 0; --12 bit per l'eventuale valore costante
    signal int_wr: std_logic;
    signal int_falling : std_logic;
    signal int_reset_n_sinc : std_logic;
    --signal int_value : integer range 4095 downto 0;
    type t_SM_Main is (s_Idle, s_write_constant, s_count_up,
                       s_count_down, s_square, s_reset);
    signal r_SM_Main, r_SM_Next : t_SM_Main := s_Idle; 
    signal int_tri : std_logic_vector (1 downto 0) :="11";
    signal int_wr_flag : std_logic := '0';
    
    begin
        sample : process (i_Clk) --campionamento degli ingressi con dei flip flop
        begin
            if rising_edge(i_Clk) then
                int_count_max_1 <= i_count_max_1;
                int_tri <= i_tri;
                int_const_value  <= i_const_value;
                int_reset_n_sinc <= int_reset_n_sinc;
            end if;
        end process sample;

        ffOUT : process(i_Clk)
        begin 
            if rising_edge(i_Clk) then
                count_1 <= std_logic_vector(to_unsigned(int_count_1, 
                                                        count_1'length));
                count_2 <= std_logic_vector(to_unsigned(int_count_2, 
                                                        count_2'length));
                o_wr <= NOT int_wr;
                o_falling <= int_falling;
            end if;
        end process ffOUT;

        write : process (i_Clk)
        begin
            if rising_edge(i_Clk) then
                if (int_reset_n_sinc /= '0') then
                    if (int_count_1 = 1 
                                    or int_count_1 = 2 
                                    or int_count_1 = 3) then
                        int_wr <= '1';
                    else
                        int_wr <= '0';                   
                    end if;
                else
                    int_wr <= '0';
                end if;
            end if;
        end process write;

        state_machine : process(i_Clk)
        begin
            if rising_edge(i_Clk) then
                r_SM_Main <= r_SM_Next;
            end if;
            if rising_edge(i_Clk) then
                case r_SM_Main is
                    when s_Idle =>
                    if (int_reset_n_sinc /= '0') then
                        int_count_1 <= 0;
                        int_count_2 <= 0;
                        int_falling <= '0';
                        if (int_tri = "00") then
                            r_SM_Next <= s_Idle;
                        elsif(int_tri = "01") then
                            r_SM_Next <= s_write_constant;
                        elsif(int_tri = "10") then
                            r_SM_Next <= s_square;
                        elsif(int_tri = "11") then
                            r_SM_Next <= s_count_up;
                        end if;
                    else 
                            r_SM_Next <= s_reset;     
                    end if;

                    when s_write_constant =>
                        if (int_reset_n_sinc /= '0') then
                            int_count_2 <= int_const_value;
                            int_count_1 <= int_count_1 + 1;
                                if(int_count_1 = 1000) then --per generare il segnale di write
                                    int_count_1 <= 0;
                                end if;
                                if(int_tri /= "01") then
                                    r_SM_Next <= s_Idle;
                                else 
                                    r_SM_Next <= s_write_constant;
                                end if;
                        else 
                            r_SM_Next <= s_reset;     
                        end if;

                    when s_count_up =>
                        if (int_reset_n_sinc /= '0') then
                            if(int_count_1 = int_count_max_1) then
                                int_count_1 <= 0;
                                int_count_2 <= int_count_2 + 1;
                                if(int_tri /= "11") then
                                    r_SM_Next <= s_Idle;
                                elsif(int_tri = "11" and int_count_2 
                                                    < g_count_max_2 - 1) then
                                    r_SM_Next <= s_count_up;
                                    int_falling  <= '0';
                                elsif(int_tri = "11" and int_count_2 
                                                    = g_count_max_2 - 1) then
                                    r_SM_Next <= s_count_down;
                                end if;
                            else 
                                int_count_1 <= int_count_1 + 1;
                            end if;
                        else 
                                r_SM_Next <= s_reset;     
                        end if;

                    when s_count_down =>
                        if (int_reset_n_sinc /= '0') then
                            if(int_count_1 = int_count_max_1) then
                                int_count_1 <= 0;
                                int_count_2 <= int_count_2 - 1;
                                if(int_tri /= "11") then
                                    r_SM_Next <= s_Idle;
                                elsif(int_tri = "11" and int_count_2 > 1) then
                                    r_SM_Next <= s_count_down;
                                    int_falling  <= '1';
                                elsif(int_tri = "11" and int_count_2 = 1) then
                                    r_SM_Next <= s_count_up;
                                end if;
                            else 
                                int_count_1 <= int_count_1 + 1;
                            end if;
                        else 
                            r_SM_Next <= s_reset;                        
                        end if;

                    when s_square => 
                        if (int_reset_n_sinc /= '0') then
                            int_count_1 <= int_count_1 + 1;
                            if(int_count_1 = int_count_max_1 - 1) then --per generare il segnale di write
                                int_count_1 <= 0;
                                if(int_count_2 = 4095) then
                                    int_count_2 <= 0;
                                else
                                    int_count_2 <= 4095;
                                end if;
                            end if;
                            if(int_tri /= "10") then
                                r_SM_Next <= s_Idle;
                            else
                                r_SM_Next <= s_square;
                            end if;
                        else 
                            r_SM_Next <= s_reset;
                        end if;

                    when s_reset =>
                        int_count_1 <= 0;
                        int_count_2 <= 0;
                        int_falling <= '0';
                        r_SM_Next <= s_Idle;

                end case;
            end if;
        end process state_machine;
  end archSM;
