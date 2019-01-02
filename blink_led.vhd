library ieee;
use IEEE.std_logic_1164.all;
use ieee.numeric_std.all;

entity blink_led is
    generic ( 
        g_period_blink : integer := 670; --numero colpi di clock in un microsecondo: 
        g_period_blink_fast : integer := 3;
        g_time_blink_fast : integer := 20
    );
    port (
        i_reset_n : in std_logic := '1';
        i_Clk : in std_logic;
        i_data_ready : in std_logic;
        o_led_control : out std_logic;
        o_count : out integer range 670 downto 0;
        o_count_DR : out integer range 67 downto 0
    );
end entity blink_led;

architecture arch_blink of blink_led is
    signal int_data_ready : std_logic := '0';
    signal int_led_control : std_logic := '0';
    signal int_led_control_1 : std_logic := '0';
    signal int_led_control_2 : std_logic := '0';
    signal int_count : integer range 67000000 downto 0 := 0;
    signal int_count_DR : integer range 140000000 downto 0 := 0; --2 secondi con clock a 66Mhz
    signal int_count_DR_blink : integer range 67000000 downto 0 := 0; --0.3 secondi
    signal int_data_ready_flag : std_logic := '0';

    begin
        -- campionamento della linea in ingresso
        sample : process (i_Clk)
        begin
            if rising_edge(i_Clk) then
                int_data_ready <= i_data_ready;
            end if;
        end process sample;

        data_ready_flag : process(i_reset_n, i_Clk, int_count, int_data_ready_flag)
        begin
            if (i_reset_n = '0') then
                int_data_ready_flag <= '0';
            elsif (i_Clk'event and i_Clk = '1') then
                if (int_data_ready = '1') then
                    int_data_ready_flag <= '1';
                elsif(int_count_DR_blink = g_time_blink_fast) then 
                    int_data_ready_flag <= '0';
                end if;
            end if;
        end process data_ready_flag;

        counter : process(i_reset_n, i_Clk, int_count,int_data_ready_flag)
        begin
            if (i_reset_n = '0' or int_data_ready_flag = '1') then
                int_count <= 0;
            elsif (i_Clk'event and i_Clk = '1') then
                if (int_count = g_period_blink) then
                    int_count <= 0;
                else
                    int_count <= int_count + 1;
                end if;                
            end if;       
        end process counter;

        periodic_blink : process(i_reset_n, i_Clk, int_count,int_data_ready_flag)
        begin
            if (i_reset_n = '0' or int_data_ready_flag = '1') then
                int_led_control_1 <= '0';
            elsif (i_Clk'event and i_Clk = '1') then
                if (int_count = g_period_blink) then
                    int_led_control_1 <= not int_led_control_1;
                end if;                
            end if;       
        end process periodic_blink;

        counter2 : process(i_reset_n, i_Clk, int_count_DR,int_data_ready_flag)
        begin
            if (i_reset_n = '0' or int_data_ready_flag = '0') then
                int_count_DR <= 0;
            elsif (i_Clk'event and i_Clk = '1') then
                    if (int_count_DR = g_period_blink_fast) then
                        int_count_DR <= 0;
                    else
                        int_count_DR <= int_count_DR + 1;
                    end if;                
            end if;      
        end process counter2;

        fast_blink : process(i_reset_n, i_Clk, int_count, int_data_ready_flag)
        begin
            if (i_reset_n = '0' or int_data_ready_flag = '0') then
                int_led_control_2 <= '0';
                int_count_DR_blink <= 0;
            elsif (i_Clk'event and i_Clk = '1') then
                if (int_count_DR = g_period_blink_fast) then
                    int_led_control_2 <= not int_led_control_2;
                    int_count_DR_blink <= int_count_DR_blink + 1;
                    if(int_count_DR_blink = g_time_blink_fast) then
                        int_count_DR_blink <= 0;
                    end if;
                end if;          
            end if;       
        end process fast_blink;
        
        led_control : process(int_led_control_1, int_led_control_2)
        begin
            int_led_control <= int_led_control_1 or int_led_control_2;
        end process led_control;
        
        blink : process (i_Clk)
        begin
            if rising_edge(i_Clk) then
                o_led_control <= not int_led_control;
                o_count <= int_count;
                o_count_DR <= int_count_DR;
            end if;
        end process blink;
end arch_blink;