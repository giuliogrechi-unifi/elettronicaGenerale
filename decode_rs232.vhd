library ieee;
use ieee.std_logic_1164.ALL;
use ieee.numeric_std.all;

entity decode_serial_data is
    port (
      i_Clk : in  std_logic;
      i_RX_Byte : in  std_logic_vector (7 downto 0);
      i_RX_DV : in std_logic;
      i_reset_n     : in  std_logic := '1';
      o_D_Ready : out std_logic;
      o_Contr_shap : out std_logic_vector (1 downto 0);
      o_Ampl : out std_logic_vector (2 downto 0);
      o_RX_Byte1 : out  std_logic_vector (7 downto 0);
      o_RX_Byte2 : out  std_logic_vector (7 downto 0);
      o_RX_Byte3 : out  std_logic_vector (7 downto 0);
      o_Value : out  std_logic_vector (11 downto 0);
      o_B1_Ready : out std_logic;
      o_B2_Ready : out std_logic
      );
end decode_serial_data;

architecture dcd of decode_serial_data is
    signal int_Rx_DV : std_logic;
    signal int_Byte1 : std_logic_vector (7 downto 0);
    signal int_Byte2 : std_logic_vector (7 downto 0);
    signal int_Byte3 : std_logic_vector (7 downto 0);
    signal int_D_ready : std_logic;
    signal int_Ampl: std_logic_vector (2 downto 0); --3 bit al massimo
    signal int_Contr_shap : std_logic_vector (1 downto 0);
    signal int_Byte1_read: std_logic;
    signal int_Byte2_read: std_logic;
    signal int_Value : std_logic_vector (11 downto 0);
    type t_SM_Main is (s_Idle, s_readByte1, s_readByte2, s_readByte3,
                       s_cleanup);
    signal r_SM_Main : t_SM_Main := s_Idle;

    begin
        p_sample : process(i_Clk)
        begin
            if rising_edge(i_Clk) then
                int_Rx_DV <= i_RX_DV;
            end if;
        end process p_sample;

        state_machine : process(i_Clk,i_reset_n)
        begin
            if (i_reset_n = '0') then
                int_D_ready <= '0';
                int_Byte1_read <= '0';
                int_Byte2_read <= '0';
                int_Ampl <= "000";
                int_Byte1 <= "00000000";
                int_Byte2 <= "00000000";
                int_Byte3 <= "00000000";
            elsif rising_edge(i_Clk) then
                case r_SM_Main is
                    when s_Idle =>
                        int_D_ready <= '0';
                        int_Byte1_read <= '0';
                        int_Byte2_read <= '0';
                        int_Ampl <= "000";
                        int_Byte1 <= "00000000";
                        int_Byte2 <= "00000000";
                        int_Byte3 <= "00000000";
                        if (int_Rx_DV = '1') then
                            r_SM_Main <= s_readByte1;
                        else
                            r_SM_Main <= s_Idle;
                        end if;

                    when s_readByte1 =>
                        if (int_Rx_DV = '0') then --legge solo quando il dataReady tornato a 0, in teoria al successivo colpo di clock questo dovrebbe essere sempre vero
                            int_Byte1 <= i_RX_Byte;
                            int_Byte1_read <= '1';
                            int_Byte2_read <= '0';
                            r_SM_Main <= s_readByte2;
                        else
                            r_SM_Main <= s_readByte1;
                        end if;
                        
                    when s_readByte2 =>
                        if (int_Rx_DV = '1') then
                            int_Byte2 <= i_RX_Byte;
                            int_Byte1_read <= '1';
                            int_Byte2_read <= '1';
                            r_SM_Main <= s_readByte3;
                        else
                            r_SM_Main <= s_readByte2;
                        end if;

                    when s_readByte3 =>
                        if (int_Rx_DV = '1') then
                            int_Byte3 <= i_RX_Byte;
                            r_SM_Main <= s_cleanup;
                            int_D_Ready <= '1';
                            int_Ampl <= int_Byte1(2 downto 0);
                            int_Contr_shap <= int_Byte1(5 downto 4);
                            --int_Value(3 downto 0) <= ("0000");
                            --int_Value(11 downto 4) <= i_RX_Byte(7 downto 0);
                            int_Value(5 downto 0) <= int_Byte2(5 downto 0);
                            int_Value(11 downto 6) <= i_RX_Byte(5 downto 0);
                        else
                            r_SM_Main <= s_readByte3;
                        end if;

                    when s_cleanup =>
                        o_RX_Byte1 <= int_Byte1;
                        o_RX_Byte2 <= int_Byte2;
                        o_RX_Byte3 <= int_Byte3;
                        o_Value <= int_value;
                        o_Contr_shap <= int_Contr_shap;
                        o_Ampl <= int_Ampl;
                        int_D_Ready <= '0';
                        r_SM_Main <= s_Idle;
                end case;
            end if;
        end process state_machine;

        ff_D_ready: process(i_Clk)
        begin
            if rising_edge(i_Clk) then
                o_D_Ready <= int_D_Ready;
                o_B1_Ready <= int_Byte1_read;
                o_B2_Ready <= int_Byte2_read;
            end if;
        end process ff_D_ready;
end dcd;
