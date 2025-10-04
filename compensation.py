import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from smps import *
import control as ct
import numpy as np

if __name__ == '__main__':

    # --- Main Window Setup ---
    # Create the main application window
    root = tk.Tk()
    root.title("Data Entry and Computation with Plot")
    root.geometry("1000x850") # Increased window size to accommodate the plot

    # --- Main Layout Frames ---
    control_frame = tk.Frame(root, padx=10, pady=10)
    control_frame.pack(side=tk.LEFT, fill=tk.Y)

    plot_frame = tk.Frame(root, padx=10, pady=10)
    plot_frame.pack(side=tk.RIGHT, expand=True, fill=tk.BOTH)

    # ----- Control Widgets (Left Frame) -----
    grid_row = 0

    # --- Dropdown Menu for Converter Type Units ---
    tk.Label(control_frame, text="Converter Type:").grid(row=grid_row, column=0, sticky="w", pady=5)
    converter_type_var = tk.StringVar(root)
    converter_options = ['buck, voltage mode', 'buck, current mode', 'boost, voltage mode', 'boost, current mode']
    converter_type_var.set(converter_options[0]) # Set default value
    conv_type_dropdown = tk.OptionMenu(control_frame, converter_type_var, *converter_options)
    conv_type_dropdown.grid(row=grid_row, column=1, sticky="ew", pady=5, padx=5)
    grid_row += 1

    tk.Label(control_frame, text="Switching Freq (f_sw):").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_fsw = tk.Entry(control_frame, width=30)
    entry_fsw.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_fsw.insert(0, "100e3") # Add default value for a simple system
    grid_row += 1

    tk.Label(control_frame, text="Input Voltage (Vin):").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_vin = tk.Entry(control_frame, width=30)
    entry_vin.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_vin.insert(0, "12.0") # Add default value for a simple system
    grid_row += 1
    
    tk.Label(control_frame, text="Output Voltage (Vout):").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_vout = tk.Entry(control_frame, width=30)
    entry_vout.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_vout.insert(0, "5.0") # Add default value for a simple system
    grid_row += 1
    
    tk.Label(control_frame, text="Typical Output Current (Iout):").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_iout = tk.Entry(control_frame, width=30)
    entry_iout.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_iout.insert(0, "1.0") # Add default value for a simple system
    grid_row += 1
    
    tk.Label(control_frame, text="Target/Approx. Efficiency:").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_eff = tk.Entry(control_frame, width=30)
    entry_eff.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_eff.insert(0, "0.90") # Add default value for a simple system
    grid_row += 1
    
    tk.Label(control_frame, text="Output Inductor (Lout):").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_lout = tk.Entry(control_frame, width=30)
    entry_lout.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_lout.insert(0, "10.0e-6") # Add default value for a simple system
    grid_row += 1
    
    tk.Label(control_frame, text="Output Capacitor (Cout):").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_cout = tk.Entry(control_frame, width=30)
    entry_cout.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_cout.insert(0, "100.0e-6") # Add default value for a simple system
    grid_row += 1
    
    tk.Label(control_frame, text="Output Cap ESR (R_cout_esr):").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_cout_esr = tk.Entry(control_frame, width=30)
    entry_cout_esr.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_cout_esr.insert(0, "0.1") # Add default value for a simple system
    grid_row += 1
    
    tk.Label(control_frame, text="Inductor Current Sense R:").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_rsns = tk.Entry(control_frame, width=30)
    entry_rsns.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_rsns.insert(0, "0.05") # Add default value for a simple system
    grid_row += 1
    
    tk.Label(control_frame, text="I_sns Gain:").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_isns_gain = tk.Entry(control_frame, width=30)
    entry_isns_gain.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_isns_gain.insert(0, "20.0") # Add default value for a simple system
    grid_row += 1

    tk.Label(control_frame, text="Vslope:").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_vslope = tk.Entry(control_frame, width=30)
    entry_vslope.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_vslope.insert(0, "5.0") # Add default value for a simple system
    grid_row += 1

    # --- SEPARATOR ---
    tk.Frame(control_frame, height=2, bg="grey", width=400).grid(row=grid_row, columnspan=2, pady=10)
    grid_row += 1

    # --- Dropdown Menu for Compensation Type ---
    tk.Label(control_frame, text="Comp Type:").grid(row=grid_row, column=0, sticky="w", pady=5)
    comp_type_var = tk.StringVar(root)
    comp_options = ['Type 1', 'Type 2', 'Type 3']
    comp_type_var.set(comp_options[0]) # Set default value
    comp_type_dropdown = tk.OptionMenu(control_frame, comp_type_var, *comp_options)
    comp_type_dropdown.grid(row=grid_row, column=1, sticky="ew", pady=5, padx=5)
    grid_row += 1

    # --- Compensation Value Inputs ---
    tk.Label(control_frame, text="Rfbt:").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_rfbt = tk.Entry(control_frame, width=30)
    entry_rfbt.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_rfbt.insert(0, "10e3")
    grid_row += 1

    tk.Label(control_frame, text="Rcomp:").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_rcomp = tk.Entry(control_frame, width=30)
    entry_rcomp.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_rcomp.insert(0, "10e3")
    grid_row += 1

    tk.Label(control_frame, text="Ccomp:").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_ccomp = tk.Entry(control_frame, width=30)
    entry_ccomp.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_ccomp.insert(0, "10e-9")
    grid_row += 1

    tk.Label(control_frame, text="Chf:").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_chf = tk.Entry(control_frame, width=30)
    entry_chf.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_chf.insert(0, "100e-12")
    grid_row += 1

    tk.Label(control_frame, text="Cff:").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_cff = tk.Entry(control_frame, width=30)
    entry_cff.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_cff.insert(0, "0")
    grid_row += 1

    tk.Label(control_frame, text="Rff:").grid(row=grid_row, column=0, sticky="w", pady=5)
    entry_rff = tk.Entry(control_frame, width=30)
    entry_rff.grid(row=grid_row, column=1, pady=5, padx=5)
    entry_rff.insert(0, "0")
    grid_row += 1

    # --- SEPARATOR ---
    tk.Frame(control_frame, height=2, bg="grey", width=400).grid(row=grid_row, columnspan=2, pady=10)
    grid_row += 1

    # --- Result Display ---
    tk.Label(control_frame, text="Steady State Duty Cycle:").grid(row=grid_row, column=0, sticky="w", pady=5)
    duty_label = tk.Label(control_frame, text="---", wraplength=300, justify="left")
    duty_label.grid(row=grid_row, column=1, pady=5, sticky="w")
    grid_row += 1

    tk.Label(control_frame, text="Ind. Ripple Current Pk-Pk:").grid(row=grid_row, column=0, sticky="w", pady=5)
    ind_rip_label = tk.Label(control_frame, text="--- A", wraplength=300, justify="left")
    ind_rip_label.grid(row=grid_row, column=1, pady=5, sticky="w")
    grid_row += 1

    tk.Label(control_frame, text="Peak Inductor Current:").grid(row=grid_row, column=0, sticky="w", pady=5)
    ind_peak_label = tk.Label(control_frame, text="--- A", wraplength=300, justify="left")
    ind_peak_label.grid(row=grid_row, column=1, pady=5, sticky="w")
    grid_row += 1
    
    # --- Perform Computation function ---
    def perform_computation():
        """
        This function is executed when the 'COMPUTE' button is pressed.
        It creates a transfer function from the user's input and generates
        a Bode plot on a single chart with two y-axes.
        """
        try:
            # basic power stage computations

            # --- Get and Parse Transfer Function Coefficients ---
            fsw = float(entry_fsw.get())
            vin = float(entry_vin.get())
            vout = float(entry_vout.get())
            iout = float(entry_iout.get())
            lout = float(entry_lout.get())
            cout = float(entry_cout.get())
            cout_esr = float(entry_cout_esr.get())
            rsns = float(entry_rsns.get())
            isns_gain = float(entry_isns_gain.get())
            vslope = float(entry_vslope.get())
            eff = float(entry_eff.get())
            conv_type_str = converter_type_var.get()
            if 'buck' in conv_type_str:
                if 'current mode' in conv_type_str:
                    conv_type = SmpsType.BUCK_CURRENT_MODE
                else:
                    conv_type = SmpsType.BUCK_VOLTAGE_MODE
            else:
                if 'current mode' in conv_type_str:
                    conv_type = SmpsType.BOOST_CURRENT_MODE
                else:
                    conv_type = SmpsType.BOOST_VOLTAGE_MODE
            
            info = SmpsInfo(conv_type, fsw, vin, vout, iout, lout, cout, cout_esr, rsns, isns_gain, vslope, eff)

            # --- basic power stage calcs ---
            duty = calc_duty_cycle(info)
            I_rip = calc_ind_ripple(info)
            I_max = calc_ind_ipeak(info)
            duty_label.config(text=f'{duty:.3f}')
            ind_rip_label.config(text=f'{I_rip:.2f} A')
            ind_peak_label.config(text=f'{I_max:.2f} A')

            # --- Create Control System Model ---
            if conv_type == SmpsType.BUCK_VOLTAGE_MODE:
                power_stage_xfer = buck_vmode(info)
            elif conv_type == SmpsType.BUCK_CURRENT_MODE:
                power_stage_xfer = buck_imode(info)
            elif conv_type == SmpsType.BOOST_VOLTAGE_MODE:
                power_stage_xfer = boost_imode(info)
            elif conv_type == SmpsType.BOOST_CURRENT_MODE:
                power_stage_xfer = boost_imode(info)

            # --- Compensation ---
            Rfbt = float(entry_rfbt.get())
            Ccomp = float(entry_ccomp.get())
            Rcomp = float(entry_rcomp.get())
            Chf = float(entry_chf.get())
            Rff = float(entry_rff.get())
            Cff = float(entry_cff.get())

            comp_type = comp_type_var.get()
            if comp_type == 'Type 1':
                comp_xfer = error_amp_type1(Rfbt, Ccomp)
            elif comp_type == 'Type 2':
                comp_xfer = error_amp_type2(Rfbt, Rcomp, Ccomp, Chf)
            elif comp_type == 'Type 3':
                comp_xfer = error_amp_type3(Rfbt, Rcomp, Ccomp, Chf, Rff, Cff)

            G = power_stage_xfer * comp_xfer #   
            # --- Generate Bode Plot Data ---
            w = np.logspace(1, 7, 500)
            # Call control.bode to get the data without plotting it automatically
            mag, phase_rad, omega = ct.bode(G, w, plot=False, Hz=True, wrap_phase=True)

            # Convert magnitude to decibels (dB) and phase to degrees
            mag_db = 20 * np.log10(mag)
            phase_deg = np.rad2deg(phase_rad)

            # --- Update the Plot ---
            # Clear the entire figure to ensure a clean slate for the twin axes
            fig.clear()

            # Create the primary axis for the magnitude plot
            ax1 = fig.add_subplot(111)

            # Create a secondary axis that shares the same x-axis for the phase plot
            ax2 = ax1.twinx()

            # Plot Magnitude on the left axis (ax1)
            ax1.semilogx(omega, mag_db, color='royalblue', label='Magnitude (dB)')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel('Magnitude (dB)', color='royalblue')
            ax1.set_ylim(-60, 60)
            ax1.tick_params(axis='y', labelcolor='royalblue')

            # Plot Phase on the right axis (ax2)
            ax2.semilogx(omega, phase_deg, color='orangered', label='Phase (deg)')
            ax2.set_ylabel('Phase (deg)', color='orangered')
            ax2.set_ylim(-90, 90)
            ax2.tick_params(axis='y', labelcolor='orangered')

            # Add grid and title
            ax1.grid(True, which='both', linestyle=':')
            ax1.set_title('Bode Plot')

            # Adjust layout to prevent labels from overlapping
            fig.tight_layout()
            canvas.draw()

        except ValueError:
            messagebox.showerror("ValueError parsing string to float, probably.", "Please enter valid numbers.")

    # --- COMPUTE Button ---
    compute_button = tk.Button(control_frame, text="COMPUTE", command=perform_computation, width=15)
    compute_button.grid(row=grid_row, column=0, columnspan=2, pady=15)
    grid_row += 1
    
    # --- Matplotlib Plot Setup (Right Frame) ---
    fig = Figure(figsize=(6, 5), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_title("Enter Transfer Function to Generate Bode Plot")
    ax.set_xlabel("Frequency (rad/s)")
    ax.grid(True, which='both', linestyle=':')
    fig.tight_layout()

    # Create a Tkinter canvas to hold the Matplotlib figure
    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    # --- Start the Application ---
    root.mainloop()