import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import random
import time

# Variables to control the simulation
running = False
variable = 10
iterations = []
values = []

# Function to randomly adjust the variable
def random_adjust(variable):
    adjustment = random.choice([-1, 1])
    return variable + adjustment

# Function to start the simulation
def start_simulation():
    global running
    running = True
    run_simulation()

# Function to pause the simulation
def pause_simulation():
    global running
    running = False

# Function to reset the simulation
def reset_simulation():
    global variable, iterations, values, running
    running = False
    variable = 10
    iterations.clear()
    values.clear()
    update_plot()

# Main simulation function
def run_simulation():
    global variable
    iteration = 0
    while running:
        variable = random_adjust(variable)
        iterations.append(iteration)
        values.append(variable)
        update_plot()
        iteration += 1
        time.sleep(0.5)  # Slow down the simulation loop
        window.update()  # Update tkinter window to reflect changes

# Function to update the matplotlib plot
def update_plot():
    ax.clear()
    ax.plot(iterations, values, marker='o', linestyle='-')
    ax.set_title("Simulation Control with Tkinter")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Value")
    canvas.draw()

# Set up the tkinter window
window = tk.Tk()
window.title("Simulation Control Panel")

# Set up matplotlib figure and axis
fig, ax = plt.subplots(figsize=(5, 4))
canvas = FigureCanvasTkAgg(fig, master=window)
canvas.get_tk_widget().pack()

# Buttons for controlling the simulation
start_button = ttk.Button(window, text="Start", command=start_simulation)
start_button.pack(side=tk.LEFT, padx=10, pady=10)

pause_button = ttk.Button(window, text="Pause", command=pause_simulation)
pause_button.pack(side=tk.LEFT, padx=10, pady=10)

reset_button = ttk.Button(window, text="Reset", command=reset_simulation)
reset_button.pack(side=tk.LEFT, padx=10, pady=10)

# Run the tkinter event loop
window.mainloop()
