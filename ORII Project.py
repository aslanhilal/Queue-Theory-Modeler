import tkinter as tk
from tkinter import ttk
from tkinter import StringVar
from tkinter import messagebox
import math

pn_dict = {}


def queuing_model_1(𝜆, µ):  # "M/M/1:GD/∞/∞"
    if 𝜆 >= µ:
        messagebox.showerror("Input Error",
                             "λ (arrival rate) should be less than µ (service rate). Please enter valid values.")
        return None



    else:
        p0 = (µ - 𝜆) / µ
        Ls = 𝜆 / (µ - 𝜆)
        Ws = Ls / 𝜆
        Wq = Ws - (1 / µ)
        Lq = 𝜆 * Wq
        c_bar = 𝜆 / µ
        lam_eff = 𝜆
        lam_lost = 0

        # pn values
        global pn_dict
        for n in range(1, 21):
            pn_dict[f"p_{n}"] = (𝜆 / µ) ** n * p0


        return {"p0": p0, "Ls": Ls, "Ws": Ws, "Wq": Wq, "Lq": Lq, "c_bar": c_bar, "lam_eff": lam_eff,
                "lam_lost": lam_lost}


def queuing_model_2(𝜆, µ, N):  # "M/M/1:GD/N/∞"
    print(f"Queuing Model 2: 𝜆={𝜆}, µ={µ}, N={N}")
    if 𝜆 >= µ:
        messagebox.showerror("Input Error",
                             "λ (arrival rate) should be less than µ (service rate). Please enter valid values.")
        return None
    else:
        rho = 𝜆 / µ
        if rho != 1:
            p0 = (1 - rho) / (1 - rho ** (N + 1))
            pN = (rho ** N) * p0
        else:
            p0 = 1 / (N + 1)
            pN = p0

        global pn_dict
        for n in range(1, 21):
            if n <= N:
                pn = (rho ** n) * p0
            else:
                pn = 0
            pn_dict[f"p_{n}"] = pn

        Ls = sum(n * pn_dict[f"p_{n}"] for n in range(1, N + 1))
        Lq = Ls - (1 - p0)
        lam_loss = 𝜆 * pN
        lam_eff = 𝜆 - lam_loss
        Ws = Ls / lam_eff
        Wq = Lq / lam_eff
        c_bar = Ls - Lq


        return {"p0": p0, "Ls": Ls, "Ws": Ws, "Wq": Wq, "Lq": Lq, "c_bar": c_bar, "lam_eff": lam_eff,
                "lam_lost": lam_loss}


def queuing_model_3(𝜆, µ, c):  # "M/M/c:GD/∞/∞"
    print(f"Queuing Model 3: 𝜆={𝜆}, µ={µ}, c={c}")
    if 𝜆 >= c*µ:
        messagebox.showerror("Input Error",
                             "λ (arrival rate) should be less than µ * c. Please enter valid values.")
        return None
    else:
        summation = sum((𝜆 ** n) / (math.factorial(n) * (µ ** n)) for n in range(c))
        rho = 𝜆 / (c * µ)
        second_term = ((𝜆 ** c) / (math.factorial(c) * (µ ** c))) * (1 / (1 - rho))
        p0 = 1 / (summation + second_term)
        Lq = ((𝜆 ** (c + 1)) / (((c - (𝜆 / µ)) ** 2) * math.factorial(c - 1) * (µ ** (c + 1)))) * p0
        Ls = Lq + (λ / μ)
        Wq = Lq / λ
        Ws = Wq + (1 / μ)
        c_bar = Ls - Lq
        global pn_dict
        for n in range(1, 21):
            if n <= c:
                pn_dict[f"p_{n}"] = ((𝜆) ** n / (math.factorial(n) * (µ ** n))) * p0
            else:
                pn_dict[f"p_{n}"] = (𝜆 ** n / ((c ** (n - c)) * math.factorial(c) * (µ ** n))) * p0
        lam_eff = λ
        lam_lost = 0

        return {"p0": p0, "Ls": Ls, "Ws": Ws, "Wq": Wq, "Lq": Lq, "c_bar": c_bar, "lam_eff": lam_eff, "lam_lost": lam_lost}


def queuing_model_4(𝜆, µ, c, N):  # "M/M/c:GD/N/∞"
    print(f"Queuing Model 4: 𝜆={𝜆}, µ={µ}, c={c}, N={N}")
    summation_term = sum((𝜆 ** n) / (math.factorial(n) * (µ ** n)) for n in range(c))
    second_term_numerator = (𝜆 ** c) * (1 - (𝜆 / (c * µ)) ** (N - c + 1))
    second_term_denominator = (math.factorial(c) * (µ ** c) * (1 - 𝜆 / (c * µ)))
    second_term = second_term_numerator / second_term_denominator
    p0 = 1 / (summation_term + second_term)

    term1 = (𝜆) ** (c + 1)
    term2 = ((c - 𝜆 / µ) ** 2) * math.factorial(c - 1) * µ ** (c + 1)
    term3 = 1 - (𝜆 / (c * µ)) ** (N - c + 1)
    term4 = (N - c + 1) * (1 - (𝜆 / (c * µ))) * (𝜆 / (c * µ)) ** (N - c)

    # Combine the terms to get Lq
    Lq = (term1 / term2) * (term3 - term4) * p0
    global pn_dict
    for n in range(N + 1):
        if n <= c:
            pn = ((𝜆) ** n / (math.factorial(n) * (µ ** n))) * p0
        elif N >= n:
            pn = (𝜆 ** n / ((c ** (n - c)) * math.factorial(c) * (µ ** n))) * p0
        else:
            pn = 0
        pn_dict[f"p_{n}"] = pn
    key = f"p_{n}"
    lam_lost = 𝜆 * pn_dict[key]
    lam_eff = 𝜆 - lam_lost
    print(lam_lost)
    print(lam_eff)
    Ls = Lq + (lam_eff / µ)
    if lam_eff > 0:
        Ws = Ls / lam_eff
        Wq = Lq / lam_eff
    else:
        Ws = 0
        Wq = 0
    c_bar = Ls - Lq
    return {"p0": p0, "Ls": Ls, "Ws": Ws, "Wq": Wq, "Lq": Lq, "c_bar": c_bar, "lam_eff": lam_eff, "lam_lost": lam_lost}


def queuing_model_5(𝜆, µ):
    print(f"Queuing Model 5: 𝜆={𝜆}, µ={µ}")
    if µ == 0:
        print("You can not write µ as 0.")
    else:
        Lq = 0
        Wq = 0
        Ls = 𝜆 / µ
        Ws = 1 / µ
        c_bar = Ls
        lam_eff = 𝜆
        lam_loss = 0
        p0 = math.exp(-𝜆 / µ)
        global pn_dict
        for n in range(1, 21):
            pn_dict[f"p_{n}"] = ((𝜆) ** n / (math.factorial(n) * (µ ** n))) * p0
    return {"p0": p0, "Ls": Ls, "Ws": Ws, "Wq": Wq, "Lq": Lq, "c_bar": c_bar, "lam_eff": lam_eff, "lam_lost": lam_loss}


def hide_or_show(selected, c_label, c_entry, N_label, N_entry, models):
    print(selected)
    if selected == models[0]:
        c_label.grid_forget()
        c_entry.grid_forget()
        N_label.grid_forget()
        N_entry.grid_forget()
    elif selected == models[1]:
        c_label.grid_forget()
        c_entry.grid_forget()
        N_label.grid(row=0, column=2, padx=20, pady=5, sticky="e")
        N_entry.grid(row=0, column=3, padx=10, pady=5, sticky="w")
    elif selected == models[2]:
        c_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")
        c_entry.grid(row=0, column=1, padx=10, pady=5, sticky="w")
        N_label.grid_forget()
        N_entry.grid_forget()
    elif selected == models[3]:
        c_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")
        c_entry.grid(row=0, column=1, padx=10, pady=5, sticky="w")
        N_label.grid(row=0, column=2, padx=20, pady=5, sticky="e")
        N_entry.grid(row=0, column=3, padx=10, pady=5, sticky="w")
    elif selected == models[4]:
        c_label.grid_forget()
        c_entry.grid_forget()
        N_label.grid_forget()
        N_entry.grid_forget()
    else:
        print("Invalid model selection!")


def show_warning():
    messagebox.showwarning("Warning", "Error: Some variables are not filled or non-numerical character has been entered")


def update_values_table(update_values_table, metrics):
    for widget in update_values_table.winfo_children():
        widget.destroy()

    tree = ttk.Treeview(update_values_table, columns=list(metrics.keys()), show="headings")

    for key in metrics.keys():
        tree.heading(key, text=key)
        tree.column(key, width=60, anchor=tk.CENTER)

    tree.insert("", tk.END, values=list(metrics.values()))
    tree.pack(fill=tk.BOTH, expand=True)


def update_p_table(update_p_table, pn_dict):
    for widget in update_p_table.winfo_children():
        widget.destroy()

    tree = ttk.Treeview(update_p_table, columns=("N", "P(N)"), show="headings")
    tree.heading("N", text="N")
    tree.heading("P(N)", text="P(N)")
    tree.column("N", width=100, anchor=tk.CENTER)
    tree.column("P(N)", width=200, anchor=tk.CENTER)

    for key, value in pn_dict.items():
        tree.insert("", tk.END, values=(key, f"{value:.5f}"))

    tree.pack(fill=tk.BOTH, expand=True)


def main():
    global pn_dict
    root = tk.Tk()
    root.title("Queuing Models")
    root.geometry("600x700")
    root.configure(bg="lightblue")
    # Seçenek listesi
    models = ["M/M/1:GD/∞/∞", "M/M/1:GD/N/∞", "M/M/c:GD/∞/∞", "M/M/c:GD/N/∞", "M/M/∞:GD/∞/∞"]

    # Variable to hold the selected value
    selected_option = StringVar()
    selected_option.set(models[0])

    tk.Label(root, text="Operation Research II  Queuing Models Project", bg="lightblue", wraplength=300,
             font=("Arial", 16, "bold")).pack(pady=10)
    dropdown = ttk.Combobox(root, textvariable=selected_option, values=models,
                            state="readonly")  # "readonly" prevents editing
    dropdown.pack(pady=20)

    label_font = ("Arial", 13, "bold")
    frame_lambda_mu = tk.Frame(root, bg="lightblue")
    frame_lambda_mu.pack(padx=20, pady=20)
    tk.Label(frame_lambda_mu, text="𝜆 Value:", bg="lightblue", font=label_font).grid(row=0, column=0, padx=10, pady=5)
    lambda_var = tk.StringVar()
    lambda_entry = tk.Entry(frame_lambda_mu, textvariable=lambda_var, width=20)
    lambda_entry.grid(row=0, column=1, padx=10, pady=5)
    tk.Label(frame_lambda_mu, text="µ Value:", bg="lightblue", font=label_font).grid(row=0, column=2, padx=20,
                                                                                     pady=5)
    mu_var = tk.StringVar()
    mu_entry = tk.Entry(frame_lambda_mu, textvariable=mu_var, width=20)
    mu_entry.grid(row=0, column=3, padx=10, pady=5)

    frame_c_N = tk.Frame(root, bg="lightblue")
    frame_c_N.pack(pady=5)
    c_label = tk.Label(frame_c_N, text="c Value:", bg="lightblue", font=label_font)
    c_var = tk.StringVar()
    c_entry = tk.Entry(frame_c_N, textvariable=c_var, width=20)
    N_label = tk.Label(frame_c_N, text="N Value:", bg="lightblue", font=label_font)
    N_var = tk.StringVar()
    N_entry = tk.Entry(frame_c_N, textvariable=N_var, width=20)

    def on_selection_change(*args):  # This function checks the change on dropdawn menu selection
        hide_or_show(selected_option.get(), c_label, c_entry, N_label, N_entry, models)

    selected_option.trace("w", on_selection_change)

    # calculating function
    def calculate(p_table_frame, values_table_frame):
        try:
            # Receive and transform user inputs
            lambda_value = float(lambda_var.get())
            mu_value = float(mu_var.get())
            c_value = int(c_var.get()) if c_var.get() else None
            N_value = int(N_var.get()) if N_var.get() else None

            # selected model
            selected = selected_option.get()

            # Model selection
            if selected == models[0]:
                metrics = queuing_model_1(lambda_value, mu_value)
            elif selected == models[1]:
                metrics = queuing_model_2(lambda_value, mu_value, N_value)
            elif selected == models[2]:
                metrics = queuing_model_3(lambda_value, mu_value, c_value)
            elif selected == models[3]:
                metrics = queuing_model_4(lambda_value, mu_value, c_value, N_value)
            elif selected == models[4]:
                metrics = queuing_model_5(lambda_value, mu_value)
            else:
                print("!")

            update_p_table(p_table_frame, pn_dict)
            update_values_table(values_table_frame, metrics)

        except ValueError:
            print(ValueError)
            show_warning()

    # Calculation button
    button = ttk.Button(root, text="Calculate", command=lambda: calculate(p_table_frame, values_table_frame))
    button.pack(pady=10)
    p_table_frame = ttk.Frame(root)
    p_table_frame.pack(fill=tk.BOTH, expand=False, padx=30, pady=30)

    values_table_frame = ttk.Frame(root)
    values_table_frame.pack(fill=tk.BOTH, expand=False, padx=40, pady=40)
    root.mainloop()


if __name__ == "__main__":
    main()
