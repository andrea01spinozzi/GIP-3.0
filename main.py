import sys
import traceback


def main():
    try:
        import core_calcolo as engine
    except Exception as e:
        print("Errore critico: impossibile importare core_calcolo.py:", e)
        traceback.print_exc()
        sys.exit(1)

    problems = engine.check_environment()
    if problems:
        print("Attenzione, sono stati rilevati i seguenti problemi di ambiente:")
        for p in problems:
            print(" - " + p)
        print("L'app si avviera' comunque: potrai correggere l'ambiente e riprovare dalla GUI.\n")
        if getattr(sys, "frozen", False):   # app .app senza terminale: mostra una finestra
            try:
                import tkinter as tk
                from tkinter import messagebox
                r = tk.Tk(); r.withdraw()
                messagebox.showwarning("DEA Explorer - environment", "\n\n".join(problems))
                r.destroy()
            except Exception:
                pass

    try:
        import gui
    except Exception as e:
        print("Errore critico: impossibile importare gui.py:", e)
        print("Verifica che tkinter, matplotlib e pandas siano installati "
              "(pip install matplotlib pandas; tkinter e' incluso in Python standard "
              "su Windows/macOS, su Linux potrebbe servire 'sudo apt install python3-tk').")
        traceback.print_exc()
        sys.exit(1)

    gui.launch()


if __name__ == "__main__":
    main()
