import sys
import traceback

def main():
    try:
        import core_calcolo as engine
    except Exception as e:
        print("Critical error: unable to import core_calcolo.py:", e)
        traceback.print_exc()
        sys.exit(1)

    problems = engine.check_environment()
    if problems:
        print("Attention, the following environmental issues have been detected:")
        for p in problems:
            print(" - " + p)
        print("The app will launch anyway; you can correct the environment and try again from the GUI.\n")
        if getattr(sys, "frozen", False):   
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
        print("Critical error: could not import gui.py:", e)
        print("Verify that tkinter, matplotlib, and pandas are installed "
              "(pip install matplotlib pandas; tkinter is included in standard Python "
              "on Windows/macOS).")
        traceback.print_exc()
        sys.exit(1)

    gui.launch()


if __name__ == "__main__":
    main()
