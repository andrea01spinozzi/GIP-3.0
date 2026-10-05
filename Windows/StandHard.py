from __future__ import annotations

import gzip
import io
import os
import queue
import subprocess
import sys
import threading
import tkinter as tk
from concurrent.futures import ThreadPoolExecutor, as_completed
from tkinter import filedialog, messagebox, scrolledtext, simpledialog, ttk

import pandas as pd
import requests

# --------------------------------------------------------------------------
#  LOGICA GDC (invariata rispetto alla versione precedente)
# --------------------------------------------------------------------------
API = "https://api.gdc.cancer.gov"
TIMEOUT = 60
ANY = "(qualsiasi)"

# (etichetta, campo GDC, suggerimento, avanzato?)
FACETS = [
    ("Progetto", "cases.project.project_id",
     "Lo studio da cui provengono i dati, ad esempio TCGA-BRCA (tumore della mammella).", False),
    ("Sito primario (organo)", "cases.primary_site",
     "L'organo in cui il tumore è nato (es. Breast, Lung, Pancreas).", False),
    ("Tipo di campione", "cases.samples.sample_type",
     "Ad esempio 'Primary Tumor' (tumore) oppure 'Solid Tissue Normal' (tessuto sano).", False),
    ("Tipo di malattia", "cases.disease_type",
     "La categoria generale della malattia.", False),
    ("Genere", "cases.demographic.gender",
     "Genere dichiarato del paziente.", False),
    ("Stato vitale", "cases.demographic.vital_status",
     "Se il paziente risulta vivo o deceduto.", False),
    ("Diagnosi primaria", "cases.diagnoses.primary_diagnosis",
     "La diagnosi istologica dettagliata.", True),
    ("Organo di origine (dettaglio)", "cases.diagnoses.tissue_or_organ_of_origin",
     "Il tessuto o organo di origine, in forma più dettagliata.", True),
    ("Metastasi (AJCC M)", "cases.diagnoses.ajcc_pathologic_m",
     "M0 = nessuna metastasi a distanza, M1 = metastasi presenti.", True),
    ("Stadio patologico (AJCC)", "cases.diagnoses.ajcc_pathologic_stage",
     "Lo stadio del tumore (Stage I, II, III, IV...).", True),
    ("Tipo di terapia", "cases.diagnoses.treatments.treatment_type",
     "Il tipo di trattamento ricevuto (es. chemioterapia, radioterapia).", True),
    ("Esito della terapia", "cases.diagnoses.treatments.treatment_outcome",
     "Il risultato registrato del trattamento.", True),
]
FIELDS = [f for _, f, _, _ in FACETS]
LABEL = {f: l for l, f, _, _ in FACETS}
SAMPLE_TYPE_FIELD = "cases.samples.sample_type"


def _in(field, values):
    return {"op": "in", "content": {"field": field, "value": list(values)}}


def make_filters(selection, exclude=None):
    content = [
        _in("data_category", ["Transcriptome Profiling"]),
        _in("data_type", ["Gene Expression Quantification"]),
        _in("experimental_strategy", ["RNA-Seq"]),
        _in("analysis.workflow_type", ["STAR - Counts"]),
        _in("access", ["open"]),
    ]
    content += [_in(f, [v]) for f, v in selection.items() if v and f != exclude]
    return {"op": "and", "content": content}


def gdc_post(endpoint, payload):
    r = requests.post(f"{API}/{endpoint}", json=payload, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()["data"]


def fetch_facet(field, filters):
    data = gdc_post("files", {"filters": filters, "facets": field, "size": 0})
    return {b["key"]: b["doc_count"] for b in data["aggregations"][field]["buckets"]
            if b["key"] != "_missing"}


def fetch_total(filters):
    return gdc_post("files", {"filters": filters, "size": 0})["pagination"]["total"]


# Informazioni cliniche da riportare nel metadata (colonna di output, in ordine).
# Una colonna viene scritta SOLO se almeno un campione ha un valore; i campioni senza valore
# ricevono NA_LABEL.
NA_LABEL = "not available"
CLINICAL_COLUMNS = [
    "project_id", "primary_site", "disease_type",
    "gender", "age_at_index", "age_at_diagnosis_years", "race", "ethnicity", "vital_status",
    "primary_diagnosis", "ajcc_pathologic_stage", "ajcc_pathologic_m",
    "treatment_type", "treatment_outcome",
]
CLINICAL_FIELDS = [
    "cases.project.project_id", "cases.primary_site", "cases.disease_type",
    "cases.demographic.gender", "cases.demographic.age_at_index", "cases.demographic.race",
    "cases.demographic.ethnicity", "cases.demographic.vital_status",
    "cases.diagnoses.age_at_diagnosis", "cases.diagnoses.primary_diagnosis",
    "cases.diagnoses.ajcc_pathologic_stage", "cases.diagnoses.ajcc_pathologic_m",
    "cases.diagnoses.treatments.treatment_type", "cases.diagnoses.treatments.treatment_outcome",
]
# Valori con cui il GDC indica "dato mancante": per noi equivalgono a "non presente".
_MISSING_MARKERS = {"", "not reported", "unknown", "not allowed to collect", "nan", "none"}


def _clean_value(v):
    """Restituisce il valore come testo pulito, oppure None se mancante."""
    if v is None:
        return None
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    t = " ".join(str(v).replace("\t", " ").split())
    return None if t.lower() in _MISSING_MARKERS else t


def extract_clinical(case):
    """Estrae dal 'case' restituito dal GDC le informazioni cliniche principali."""
    demo = case.get("demographic") or {}
    diags = case.get("diagnoses") or []
    diag = diags[0] if diags else {}
    treatments = [t for d in diags for t in (d.get("treatments") or [])]

    def joined(key):
        vals = []
        for t in treatments:
            v = _clean_value(t.get(key))
            if v and v not in vals:
                vals.append(v)
        return "; ".join(vals) if vals else None

    age_dx = None
    try:
        days = diag.get("age_at_diagnosis")
        if days is not None:
            age_dx = str(int(round(float(days) / 365.25)))   # il GDC lo fornisce in giorni
    except (TypeError, ValueError):
        pass

    project = case.get("project") or {}
    return {
        "project_id": _clean_value(project.get("project_id")),
        "primary_site": _clean_value(case.get("primary_site")),
        "disease_type": _clean_value(case.get("disease_type")),
        "gender": _clean_value(demo.get("gender")),
        "age_at_index": _clean_value(demo.get("age_at_index")),
        "age_at_diagnosis_years": _clean_value(age_dx),
        "race": _clean_value(demo.get("race")),
        "ethnicity": _clean_value(demo.get("ethnicity")),
        "vital_status": _clean_value(demo.get("vital_status")),
        "primary_diagnosis": _clean_value(diag.get("primary_diagnosis")),
        "ajcc_pathologic_stage": _clean_value(diag.get("ajcc_pathologic_stage")),
        "ajcc_pathologic_m": _clean_value(diag.get("ajcc_pathologic_m")),
        "treatment_type": joined("treatment_type"),
        "treatment_outcome": joined("treatment_outcome"),
    }


def add_clinical_columns(meta):
    """Aggiunge a 'meta' le colonne cliniche: scarta quelle senza alcun valore e scrive
    NA_LABEL dove il dato manca solo per alcuni campioni."""
    meta = meta.copy()
    kept = []
    for col in CLINICAL_COLUMNS:
        if col not in meta.columns:
            continue
        s = meta[col].map(_clean_value)
        if s.isna().all():
            meta = meta.drop(columns=[col])
            continue
        meta[col] = s.fillna(NA_LABEL)
        kept.append(col)
    return meta, kept


def list_files(filters, size):
    fields = ",".join(["file_id", "cases.case_id", "cases.samples.submitter_id",
                       "cases.samples.sample_type", "associated_entities.entity_submitter_id"]
                      + CLINICAL_FIELDS)
    data = gdc_post("files", {"filters": filters, "fields": fields, "size": size, "sort": "file_id"})
    rows = []
    for h in data["hits"]:
        case = (h.get("cases") or [{}])[0]
        barcode = (h.get("associated_entities") or [{}])[0].get("entity_submitter_id", "")
        samples = case.get("samples") or [{}]
        sample = next((s for s in samples if barcode.startswith(s.get("submitter_id", "\0"))), samples[0])
        rows.append({"file_id": h["file_id"], "case_id": case.get("case_id", ""),
                     "original_barcode": barcode, "sample_type": sample.get("sample_type", ""),
                     **extract_clinical(case)})
    return rows


def download_counts(file_id):
    last = None
    for _ in range(3):
        try:
            r = requests.get(f"{API}/data/{file_id}", timeout=TIMEOUT * 2)
            r.raise_for_status()
            raw = r.content
            if raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            df = pd.read_csv(io.BytesIO(raw), sep="\t", comment="#", usecols=["gene_id", "unstranded"])
            df = df[~df["gene_id"].str.startswith("N_")]
            return file_id, df.set_index("gene_id")["unstranded"]
        except Exception as e:
            last = e
    raise RuntimeError(f"{file_id}: {last}")


# --------------------------------------------------------------------------
#  TEMA GRAFICO (viola)
# --------------------------------------------------------------------------
PALETTE = {
    "bg": "#f8f5fc",
    "bg_card": "#ffffff",
    "accent": "#4a1d7a",
    "accent_light": "#7b3fb8",
    "accent_pale": "#e4d6f2",
    "accent_pale_2": "#f1eaf9",
    "accent_soft": "#a97fd6",
    "text": "#2a1b3a",
    "muted": "#8b7f99",
    "border": "#e0d4ee",
    "border_strong": "#c9b3e2",
    "ok": "#1e8449",
    "warn": "#b9770e",
    "err": "#c0392b",
}

FONT_BASE = ("Trebuchet MS", 10)
FONT_BOLD = ("Trebuchet MS", 10, "bold")
FONT_SMALL = ("Trebuchet MS", 9)
FONT_HEADER = ("Trebuchet MS", 16, "bold")
FONT_BIG = ("Trebuchet MS", 30, "bold")
FONT_MONO = ("Consolas", 9)

HELP_CONTENT = {
    "filtri": {
        "title": "Aiuto — Passo 1: scegli i pazienti",
        "meaning": (
            "In questo programma scarichi dal GDC (Genomic Data Commons) i dati di espressione "
            "genica RNA-Seq (file 'STAR - Counts') di pazienti oncologici.\n\n"
            "Ogni file corrisponde a UN campione (per esempio un pezzo di tumore di un paziente). "
            "Con i menu a tendina descrivi che tipo di campioni ti interessano: "
            "per esempio solo tumori della mammella, oppure solo tessuto sano.\n\n"
            "Come funzionano i menu:\n"
            "- Il numero tra parentesi quadre accanto a ogni voce è quanti file esistono "
            "per quella voce, tenendo conto delle altre scelte già fatte.\n"
            "- I menu sono collegati: scegliendo un organo, negli altri menu compaiono solo "
            "i valori che esistono davvero per quell'organo.\n"
            "- Puoi usarli in qualsiasi ordine.\n"
            "- La casella grande 'File trovati' ti dice sempre quanti campioni corrispondono "
            "alla tua selezione.\n\n"
            "Un GRUPPO è una selezione che dà un insieme di campioni con caratteristiche comuni "
            "(per esempio 'Tumore'). Di solito servono DUE gruppi da confrontare "
            "(per esempio 'Tumore' e 'Tessuto sano')."
        ),
        "usage": (
            "1. Scegli i filtri nei menu a tendina (di solito bastano Progetto e Tipo di campione). "
            "Con la ✕ accanto a un menu lo azzeri.\n"
            "2. Se ti servono criteri più fini (stadio, terapia, metastasi...) apri "
            "'Mostra filtri avanzati'.\n"
            "3. Controlla che il numero di 'File trovati' sia quello che ti aspetti.\n"
            "4. Scrivi un nome per il gruppo (facoltativo) e il numero di campioni da scaricare "
            "(il pulsante 'Tutti' prende tutti quelli disponibili).\n"
            "5. Premi '➕ Aggiungi come gruppo'.\n"
            "6. Cambia i filtri (per esempio da 'Primary Tumor' a 'Solid Tissue Normal') e "
            "aggiungi il secondo gruppo.\n"
            "7. Quando hai finito, premi 'Avanti' per passare al riepilogo dei gruppi.\n\n"
            "'Azzera tutti i filtri' riporta tutti i menu su '(qualsiasi)'."
        ),
    },
    "gruppi": {
        "title": "Aiuto — Passo 2: gruppi scelti",
        "meaning": (
            "Qui vedi l'elenco dei gruppi che hai aggiunto. Per ogni gruppo trovi:\n"
            "- Condizione: il nome che comparirà nella colonna 'condition' del file metadata.\n"
            "- Filtri: le scelte fatte nei menu.\n"
            "- Disponibili: quanti campioni esistono al GDC per quei filtri.\n"
            "- Da scaricare: quanti ne verranno scaricati davvero.\n\n"
            "Se lo stesso campione compare in due gruppi viene scaricato una volta sola "
            "e assegnato al primo gruppo."
        ),
        "usage": (
            "- Seleziona una riga e premi 'Modifica N' per cambiare quanti campioni scaricare.\n"
            "- 'Rinomina' cambia il nome della condizione (es. 'Tumore', 'Controllo').\n"
            "- 'Rimuovi' elimina il gruppo selezionato.\n"
            "- Quando i gruppi sono a posto premi 'Avanti'."
        ),
    },
    "scarica": {
        "title": "Aiuto — Passo 3: scarica e unisci",
        "meaning": (
            "Il programma scarica un file per ogni campione, li unisce in un'unica tabella e "
            "salva due file nella cartella che scegli:\n\n"
            "- counts_matrix_unito.tsv: una riga per gene e una colonna per campione "
            "(S1, S2, ...), con i conteggi grezzi.\n"
            "- metadata_unito.tsv: per ogni campione indica il codice originale, la condizione "
            "(il nome del gruppo), il paziente, il tipo di campione e il file di origine. "
            "In coda vengono aggiunte le informazioni cliniche disponibili (genere, età, "
            "stadio, terapia, ecc.): una colonna compare solo se almeno un campione ha il dato, "
            "e i campioni che non lo hanno riportano 'not available'.\n\n"
            "I due file sono già nel formato accettato da DEA Explorer."
        ),
        "usage": (
            "1. Premi sul riquadro 'Cartella di destinazione' e scegli dove salvare i file.\n"
            "2. Controlla l'elenco 'Prima di iniziare': deve essere tutto spuntato.\n"
            "3. Premi '⬇ Scarica e unisci' e attendi: la barra mostra l'avanzamento.\n"
            "4. Alla fine il programma ti propone di aprire la cartella.\n\n"
            "Il download può richiedere diversi minuti se i campioni sono molti. "
            "Serve una connessione a Internet attiva. Se qualche file non si scarica, "
            "il programma lo segnala e prosegue con gli altri."
        ),
    },
}


def apply_theme(root: tk.Tk) -> ttk.Style:
    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except tk.TclError:
        pass
    P = PALETTE
    root.configure(bg=P["bg"])

    style.configure(".", background=P["bg"], foreground=P["text"], font=FONT_BASE)
    style.configure("TFrame", background=P["bg"])
    style.configure("Card.TFrame", background=P["bg_card"])
    style.configure("TLabel", background=P["bg"], foreground=P["text"], font=FONT_BASE)
    style.configure("Card.TLabel", background=P["bg_card"], foreground=P["text"], font=FONT_BASE)
    style.configure("CardBold.TLabel", background=P["bg_card"], foreground=P["accent"], font=FONT_BOLD)
    style.configure("Muted.TLabel", background=P["bg"], foreground=P["muted"], font=FONT_SMALL)
    style.configure("CardMuted.TLabel", background=P["bg_card"], foreground=P["muted"], font=FONT_SMALL)
    style.configure("Ok.TLabel", background=P["bg_card"], foreground=P["ok"], font=FONT_BOLD)
    style.configure("Warn.TLabel", background=P["bg_card"], foreground=P["warn"], font=FONT_BOLD)
    style.configure("Err.TLabel", background=P["bg_card"], foreground=P["err"], font=FONT_BOLD)

    style.configure("TLabelframe", background=P["bg_card"], bordercolor=P["border_strong"],
                    relief="solid", borderwidth=1)
    style.configure("TLabelframe.Label", background=P["bg_card"], foreground=P["accent"],
                    font=("Trebuchet MS", 11, "bold"))

    style.configure("TButton", background=P["accent_pale_2"], foreground=P["accent"],
                    padding=(12, 7), font=FONT_BASE, borderwidth=0, relief="flat")
    style.map("TButton",
              background=[("active", P["accent_pale"]), ("disabled", "#eeebf2")],
              foreground=[("disabled", "#a59db0")])

    style.configure("Accent.TButton", background=P["accent"], foreground="white",
                    padding=(18, 10), font=FONT_BOLD, borderwidth=0, relief="flat")
    style.map("Accent.TButton",
              background=[("active", P["accent_light"]), ("disabled", "#c3b3d6")],
              foreground=[("disabled", "white")])

    style.configure("Danger.TButton", background="#f8ebe9", foreground=P["err"],
                    padding=(12, 7), font=FONT_BASE, borderwidth=0)
    style.map("Danger.TButton", background=[("active", "#f1d7d3")])

    style.configure("Help.TButton", background=P["accent_pale"], foreground=P["accent"],
                    padding=6, font=FONT_BOLD, borderwidth=0)
    style.map("Help.TButton", background=[("active", P["accent_light"])],
              foreground=[("active", "white")])

    style.configure("Clear.TButton", background=P["accent_pale_2"], foreground=P["muted"],
                    padding=(6, 4), font=FONT_SMALL, borderwidth=0)
    style.map("Clear.TButton", background=[("active", "#f1d7d3")],
              foreground=[("active", P["err"])])

    style.configure("TNotebook", background=P["bg"], borderwidth=0, tabmargins=(8, 8, 8, 0))
    style.configure("TNotebook.Tab", background=P["accent_pale_2"], foreground=P["accent"],
                    padding=(20, 10), font=FONT_BOLD, borderwidth=0)
    style.map("TNotebook.Tab",
              background=[("selected", P["accent"]), ("active", P["accent_pale"])],
              foreground=[("selected", "white")])

    style.configure("Treeview", background="white", fieldbackground="white", foreground=P["text"],
                    rowheight=28, font=FONT_BASE, borderwidth=0)
    style.configure("Treeview.Heading", background=P["accent"], foreground="white", font=FONT_BOLD,
                    relief="flat", padding=(6, 6))
    style.map("Treeview.Heading", background=[("active", P["accent_light"])])
    style.map("Treeview", background=[("selected", P["accent_pale"])],
              foreground=[("selected", P["text"])])

    style.configure("Horizontal.TProgressbar", background=P["accent_light"],
                    troughcolor=P["accent_pale_2"], borderwidth=0, thickness=16)

    style.configure("TEntry", padding=5, fieldbackground="white")
    style.configure("TSpinbox", padding=5, fieldbackground="white")
    style.configure("TCombobox", padding=5, fieldbackground="white",
                    selectbackground="white", selectforeground=P["text"])
    style.map("TCombobox", fieldbackground=[("readonly", "white")],
              selectbackground=[("readonly", "white")],
              selectforeground=[("readonly", P["text"])])
    style.configure("TScrollbar", background=P["accent_pale_2"], troughcolor=P["bg"],
                    bordercolor=P["border"])
    return style


class Tooltip:
    def __init__(self, widget, text, delay=500, wraplength=300):
        self.widget, self.text = widget, text
        self.delay, self.wraplength = delay, wraplength
        self._after_id = None
        self._tip = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, _e=None):
        self._unschedule()
        self._after_id = self.widget.after(self.delay, self._show)

    def _unschedule(self):
        if self._after_id:
            self.widget.after_cancel(self._after_id)
            self._after_id = None

    def _show(self):
        if self._tip or not self.text:
            return
        x = self.widget.winfo_rootx() + 12
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 8
        self._tip = tk.Toplevel(self.widget)
        self._tip.wm_overrideredirect(True)
        try:
            self._tip.wm_attributes("-topmost", True)
        except tk.TclError:
            pass
        self._tip.wm_geometry(f"+{x}+{y}")
        tk.Label(self._tip, text=self.text, justify="left", background="#fffee0",
                 relief="solid", borderwidth=1, wraplength=self.wraplength,
                 font=FONT_SMALL, padx=6, pady=4).pack()

    def _hide(self, _e=None):
        self._unschedule()
        if self._tip:
            self._tip.destroy()
            self._tip = None


def add_tip(widget, text):
    Tooltip(widget, text)
    return widget


class ScrollableFrame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.canvas = tk.Canvas(self, borderwidth=0, highlightthickness=0, bg=PALETTE["bg"])
        self.vbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas)
        self.canvas.configure(yscrollcommand=self.vbar.set)
        self.vbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self._win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.inner.bind("<Configure>", lambda e: self.canvas.configure(
            scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfig(self._win, width=e.width))
        self.canvas.bind("<Enter>", self._bind_wheel)
        self.canvas.bind("<Leave>", self._unbind_wheel)
        self.inner.bind("<Enter>", self._bind_wheel)

    def _bind_wheel(self, _e):
        self.canvas.bind_all("<MouseWheel>", self.on_wheel)
        self.canvas.bind_all("<Button-4>", self.on_wheel)
        self.canvas.bind_all("<Button-5>", self.on_wheel)

    def _unbind_wheel(self, _e):
        # evita di sganciare la rotellina se il mouse è ancora sopra il contenuto
        x, y = self.winfo_pointerxy()
        w = self.winfo_containing(x, y)
        while w is not None:
            if w is self:
                return
            w = getattr(w, "master", None)
        self.canvas.unbind_all("<MouseWheel>")
        self.canvas.unbind_all("<Button-4>")
        self.canvas.unbind_all("<Button-5>")

    def on_wheel(self, event):
        if event.num == 5 or event.delta < 0:
            self.canvas.yview_scroll(1, "units")
        elif event.num == 4 or event.delta > 0:
            self.canvas.yview_scroll(-1, "units")
        return "break"


class FolderZone(tk.Frame):
    """Riquadro cliccabile per scegliere la cartella di destinazione."""

    def __init__(self, parent, path_var: tk.StringVar, pick_command):
        super().__init__(parent, bg=PALETTE["bg_card"], highlightthickness=2,
                         highlightbackground=PALETTE["border_strong"], bd=0, cursor="hand2")
        self.path_var, self.pick_command = path_var, pick_command
        self._idle, self._hover = PALETTE["bg_card"], PALETTE["accent_pale_2"]

        self.icon_lbl = tk.Label(self, text="📁", font=("Segoe UI Emoji", 28),
                                 bg=self._idle, fg=PALETTE["accent"])
        self.icon_lbl.pack(pady=(14, 2))
        self.title_lbl = tk.Label(self, text="Cartella di destinazione", font=FONT_BOLD,
                                  bg=self._idle, fg=PALETTE["accent"])
        self.title_lbl.pack()
        self.sub_lbl = tk.Label(self, text="Qui verranno salvati i due file TSV",
                                font=FONT_SMALL, bg=self._idle, fg=PALETTE["muted"])
        self.sub_lbl.pack(pady=(2, 8))
        self.path_lbl = tk.Label(self, text="", font=FONT_SMALL, bg=self._idle,
                                 fg=PALETTE["muted"], wraplength=520, justify="center")
        self.path_lbl.pack(padx=10, pady=(0, 8))
        self.btn = ttk.Button(self, text="Scegli cartella...", command=self.pick_command)
        self.btn.pack(pady=(0, 14))

        self._widgets = [self.icon_lbl, self.title_lbl, self.sub_lbl, self.path_lbl]
        for w in [self] + self._widgets:
            w.bind("<Button-1>", lambda _e: self.pick_command())
            w.bind("<Enter>", lambda _e: self._set_bg(self._hover))
            w.bind("<Leave>", lambda _e: self._set_bg(self._idle))
        self.path_var.trace_add("write", lambda *a: self._refresh())
        self._refresh()

    def _set_bg(self, color):
        self.configure(bg=color)
        for w in self._widgets:
            w.configure(bg=color)

    def _refresh(self):
        p = self.path_var.get().strip()
        if not p:
            self.path_lbl.configure(text="Nessuna cartella selezionata — clicca qui",
                                    fg=PALETTE["muted"])
            self.icon_lbl.configure(fg=PALETTE["accent"])
            self.configure(highlightbackground=PALETTE["border_strong"])
        elif os.path.isdir(p):
            self.path_lbl.configure(text=f"✓  {p}", fg=PALETTE["ok"])
            self.icon_lbl.configure(fg=PALETTE["ok"])
            self.configure(highlightbackground=PALETTE["ok"])
        else:
            self.path_lbl.configure(text=f"✗  Cartella non trovata: {p}", fg=PALETTE["err"])
            self.icon_lbl.configure(fg=PALETTE["err"])
            self.configure(highlightbackground=PALETTE["err"])


def open_folder(path):
    try:
        if sys.platform.startswith("win"):
            os.startfile(path)  # type: ignore[attr-defined]
        elif sys.platform == "darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])
    except Exception:
        pass


# --------------------------------------------------------------------------
#  APPLICAZIONE
# --------------------------------------------------------------------------
class App:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title("GDC RNA-Seq Downloader")
        root.geometry("1120x820")
        root.minsize(980, 700)
        apply_theme(root)

        self.q = queue.Queue()
        self.pool = ThreadPoolExecutor(max_workers=8)
        self.gen = 0
        self.busy = False
        self.sel = {f: None for f in FIELDS}
        self.maps = {f: {} for f in FIELDS}
        self.vars, self.boxes = {}, {}
        self.total = 0
        self.groups = []

        self.cond_var = tk.StringVar()
        self.n_var = tk.IntVar(value=10)
        self.total_var = tk.StringVar(value="…")
        self.status_var = tk.StringVar(value="Pronto. Inizia scegliendo i filtri nel passo 1.")
        self.outdir = tk.StringVar()
        self.outdir.trace_add("write", lambda *a: self._update_checks())

        self._build()
        self._update_checks()
        self.refresh()
        self.root.after(100, self.poll)

    # ---------------- costruzione interfaccia ----------------
    def _build(self):
        self._build_header()

        status = tk.Frame(self.root, bg=PALETTE["accent_pale_2"])
        status.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Label(status, textvariable=self.status_var, anchor="w",
                  background=PALETTE["accent_pale_2"], foreground=PALETTE["accent"],
                  font=FONT_SMALL, padding=(10, 6)).pack(fill=tk.X)

        self.nb = ttk.Notebook(self.root)
        self.nb.pack(fill=tk.BOTH, expand=True, padx=10, pady=(8, 0))
        self.tab_filters = ttk.Frame(self.nb)
        self.tab_groups = ttk.Frame(self.nb)
        self.tab_dl = ttk.Frame(self.nb)
        self.nb.add(self.tab_filters, text="1 · Scegli i pazienti")
        self.nb.add(self.tab_groups, text="2 · Gruppi scelti")
        self.nb.add(self.tab_dl, text="3 · Scarica")

        self.sf = ScrollableFrame(self.tab_filters)
        self.sf.pack(fill=tk.BOTH, expand=True)
        self._build_filters_tab(self.sf.inner)
        self._build_groups_tab(self.tab_groups)
        sf3 = ScrollableFrame(self.tab_dl)
        sf3.pack(fill=tk.BOTH, expand=True)
        self._build_download_tab(sf3.inner)

    def _build_header(self):
        header = tk.Frame(self.root, bg=PALETTE["accent"])
        header.pack(side=tk.TOP, fill=tk.X)
        inner = tk.Frame(header, bg=PALETTE["accent"])
        inner.pack(fill=tk.X, padx=24, pady=(16, 14))
        row = tk.Frame(inner, bg=PALETTE["accent"])
        row.pack(fill="x", anchor="w")
        tk.Label(row, text="GDC Downloader", font=("Trebuchet MS", 21, "bold"),
                 fg="white", bg=PALETTE["accent"]).pack(side=tk.LEFT)
        tk.Label(row, text=" RNA-Seq ", font=("Trebuchet MS", 9, "bold"),
                 fg=PALETTE["accent"], bg=PALETTE["accent_soft"], padx=2
                 ).pack(side=tk.LEFT, padx=(10, 0), pady=(4, 0))
        tk.Label(inner, text="Scarica e unisci i dati RNA-Seq (STAR - Counts) del Genomic Data "
                             "Commons in 3 semplici passi",
                 font=("Trebuchet MS", 11), fg=PALETTE["accent_pale"],
                 bg=PALETTE["accent"]).pack(fill="x", anchor="w", pady=(4, 0))
        tk.Frame(self.root, bg=PALETTE["accent_soft"], height=3).pack(side=tk.TOP, fill=tk.X)

    def _help_button(self, parent, key):
        btn = ttk.Button(parent, text="? Aiuto", style="Help.TButton",
                         command=lambda: self._show_help(key))
        add_tip(btn, "Apre la spiegazione di questo passo e di come si usa.")
        return btn

    def _intro(self, parent, text, help_key):
        bar = ttk.Frame(parent)
        bar.pack(fill=tk.X, padx=12, pady=(12, 4))
        ttk.Label(bar, text=text, wraplength=820, justify="left").pack(side=tk.LEFT, anchor="w")
        self._help_button(bar, help_key).pack(side=tk.RIGHT, padx=(10, 0))

    # ---- Tab 1
    def _build_filters_tab(self, f):
        self._intro(f, "Descrivi i pazienti che ti interessano con i menu a tendina. "
                       "Il numero tra parentesi quadre indica quanti file esistono per ogni voce.",
                    "filtri")

        # contatore grande
        counter = tk.Frame(f, bg=PALETTE["accent_pale_2"], highlightthickness=1,
                           highlightbackground=PALETTE["border_strong"])
        counter.pack(fill=tk.X, padx=12, pady=8)
        tk.Label(counter, text="🔎", font=("Segoe UI Emoji", 26), bg=PALETTE["accent_pale_2"],
                 fg=PALETTE["accent"]).pack(side=tk.LEFT, padx=(18, 8), pady=8)
        tk.Label(counter, textvariable=self.total_var, font=FONT_BIG,
                 bg=PALETTE["accent_pale_2"], fg=PALETTE["accent"]).pack(side=tk.LEFT)
        tk.Label(counter, text="File trovati con i filtri attuali\n(1 file = 1 campione)",
                 font=FONT_BASE, justify="left", bg=PALETTE["accent_pale_2"],
                 fg=PALETTE["text"]).pack(side=tk.LEFT, padx=14)
        reset = ttk.Button(counter, text="↺ Azzera tutti i filtri", command=self.reset)
        reset.pack(side=tk.RIGHT, padx=18)
        add_tip(reset, "Riporta tutti i menu su '(qualsiasi)'.")

        basic = ttk.LabelFrame(f, text="Filtri principali")
        basic.pack(fill=tk.X, padx=12, pady=8)
        self._fill_filters(basic, [x for x in FACETS if not x[3]])

        self.adv_open = False
        self.adv_btn = ttk.Button(f, text="▸ Mostra filtri avanzati", command=self._toggle_adv)
        self.adv_btn.pack(anchor="w", padx=12, pady=(2, 4))
        self.adv_frame = ttk.LabelFrame(f, text="Filtri avanzati")
        self._fill_filters(self.adv_frame, [x for x in FACETS if x[3]])

        add = ttk.LabelFrame(f, text="Aggiungi questa selezione come gruppo")
        add.pack(fill=tk.X, padx=12, pady=8)
        self.add_frame = add
        row = ttk.Frame(add, style="Card.TFrame")
        row.pack(fill=tk.X, padx=12, pady=(10, 6))
        ttk.Label(row, text="Nome del gruppo (facoltativo):", style="Card.TLabel"
                  ).grid(row=0, column=0, sticky="w", padx=4, pady=4)
        e = ttk.Entry(row, textvariable=self.cond_var, width=34)
        e.grid(row=0, column=1, sticky="w", padx=4)
        add_tip(e, "Es. 'Tumore' oppure 'Controllo'. Comparirà nella colonna 'condition' "
                   "del metadata. Se lo lasci vuoto viene creato dai filtri scelti.")
        ttk.Label(row, text="Quanti campioni:", style="Card.TLabel"
                  ).grid(row=0, column=2, sticky="e", padx=(20, 4))
        sp = ttk.Spinbox(row, from_=1, to=100000, textvariable=self.n_var, width=8)
        sp.grid(row=0, column=3, sticky="w")
        add_tip(sp, "Numero di campioni da scaricare per questo gruppo.")
        allb = ttk.Button(row, text="Tutti", command=self._all_n)
        allb.grid(row=0, column=4, padx=6)
        add_tip(allb, "Imposta il numero massimo disponibile con i filtri attuali.")

        brow = ttk.Frame(add, style="Card.TFrame")
        brow.pack(fill=tk.X, padx=12, pady=(4, 12))
        ttk.Button(brow, text="➕ Aggiungi come gruppo", style="Accent.TButton",
                   command=self.add_group).pack(side=tk.LEFT)
        self.groups_lbl = ttk.Label(brow, text="Nessun gruppo ancora aggiunto.",
                                    style="CardMuted.TLabel")
        self.groups_lbl.pack(side=tk.LEFT, padx=14)
        ttk.Button(brow, text="Avanti ➜", command=lambda: self.nb.select(self.tab_groups)
                   ).pack(side=tk.RIGHT)

    def _fill_filters(self, parent, facets):
        grid = ttk.Frame(parent, style="Card.TFrame")
        grid.pack(fill=tk.X, padx=12, pady=10)
        grid.columnconfigure(0, weight=1)
        grid.columnconfigure(1, weight=1)
        for i, (label, field, tip, _adv) in enumerate(facets):
            r, c = divmod(i, 2)
            cell = ttk.Frame(grid, style="Card.TFrame")
            cell.grid(row=r, column=c, sticky="ew", padx=8, pady=6)
            ttk.Label(cell, text=label, style="CardBold.TLabel").pack(anchor="w")
            line = ttk.Frame(cell, style="Card.TFrame")
            line.pack(fill=tk.X, pady=(2, 0))
            var = tk.StringVar(value=ANY)
            cb = ttk.Combobox(line, textvariable=var, values=[ANY], state="readonly")
            cb.pack(side=tk.LEFT, fill=tk.X, expand=True)
            cb.bind("<<ComboboxSelected>>", lambda e, fld=field: self.on_select(fld))
            self._no_wheel(cb)
            clr = ttk.Button(line, text="✕", style="Clear.TButton", width=3,
                             command=lambda fld=field: self.clear_one(fld))
            clr.pack(side=tk.LEFT, padx=(4, 0))
            add_tip(cb, tip)
            add_tip(clr, "Azzera questo filtro")
            self.vars[field], self.boxes[field] = var, cb

    def _no_wheel(self, cb):
        """Impedisce che la rotellina cambi per sbaglio il valore del menu: scorre la pagina."""
        def handler(event):
            return self.sf.on_wheel(event)
        for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            cb.bind(seq, handler)

    def _toggle_adv(self):
        self.adv_open = not self.adv_open
        if self.adv_open:
            self.adv_frame.pack(fill=tk.X, padx=12, pady=8, before=self.add_frame)
            self.adv_btn.config(text="▾ Nascondi filtri avanzati")
        else:
            self.adv_frame.pack_forget()
            self.adv_btn.config(text="▸ Mostra filtri avanzati")

    def _all_n(self):
        if self.total > 0:
            self.n_var.set(self.total)

    # ---- Tab 2
    def _build_groups_tab(self, f):
        self._intro(f, "Ecco i gruppi che hai aggiunto. Di solito ne servono due da confrontare "
                       "(per esempio 'Tumore' e 'Controllo').", "gruppi")

        box = ttk.LabelFrame(f, text="Gruppi da scaricare")
        box.pack(fill=tk.BOTH, expand=True, padx=12, pady=8)
        inner = ttk.Frame(box, style="Card.TFrame")
        inner.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        cols = ("cond", "filtri", "disp", "n")
        self.tree = ttk.Treeview(inner, columns=cols, show="headings", height=10,
                                 selectmode="browse")
        for c, t, w in [("cond", "Condizione", 200), ("filtri", "Filtri scelti", 460),
                        ("disp", "Disponibili", 100), ("n", "Da scaricare", 110)]:
            self.tree.heading(c, text=t)
            self.tree.column(c, width=w, anchor="w" if c in ("cond", "filtri") else "center")
        sb = ttk.Scrollbar(inner, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.LEFT, fill=tk.Y)

        side = ttk.Frame(inner, style="Card.TFrame")
        side.pack(side=tk.RIGHT, fill=tk.Y, padx=(10, 0))
        for text, cmd, tip, style in [
            ("✎ Modifica N", self.edit_n, "Cambia quanti campioni scaricare per il gruppo selezionato.", "TButton"),
            ("🏷 Rinomina", self.rename_group, "Cambia il nome della condizione.", "TButton"),
            ("🗑 Rimuovi", self.remove_group, "Elimina il gruppo selezionato.", "Danger.TButton"),
        ]:
            b = ttk.Button(side, text=text, command=cmd, style=style)
            b.pack(fill=tk.X, pady=3)
            add_tip(b, tip)

        self.summary_lbl = ttk.Label(f, text="", style="Muted.TLabel")
        self.summary_lbl.pack(anchor="w", padx=14)

        nav = ttk.Frame(f)
        nav.pack(fill=tk.X, padx=12, pady=10)
        ttk.Button(nav, text="⬅ Indietro", command=lambda: self.nb.select(self.tab_filters)
                   ).pack(side=tk.LEFT)
        ttk.Button(nav, text="Avanti ➜", style="Accent.TButton",
                   command=lambda: self.nb.select(self.tab_dl)).pack(side=tk.RIGHT)

    # ---- Tab 3
    def _build_download_tab(self, f):
        self._intro(f, "Scegli dove salvare i file e premi il pulsante: il programma scarica i "
                       "campioni e li unisce in due file pronti per DEA Explorer.", "scarica")

        chk = ttk.LabelFrame(f, text="Prima di iniziare")
        chk.pack(fill=tk.X, padx=12, pady=8)
        self.chk_groups = ttk.Label(chk, text="", style="CardMuted.TLabel")
        self.chk_groups.grid(row=0, column=0, sticky="w", padx=10, pady=4)
        self.chk_folder = ttk.Label(chk, text="", style="CardMuted.TLabel")
        self.chk_folder.grid(row=1, column=0, sticky="w", padx=10, pady=4)

        zone_box = ttk.Frame(f)
        zone_box.pack(fill=tk.X, padx=12, pady=8)
        self.folder_zone = FolderZone(zone_box, self.outdir, self.pick_folder)
        self.folder_zone.pack(fill=tk.X)

        run = ttk.Frame(f)
        run.pack(fill=tk.X, padx=12, pady=(8, 4))
        self.btn = ttk.Button(run, text="⬇ Scarica e unisci", style="Accent.TButton",
                              command=self.start_download)
        self.btn.pack(side=tk.LEFT)
        self.prog = ttk.Progressbar(run, length=420, mode="determinate")
        self.prog.pack(side=tk.LEFT, padx=14)
        self.msg = ttk.Label(run, text="")
        self.msg.pack(side=tk.LEFT)

        logf = ttk.LabelFrame(f, text="Cosa sta succedendo")
        logf.pack(fill=tk.BOTH, expand=True, padx=12, pady=8)
        self.log_text = scrolledtext.ScrolledText(logf, height=9, wrap="word", font=FONT_MONO,
                                                  bg="white", fg=PALETTE["text"],
                                                  state="disabled", borderwidth=0)
        self.log_text.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        nav = ttk.Frame(f)
        nav.pack(fill=tk.X, padx=12, pady=(0, 12))
        ttk.Button(nav, text="⬅ Indietro", command=lambda: self.nb.select(self.tab_groups)
                   ).pack(side=tk.LEFT)

    # ---------------- aiuto ----------------
    def _show_help(self, key):
        content = HELP_CONTENT[key]
        dlg = tk.Toplevel(self.root)
        dlg.title(content["title"])
        dlg.geometry("660x560")
        dlg.configure(bg=PALETTE["bg_card"])
        dlg.transient(self.root)
        head = tk.Frame(dlg, bg=PALETTE["accent"])
        head.pack(side=tk.TOP, fill=tk.X)
        tk.Label(head, text=content["title"], font=FONT_HEADER, fg="white",
                 bg=PALETTE["accent"]).pack(anchor="w", padx=16, pady=12)
        nb = ttk.Notebook(dlg)
        nb.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        def tab(text):
            t = ttk.Frame(nb)
            w = scrolledtext.ScrolledText(t, wrap="word", font=FONT_BASE, bg="white",
                                          fg=PALETTE["text"], padx=10, pady=10)
            w.insert("1.0", text)
            w.configure(state="disabled")
            w.pack(fill=tk.BOTH, expand=True)
            return t

        nb.add(tab(content["meaning"]), text="Cos'è")
        nb.add(tab(content["usage"]), text="Come si usa")
        btns = ttk.Frame(dlg)
        btns.pack(fill=tk.X, padx=10, pady=(0, 10))
        ttk.Button(btns, text="Chiudi", style="Accent.TButton", command=dlg.destroy
                   ).pack(side=tk.RIGHT)

    # ---------------- filtri collegati al GDC ----------------
    def on_select(self, field):
        label = self.vars[field].get()
        self.sel[field] = None if label == ANY else self.maps[field].get(label)
        self.refresh(skip=field)

    def clear_one(self, field):
        self.sel[field] = None
        self.vars[field].set(ANY)
        self.refresh()

    def reset(self):
        for f in self.sel:
            self.sel[f] = None
            self.vars[f].set(ANY)
        self.refresh()

    def refresh(self, skip=None):
        self.gen += 1
        g = self.gen
        self.total_var.set("…")
        self.status_var.set("Aggiornamento dei dati dal GDC...")
        sel = dict(self.sel)
        for field in FIELDS:
            if field == skip:
                continue
            self.pool.submit(self._facet_job, g, field, make_filters(sel, exclude=field))
        self.pool.submit(self._total_job, g, make_filters(sel))

    def _facet_job(self, g, field, filters):
        try:
            self.q.put(("facet", g, field, fetch_facet(field, filters)))
        except Exception as e:
            self.q.put(("err", g, field, str(e)))

    def _total_job(self, g, filters):
        try:
            self.q.put(("total", g, None, fetch_total(filters)))
        except Exception as e:
            self.q.put(("err", g, "totale", str(e)))

    def poll(self):
        need_refresh = False
        try:
            while True:
                kind, g, field, payload = self.q.get_nowait()
                if kind == "prog":
                    self.prog["value"] = g
                    self.prog["maximum"] = field
                    self.msg.config(text=f"{g}/{field} file")
                    self.status_var.set(f"Download in corso: {g} di {field} file...")
                elif kind == "log":
                    self._log(payload)
                elif kind == "done":
                    self.busy = False
                    self._update_checks()
                    self.msg.config(text="Completato ✓")
                    self.status_var.set("Download completato.")
                    self._log(payload)
                    if messagebox.askyesno("Fatto", payload + "\n\nVuoi aprire la cartella?"):
                        open_folder(self.outdir.get())
                elif kind == "fail":
                    self.busy = False
                    self._update_checks()
                    self.msg.config(text="Errore")
                    self.status_var.set("Si è verificato un errore.")
                    self._log("ERRORE: " + payload)
                    messagebox.showerror("Errore", payload)
                elif g != self.gen:
                    continue
                elif kind == "facet":
                    need_refresh |= self._apply_facet(field, payload)
                elif kind == "total":
                    self.total = payload
                    self.total_var.set(f"{payload:,}".replace(",", "."))
                    self.status_var.set(f"File STAR - Counts corrispondenti ai filtri: {payload}")
                elif kind == "err":
                    self.total_var.set("—")
                    self.status_var.set("Impossibile contattare il GDC. Controlla la connessione "
                                        f"a Internet e riprova. ({field}: {payload[:80]})")
        except queue.Empty:
            pass
        if need_refresh:
            self.refresh()
        self.root.after(100, self.poll)

    def _apply_facet(self, field, buckets):
        items = sorted(buckets.items(), key=lambda kv: -kv[1])
        mapping = {f"{k}  [{n}]": k for k, n in items}
        self.maps[field] = mapping
        self.boxes[field]["values"] = [ANY] + list(mapping)
        cur = self.sel[field]
        if cur is None:
            self.vars[field].set(ANY)
            return False
        for lab, val in mapping.items():
            if val == cur:
                self.vars[field].set(lab)
                return False
        self.sel[field] = None
        self.vars[field].set(ANY)
        return True

    # ---------------- gruppi ----------------
    def add_group(self):
        if self.total <= 0:
            messagebox.showwarning("Nessun dato", "Nessun file corrisponde ai filtri correnti.\n"
                                                  "Prova ad azzerare qualche filtro.")
            return
        chosen = {f: v for f, v in self.sel.items() if v}
        desc = "; ".join(f"{LABEL[f]} = {v}" for f, v in chosen.items()) or "nessun filtro"
        name = self.cond_var.get().strip() or " | ".join(chosen.values()) or "Tutti"
        try:
            n = int(self.n_var.get())
        except (tk.TclError, ValueError):
            messagebox.showwarning("Numero non valido", "Scrivi un numero intero in 'Quanti campioni'.")
            return
        n = max(1, min(n, self.total))
        self.groups.append({"cond": name, "sel": chosen, "disp": self.total, "n": n})
        self._redraw()
        self.cond_var.set("")
        self.status_var.set(f"Gruppo '{name}' aggiunto ({n} campioni). Puoi cambiare i filtri e "
                            "aggiungerne un altro, oppure premere 'Avanti'.")

    def _redraw(self):
        self.tree.delete(*self.tree.get_children())
        for i, g in enumerate(self.groups):
            desc = "; ".join(f"{LABEL[f]} = {v}" for f, v in g["sel"].items()) or "nessun filtro"
            self.tree.insert("", "end", iid=str(i), values=(g["cond"], desc, g["disp"], g["n"]))
        self._update_checks()

    def _selected_index(self):
        s = self.tree.selection()
        if not s:
            messagebox.showinfo("Seleziona un gruppo", "Clicca prima su una riga dell'elenco.")
            return None
        return int(s[0])

    def remove_group(self):
        i = self._selected_index()
        if i is None:
            return
        del self.groups[i]
        self._redraw()

    def edit_n(self):
        i = self._selected_index()
        if i is None:
            return
        g = self.groups[i]
        n = simpledialog.askinteger("Modifica N", f"Campioni da scaricare (massimo {g['disp']}):",
                                    minvalue=1, maxvalue=g["disp"], initialvalue=g["n"],
                                    parent=self.root)
        if n:
            g["n"] = n
            self._redraw()

    def rename_group(self):
        i = self._selected_index()
        if i is None:
            return
        name = simpledialog.askstring("Rinomina", "Nuovo nome della condizione:",
                                      initialvalue=self.groups[i]["cond"], parent=self.root)
        if name and name.strip():
            self.groups[i]["cond"] = name.strip()
            self._redraw()

    # ---------------- controlli di stato ----------------
    def _update_checks(self):
        if not hasattr(self, "chk_groups"):
            return
        ng = len(self.groups)
        tot = sum(g["n"] for g in self.groups)
        if ng == 0:
            self.chk_groups.config(text="○ Nessun gruppo aggiunto (torna al passo 1)",
                                   style="CardMuted.TLabel")
            self.groups_lbl.config(text="Nessun gruppo ancora aggiunto.")
            self.summary_lbl.config(text="")
        elif ng == 1:
            self.chk_groups.config(text=f"⚠ Hai un solo gruppo ({tot} campioni): di solito per un "
                                        "confronto ne servono due", style="Warn.TLabel")
            self.groups_lbl.config(text="1 gruppo aggiunto")
            self.summary_lbl.config(text=f"Totale: {tot} campioni da scaricare.")
        else:
            self.chk_groups.config(text=f"✓ {ng} gruppi pronti ({tot} campioni in totale)",
                                   style="Ok.TLabel")
            self.groups_lbl.config(text=f"{ng} gruppi aggiunti")
            self.summary_lbl.config(text=f"Totale: {tot} campioni da scaricare.")
        d = self.outdir.get().strip()
        if d and os.path.isdir(d):
            self.chk_folder.config(text="✓ Cartella di destinazione scelta", style="Ok.TLabel")
        elif d:
            self.chk_folder.config(text="✗ La cartella scelta non esiste", style="Err.TLabel")
        else:
            self.chk_folder.config(text="○ Cartella di destinazione non ancora scelta",
                                   style="CardMuted.TLabel")
        self.btn.config(state="disabled" if (self.busy or ng == 0) else "normal")

    def _log(self, text):
        self.log_text.configure(state="normal")
        self.log_text.insert("end", text + "\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    # ---------------- download ----------------
    def pick_folder(self):
        d = filedialog.askdirectory(title="Cartella di destinazione")
        if d:
            self.outdir.set(d)

    def start_download(self):
        if not self.groups:
            messagebox.showwarning("Nessun gruppo", "Aggiungi almeno un gruppo (passo 1).")
            return
        outdir = self.outdir.get().strip()
        if not outdir or not os.path.isdir(outdir):
            self.pick_folder()
            outdir = self.outdir.get().strip()
            if not outdir or not os.path.isdir(outdir):
                return
        self.busy = True
        self._update_checks()
        self.prog["value"] = 0
        self.msg.config(text="Seleziono i file...")
        self.status_var.set("Seleziono i file da scaricare...")
        self._log("— Avvio —")
        threading.Thread(target=self._download, args=(list(self.groups), outdir),
                         daemon=True).start()

    def _download(self, groups, outdir):
        log = lambda t: self.q.put(("log", 0, None, t))
        try:
            seen, chosen = set(), []
            for g in groups:
                log(f"Cerco i file del gruppo '{g['cond']}'...")
                rows = list_files(make_filters(g["sel"]), g["disp"] + 500)
                taken = 0
                for r in rows:
                    if taken >= g["n"]:
                        break
                    if r["file_id"] in seen:
                        continue
                    st_sel = g["sel"].get(SAMPLE_TYPE_FIELD)
                    if st_sel and r["sample_type"].lower() != st_sel.lower():
                        continue
                    seen.add(r["file_id"])
                    chosen.append({**r, "condition": g["cond"]})
                    taken += 1
                log(f"  → {taken} file selezionati")
            if not chosen:
                raise RuntimeError("Nessun file trovato.")

            log(f"Scarico {len(chosen)} file (può richiedere alcuni minuti)...")
            series, errors = {}, []
            with ThreadPoolExecutor(max_workers=8) as ex:
                futs = [ex.submit(download_counts, c["file_id"]) for c in chosen]
                for done, fut in enumerate(as_completed(futs), 1):
                    try:
                        fid, s = fut.result()
                        series[fid] = s
                    except Exception as e:
                        errors.append(str(e))
                        log(f"  ! file non scaricato: {str(e)[:120]}")
                    self.q.put(("prog", done, len(chosen), None))

            ok = [c for c in chosen if c["file_id"] in series]
            if not ok:
                raise RuntimeError("Nessun file è stato scaricato correttamente. "
                                   "Controlla la connessione a Internet.")
            log("Unisco i dati e salvo i file...")
            meta = pd.DataFrame(ok)
            meta.insert(0, "sample_id", [f"S{i}" for i in range(1, len(ok) + 1)])
            # Le colonne di base restano nelle stesse posizioni di prima (condition = 3ª colonna);
            # quelle cliniche vengono aggiunte in coda.
            meta, clinical_kept = add_clinical_columns(meta)
            meta = meta[["sample_id", "original_barcode", "condition", "case_id",
                         "sample_type", "file_id"] + clinical_kept]
            counts = pd.concat([series[c["file_id"]] for c in ok], axis=1)
            counts.columns = meta["sample_id"].tolist()
            counts.index.name = "gene_id"

            counts.reset_index().to_csv(os.path.join(outdir, "counts_matrix_unito.tsv"),
                                        sep="\t", index=False)
            meta.to_csv(os.path.join(outdir, "metadata_unito.tsv"), sep="\t", index=False)
            txt = (f"Salvati in {outdir}:\n"
                   f"- counts_matrix_unito.tsv ({counts.shape[0]} geni × {len(ok)} sample)\n"
                   f"- metadata_unito.tsv")
            if clinical_kept:
                txt += f" (con le colonne cliniche: {', '.join(clinical_kept)})"
            else:
                txt += " (nessuna informazione clinica disponibile per questi campioni)"
            if errors:
                txt += f"\n\n{len(errors)} file non scaricati."
            self.q.put(("done", 0, None, txt))
        except Exception as e:
            self.q.put(("fail", 0, None, str(e)))


def main():
    if sys.platform.startswith("win"):
        try:
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except Exception:
            pass
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
