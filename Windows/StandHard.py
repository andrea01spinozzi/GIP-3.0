import gzip
import io
import os
import queue
import threading
import tkinter as tk
from concurrent.futures import ThreadPoolExecutor, as_completed
from tkinter import filedialog, messagebox, simpledialog, ttk
import pandas as pd
import requests

API = "https://api.gdc.cancer.gov"
TIMEOUT = 60
ANY = "(qualsiasi)"
FACETS = [
    ("Progetto", "cases.project.project_id"),
    ("Sito primario (organo)", "cases.primary_site"),
    ("Tipo di campione", "cases.samples.sample_type"),
    ("Tipo di malattia", "cases.disease_type"),
    ("Diagnosi primaria", "cases.diagnoses.primary_diagnosis"),
    ("Organo di origine (dettaglio)", "cases.diagnoses.tissue_or_organ_of_origin"),
    ("Metastasi (AJCC M)", "cases.diagnoses.ajcc_pathologic_m"),
    ("Stadio patologico (AJCC)", "cases.diagnoses.ajcc_pathologic_stage"),
    ("Genere", "cases.demographic.gender"),
    ("Stato vitale", "cases.demographic.vital_status"),
    ("Tipo di terapia", "cases.diagnoses.treatments.treatment_type"),
    ("Esito della terapia", "cases.diagnoses.treatments.treatment_outcome"),
]
LABEL = dict((f, l) for l, f in FACETS)
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


def list_files(filters, size):
    fields = ",".join(["file_id", "cases.case_id", "cases.samples.submitter_id",
                       "cases.samples.sample_type", "associated_entities.entity_submitter_id"])
    data = gdc_post("files", {"filters": filters, "fields": fields, "size": size, "sort": "file_id"})
    rows = []
    for h in data["hits"]:
        case = (h.get("cases") or [{}])[0]
        barcode = (h.get("associated_entities") or [{}])[0].get("entity_submitter_id", "")
        samples = case.get("samples") or [{}]
        sample = next((s for s in samples if barcode.startswith(s.get("submitter_id", "\0"))), samples[0])
        rows.append({"file_id": h["file_id"], "case_id": case.get("case_id", ""),
                     "original_barcode": barcode, "sample_type": sample.get("sample_type", "")})
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


class App:
    def __init__(self, root):
        self.root = root
        root.title("GDC RNA-Seq Downloader (STAR - Counts)")
        root.geometry("1100x720")
        self.q = queue.Queue()
        self.pool = ThreadPoolExecutor(max_workers=8)
        self.gen = 0
        self.sel = {f: None for _, f in FACETS}
        self.maps = {f: {} for _, f in FACETS}      
        self.vars, self.boxes = {}, {}
        self.total = 0
        self.groups = []                           
        self._build()
        self.refresh()
        self.root.after(100, self.poll)


    def _build(self):
        top = ttk.LabelFrame(self.root, text="1. Filtri collegati al GDC (puoi usarli in qualsiasi ordine)")
        top.pack(fill="x", padx=10, pady=8)
        for i, (label, field) in enumerate(FACETS):
            r, c = divmod(i, 2)
            ttk.Label(top, text=label + ":").grid(row=r, column=c * 2, sticky="e", padx=6, pady=3)
            var = tk.StringVar(value=ANY)
            cb = ttk.Combobox(top, textvariable=var, values=[ANY], state="readonly", width=48)
            cb.grid(row=r, column=c * 2 + 1, sticky="w", padx=6, pady=3)
            cb.bind("<<ComboboxSelected>>", lambda e, f=field: self.on_select(f))
            self.vars[field], self.boxes[field] = var, cb

        bar = ttk.Frame(self.root)
        bar.pack(fill="x", padx=10)
        self.status = ttk.Label(bar, text="", font=("TkDefaultFont", 10, "bold"))
        self.status.pack(side="left")
        ttk.Button(bar, text="Azzera filtri", command=self.reset).pack(side="right")

        add = ttk.LabelFrame(self.root, text="2. Aggiungi la selezione corrente come gruppo/condizione")
        add.pack(fill="x", padx=10, pady=8)
        ttk.Label(add, text="Nome condizione (opzionale):").pack(side="left", padx=6)
        self.cond_var = tk.StringVar()
        ttk.Entry(add, textvariable=self.cond_var, width=35).pack(side="left", padx=6)
        ttk.Label(add, text="N sample:").pack(side="left", padx=6)
        self.n_var = tk.IntVar(value=10)
        ttk.Spinbox(add, from_=1, to=100000, textvariable=self.n_var, width=7).pack(side="left")
        ttk.Button(add, text="➕ Aggiungi gruppo", command=self.add_group).pack(side="left", padx=10)

        mid = ttk.LabelFrame(self.root, text="3. Gruppi da scaricare")
        mid.pack(fill="both", expand=True, padx=10, pady=4)
        cols = ("cond", "filtri", "disp", "n")
        self.tree = ttk.Treeview(mid, columns=cols, show="headings", height=8)
        for c, t, w in [("cond", "Condizione", 220), ("filtri", "Filtri", 560),
                        ("disp", "Disponibili", 90), ("n", "Da scaricare", 100)]:
            self.tree.heading(c, text=t)
            self.tree.column(c, width=w, anchor="w")
        self.tree.pack(side="left", fill="both", expand=True)
        side = ttk.Frame(mid)
        side.pack(side="right", fill="y", padx=6)
        ttk.Button(side, text="Modifica N", command=self.edit_n).pack(fill="x", pady=2)
        ttk.Button(side, text="Rimuovi", command=self.remove_group).pack(fill="x", pady=2)

        bot = ttk.Frame(self.root)
        bot.pack(fill="x", padx=10, pady=8)
        self.btn = ttk.Button(bot, text="⬇ Scarica e unisci", command=self.start_download)
        self.btn.pack(side="left")
        self.prog = ttk.Progressbar(bot, length=500, mode="determinate")
        self.prog.pack(side="left", padx=10)
        self.msg = ttk.Label(bot, text="")
        self.msg.pack(side="left")


    def on_select(self, field):
        label = self.vars[field].get()
        self.sel[field] = None if label == ANY else self.maps[field].get(label)
        self.refresh(skip=field)

    def reset(self):
        for f in self.sel:
            self.sel[f] = None
            self.vars[f].set(ANY)
        self.refresh()

    def refresh(self, skip=None):
        """Per ogni tendina chiede al GDC i valori compatibili con TUTTE le altre tendine."""
        self.gen += 1
        g = self.gen
        self.status.config(text="Aggiornamento dal GDC...")
        sel = dict(self.sel)
        for _, field in FACETS:
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
                elif kind == "done":
                    self.btn.config(state="normal")
                    self.msg.config(text="Completato")
                    messagebox.showinfo("Fatto", payload)
                elif kind == "fail":
                    self.btn.config(state="normal")
                    self.msg.config(text="Errore")
                    messagebox.showerror("Errore", payload)
                elif g != self.gen:
                    continue                            
                elif kind == "facet":
                    need_refresh |= self._apply_facet(field, payload)
                elif kind == "total":
                    self.total = payload
                    self.status.config(text=f"File STAR - Counts corrispondenti ai filtri: {payload}")
                elif kind == "err":
                    self.status.config(text=f"Errore di rete ({field}): {payload}")
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
        self.sel[field], _ = None, self.vars[field].set(ANY)   
        return True


    def add_group(self):
        if self.total <= 0:
            messagebox.showwarning("Nessun dato", "Nessun file corrisponde ai filtri correnti.")
            return
        chosen = {f: v for f, v in self.sel.items() if v}
        desc = "; ".join(f"{LABEL[f]} = {v}" for f, v in chosen.items()) or "nessun filtro"
        name = self.cond_var.get().strip() or " | ".join(chosen.values()) or "Tutti"
        n = max(1, min(int(self.n_var.get()), self.total))
        self.groups.append({"cond": name, "sel": chosen, "disp": self.total, "n": n})
        self.tree.insert("", "end", iid=str(len(self.groups) - 1), values=(name, desc, self.total, n))
        self.cond_var.set("")

    def _redraw(self):
        self.tree.delete(*self.tree.get_children())
        for i, g in enumerate(self.groups):
            desc = "; ".join(f"{LABEL[f]} = {v}" for f, v in g["sel"].items()) or "nessun filtro"
            self.tree.insert("", "end", iid=str(i), values=(g["cond"], desc, g["disp"], g["n"]))

    def remove_group(self):
        for iid in sorted(self.tree.selection(), key=int, reverse=True):
            del self.groups[int(iid)]
        self._redraw()

    def edit_n(self):
        s = self.tree.selection()
        if not s:
            return
        g = self.groups[int(s[0])]
        n = simpledialog.askinteger("Modifica N", f"Sample da scaricare (max {g['disp']}):",
                                    minvalue=1, maxvalue=g["disp"], initialvalue=g["n"])
        if n:
            g["n"] = n
            self._redraw()


    def start_download(self):
        if not self.groups:
            messagebox.showwarning("Nessun gruppo", "Aggiungi almeno un gruppo.")
            return
        outdir = filedialog.askdirectory(title="Cartella di destinazione")
        if not outdir:
            return
        self.btn.config(state="disabled")
        self.msg.config(text="Seleziono i file...")
        threading.Thread(target=self._download, args=(list(self.groups), outdir), daemon=True).start()

    def _download(self, groups, outdir):
        try:
            seen, chosen = set(), []
            for g in groups:
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
            if not chosen:
                raise RuntimeError("Nessun file trovato.")

            series, errors = {}, []
            with ThreadPoolExecutor(max_workers=8) as ex:
                futs = [ex.submit(download_counts, c["file_id"]) for c in chosen]
                for done, f in enumerate(as_completed(futs), 1):
                    try:
                        fid, s = f.result()
                        series[fid] = s
                    except Exception as e: 
                        errors.append(str(e))
                    self.q.put(("prog", done, len(chosen), None))

            ok = [c for c in chosen if c["file_id"] in series]
            meta = pd.DataFrame(ok)
            meta.insert(0, "sample_id", [f"S{i}" for i in range(1, len(ok) + 1)])
            meta = meta[["sample_id", "original_barcode", "condition", "case_id", "sample_type", "file_id"]]
            counts = pd.concat([series[c["file_id"]] for c in ok], axis=1)
            counts.columns = meta["sample_id"].tolist()
            counts.index.name = "gene_id"

            counts.reset_index().to_csv(os.path.join(outdir, "counts_matrix_unito.tsv"), sep="\t", index=False)
            meta.to_csv(os.path.join(outdir, "metadata_unito.tsv"), sep="\t", index=False)
            txt = f"Salvati in {outdir}:\n- counts_matrix_unito.tsv ({counts.shape[0]} geni × {len(ok)} sample)\n- metadata_unito.tsv"
            if errors:
                txt += f"\n\n{len(errors)} file non scaricati."
            self.q.put(("done", 0, None, txt))
        except Exception as e: 
            self.q.put(("fail", 0, None, str(e)))


if __name__ == "__main__":
    root = tk.Tk()
    App(root)
    root.mainloop()
