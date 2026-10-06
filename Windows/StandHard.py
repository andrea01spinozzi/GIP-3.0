from __future__ import annotations

import csv
import gzip
import io
import os
import queue
import re
import subprocess
import sys
import tarfile
import threading
import time
import tkinter as tk
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from tkinter import filedialog, messagebox, scrolledtext, simpledialog, ttk

import numpy as np
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
#  LOGICA GEO (RNA-seq di ratto, tutte le malattie)
# --------------------------------------------------------------------------
# NB: NCBI genera i conteggi uniformi (pipeline unica) SOLO per umano e topo. Per il ratto
# bisogna usare i file di conteggi caricati dagli autori di ogni studio (formati e ID dei geni
# diversi da studio a studio): qui vengono riconosciuti automaticamente quelli più comuni.
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
ORGANISM = "Rattus norvegicus"
GEO_MAX_FILTER_VALUES = 40          # una caratteristica con più valori diversi non è un filtro utile
DISEASE_HINTS = [
    "", "diabetes", "obesity", "hypertension", "stroke", "ischemia", "myocardial infarction",
    "heart failure", "Alzheimer", "Parkinson", "epilepsy", "depression", "spinal cord injury",
    "traumatic brain injury", "neuropathic pain", "kidney injury", "fibrosis", "liver injury",
    "NAFLD", "arthritis", "osteoporosis", "sepsis", "inflammation", "cancer", "pulmonary hypertension",
]
_GEO_MISSING = {"", "na", "n/a", "nan", "not available", "not reported", "not applicable",
                "unknown", "--", "-"}
_COUNT_COLS = ["count", "counts", "raw_count", "raw_counts", "readcount", "read_count",
               "numreads", "expected_count", "unstranded", "htseq_count"]
_SUPP_OK = re.compile(r"\.(txt|tsv|csv|tab|xlsx|xls|counts?)(\.gz)?$", re.I)
_SUPP_BAD = re.compile(r"(?<![a-z])(fpkm|rpkm|tpm|cpm|normali[sz]ed|norm|rlog|vst|log2|degs?|diff\w*|"
                       r"deseq2?_res\w*|filelist|readme|md5)(?![a-z])", re.I)

_net_lock = threading.Lock()
_last_eutils = [0.0]


def _throttle_eutils():
    """E-utilities: massimo ~3 richieste al secondo senza chiave API."""
    with _net_lock:
        wait = 0.35 - (time.time() - _last_eutils[0])
        if wait > 0:
            time.sleep(wait)
        _last_eutils[0] = time.time()


def geo_get(url, params=None, stream=False, timeout=(20, 120)):
    last = "nessuna risposta"
    for attempt in range(4):
        try:
            if url.startswith(EUTILS):
                _throttle_eutils()
            r = requests.get(url, params=params, timeout=timeout, stream=stream)
            if r.status_code == 404:
                raise FileNotFoundError(url)
            if r.status_code in (429, 500, 502, 503, 504):
                last = f"HTTP {r.status_code}"
                time.sleep(1.5 * (attempt + 1))
                continue
            r.raise_for_status()
            return r
        except FileNotFoundError:
            raise
        except requests.RequestException as e:
            last = str(e)
            time.sleep(1 + attempt)
    raise RuntimeError(f"NCBI non raggiungibile ({last})")


def search_geo_series(text, only_with_files=True, retmax=150):
    """Cerca studi GEO RNA-seq di ratto. Ritorna ([{gse,title,n_samples,year,suppfile}], totale)."""
    term = (f'"{ORGANISM}"[Organism] AND "expression profiling by high throughput sequencing"'
            '[DataSet Type] AND gse[Entry Type]')
    text = (text or "").strip()
    if text:
        term += f" AND ({text})"
    r = geo_get(f"{EUTILS}/esearch.fcgi",
                params={"db": "gds", "term": term, "retmax": retmax, "retmode": "json"})
    res = r.json()["esearchresult"]
    ids = res.get("idlist", [])
    total = int(res.get("count", 0))
    if not ids:
        return [], total
    r = geo_get(f"{EUTILS}/esummary.fcgi",
                params={"db": "gds", "id": ",".join(ids), "retmode": "json"})
    data = r.json()["result"]
    rows = []
    for uid in data.get("uids", []):
        d = data[uid]
        acc = str(d.get("accession", ""))
        if not acc.startswith("GSE"):
            continue
        supp = str(d.get("suppfile", "") or "").strip()
        if only_with_files and not supp:
            continue
        rows.append({"gse": acc, "title": str(d.get("title", "")),
                     "n_samples": int(d.get("n_samples") or 0),
                     "year": str(d.get("pdat", ""))[:4], "suppfile": supp})
    rows.sort(key=lambda x: x["year"], reverse=True)
    return rows, total


def geo_series_url(gse):
    num = gse[3:]
    return f"{GEO_FTP}/GSE{num[:-3]}nnn/{gse}"


def geo_list_dir(url):
    """Nomi dei file in una cartella dell'FTP di GEO (lista vuota se la cartella non esiste)."""
    try:
        html = geo_get(url.rstrip("/") + "/").text
    except FileNotFoundError:
        return []
    names = re.findall(r'href="([^"?#/][^"]*)"', html)
    return [urllib.parse.unquote(n) for n in names
            if not n.endswith("/") and not n.startswith(("http", "mailto"))]


def _snake(key):
    s = re.sub(r"[^0-9A-Za-z]+", "_", str(key)).strip("_").lower()
    return s or "campo"


def _clean_geo_value(v):
    t = " ".join(str(v).replace("\t", " ").split())
    return NA_LABEL if t.lower() in _GEO_MISSING else t


def parse_series_matrix(text):
    """Legge le righe !Sample_* di un file series_matrix. Ritorna (DataFrame, {colonna: nome originale}).
    Colonne: gsm, title, source_name + una colonna per ogni caratteristica ('chiave: valore')."""
    rows = {}
    for line in text.splitlines():
        if line.startswith("!series_matrix_table_begin"):
            break
        if not line.startswith("!Sample_"):
            continue
        parts = next(csv.reader([line], delimiter="\t", quotechar='"'))
        rows.setdefault(parts[0][len("!Sample_"):], []).append(parts[1:])
    gsm = [g.strip() for g in rows.get("geo_accession", [[]])[0]]
    n = len(gsm)
    if n == 0:
        raise ValueError("il file series matrix non contiene campioni")

    def first(key):
        r = rows.get(key)
        return ((r[0] if r else []) + [""] * n)[:n]

    titles, sources = first("title"), first("source_name_ch1")
    recs = [{"gsm": gsm[i], "title": titles[i].strip(), "source_name": sources[i].strip()}
            for i in range(n)]
    display = {"source_name": "source name"}
    for row in rows.get("characteristics_ch1", []):
        row = (row + [""] * n)[:n]
        for i, cell in enumerate(row):
            cell = cell.strip()
            if not cell:
                continue
            k, v = cell.split(":", 1) if ":" in cell else ("characteristics", cell)
            col = _snake(k)
            if col in ("gsm", "title", "source_name"):
                col += "_char"
            display.setdefault(col, k.strip())
            old = recs[i].get(col)
            recs[i][col] = v.strip() if old is None else f"{old}; {v.strip()}"
    df = pd.DataFrame(recs)
    for c in df.columns:
        if c not in ("gsm", "title"):
            df[c] = df[c].map(_clean_geo_value)
    return df, display


def fetch_series_samples(gse):
    """Scarica i metadati dei campioni di uno studio. Ritorna (DataFrame, {colonna: nome originale})."""
    base = geo_series_url(gse)
    names = [n for n in geo_list_dir(f"{base}/matrix") if n.endswith("series_matrix.txt.gz")]
    if not names:
        raise RuntimeError(f"{gse}: file dei metadati (series matrix) non trovato su GEO.")
    frames, display = [], {}
    for n in names:
        raw = geo_get(f"{base}/matrix/{n}").content
        df, disp = parse_series_matrix(gzip.decompress(raw).decode("utf-8", errors="replace"))
        frames.append(df)
        for k, v in disp.items():
            display.setdefault(k, v)
    df = pd.concat(frames, ignore_index=True).drop_duplicates("gsm").reset_index(drop=True)
    for c in df.columns:
        if c not in ("gsm", "title"):
            df[c] = df[c].fillna(NA_LABEL)
    return df, display


def rank_counts_files(names):
    """File supplementari ordinati per probabilità di contenere i conteggi grezzi: [(nome, punteggio)]."""
    out = []
    for n in names:
        low = n.lower()
        if low.endswith("_raw.tar"):
            out.append((n, 4))          # un file per campione dentro il tar
            continue
        if not _SUPP_OK.search(low):
            continue
        score = 1
        if "count" in low:
            score += 10
        if "raw" in low:
            score += 3
        if re.search(r"matrix|merged|all|gene", low):
            score += 2
        if _SUPP_BAD.search(low):
            score -= 20                  # FPKM/TPM/normalizzati/risultati DE: non sono conteggi grezzi
        out.append((n, score))
    out.sort(key=lambda x: -x[1])
    return out


def download_bytes(url):
    last = None
    for _ in range(3):
        try:
            r = geo_get(url, stream=True, timeout=(20, 300))
            buf = io.BytesIO()
            for chunk in r.iter_content(1 << 20):
                buf.write(chunk)
            return buf.getvalue()
        except FileNotFoundError:
            raise
        except Exception as e:
            last = e
            time.sleep(2)
    raise RuntimeError(f"download non riuscito: {last}")


def _read_table(name, raw):
    """Legge un file di testo/Excel (anche .gz) in un DataFrame con intestazione."""
    low = name.lower()
    if low.endswith(".gz") or raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
        low = low[:-3] if low.endswith(".gz") else low
    if low.endswith((".xlsx", ".xls")):
        return pd.read_excel(io.BytesIO(raw), sheet_name=0)
    lines = raw.decode("utf-8-sig", errors="replace").splitlines()
    k = 0
    while k < len(lines) and lines[k].startswith("#"):      # commenti (es. featureCounts)
        k += 1
    lines = lines[k:]
    if not lines:
        raise ValueError("file vuoto")
    head = lines[0]
    if low.endswith(".csv"):
        sep = ","
    elif "\t" in head:
        sep = "\t"
    elif ";" in head and "," not in head:
        sep = ";"
    elif "," in head:
        sep = ","
    else:
        sep = r"\s+"
    return pd.read_csv(io.StringIO("\n".join(lines)), sep=sep,
                       engine="python" if sep == r"\s+" else "c")


def _norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def match_columns(columns, samples):
    """Abbina le colonne di una matrice ai campioni GEO. samples = [(gsm, titolo)].
    Ritorna {gsm: nome colonna}. Regole (in ordine): codice GSM nel nome; titolo identico;
    titolo contenuto nel nome (solo se l'abbinamento è univoco)."""
    cols = [str(c) for c in columns]
    gsms = {g for g, _ in samples}
    mapping, used = {}, set()
    for c in cols:
        m = re.search(r"GSM\d+", c)
        if m and m.group(0) in gsms and m.group(0) not in mapping:
            mapping[m.group(0)] = c
            used.add(c)
    ncols = {c: _norm(c) for c in cols if c not in used}
    for g, t in samples:
        if g in mapping or not _norm(t):
            continue
        hit = [c for c, n in ncols.items() if n == _norm(t) and c not in used]
        if len(hit) == 1:
            mapping[g] = hit[0]
            used.add(hit[0])
    rest = [(g, _norm(t)) for g, t in samples if g not in mapping and len(_norm(t)) >= 3]
    for g, nt in rest:
        hit = [c for c, n in ncols.items() if c not in used and (nt in n or (len(n) >= 3 and n in nt))]
        if len(hit) == 1:
            n = ncols[hit[0]]
            if not [g2 for g2, nt2 in rest if g2 != g and (nt2 in n or n in nt2)]:
                mapping[g] = hit[0]
                used.add(hit[0])
    return mapping


def counts_from_table(df, samples):
    """Matrice geni × campioni (colonne = codici GSM) da una tabella con i geni nella prima colonna."""
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    mapping = match_columns(list(df.columns[1:]), samples)
    if not mapping:
        raise ValueError("nessuna colonna del file corrisponde ai campioni (né per codice GSM né per "
                         f"titolo). Prime colonne del file: {', '.join(list(df.columns[1:9]))}")
    gene = df.iloc[:, 0].astype(str).str.strip()
    sub = df[[mapping[g] for g in mapping]].apply(pd.to_numeric, errors="coerce")
    sub.columns = list(mapping)
    sub.index = gene
    sub = sub[~sub.index.isin(["", "nan"])]
    if sub.index.has_duplicates:
        sub = sub.groupby(level=0, sort=False).sum(min_count=1)
    return sub


def parse_single_counts(name, raw):
    """File di conteggi di un solo campione (HTSeq, STAR ReadsPerGene, tabella gene/conteggio)."""
    if name.lower().endswith(".gz") or raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    lines = [l for l in raw.decode("utf-8", errors="replace").splitlines()
             if l.strip() and not l.startswith("#")]
    if not lines:
        raise ValueError("file vuoto")
    sep = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else r"\s+")
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep=sep, header=None, dtype=str,
                     engine="python" if sep == r"\s+" else "c")
    first = df.iloc[0]
    if pd.to_numeric(first.iloc[1:], errors="coerce").isna().all():      # riga di intestazione
        names = [str(x).strip().lower() for x in first]
        df = df.iloc[1:]
        pick = next((names.index(c) for c in _COUNT_COLS if c in names[1:]), None)
        if pick is None:
            c2 = [i for i, n in enumerate(names) if i > 0 and "count" in n]
            if c2:
                pick = c2[0]
            elif df.shape[1] == 2:
                pick = 1
            else:
                raise ValueError("colonna dei conteggi non riconosciuta")
    else:
        pick = 1                                                          # senza intestazione
    idx = df.iloc[:, 0].astype(str).str.strip()
    s = pd.Series(pd.to_numeric(df.iloc[:, pick], errors="coerce").values, index=idx.values)
    s = s[~s.index.str.startswith("__")
          & ~s.index.isin(["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"])]
    s = s.dropna()
    if s.index.has_duplicates:
        s = s.groupby(level=0, sort=False).sum()
    return s


def read_raw_tar(raw, samples):
    """Estrae dal tar *_RAW.tar i file dei campioni richiesti. Ritorna (matrice, [errori])."""
    wanted = {g for g, _ in samples}
    series, errors = {}, []
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as tf:
        for m in tf:
            if not m.isfile():
                continue
            h = re.search(r"GSM\d+", m.name)
            if not h or h.group(0) not in wanted or h.group(0) in series:
                continue
            try:
                series[h.group(0)] = parse_single_counts(m.name, tf.extractfile(m).read())
            except Exception as e:
                errors.append(f"{m.name}: {e}")
    if not series:
        raise ValueError("nel file .tar non ci sono file leggibili per i campioni scelti"
                         + (f" ({errors[0]})" if errors else ""))
    return pd.concat(series, axis=1), errors


def fetch_series_counts(gse, fname, samples, log=None):
    """Scarica e legge i conteggi di uno studio. Ritorna una matrice geni × campioni (colonne = GSM)."""
    log = log or (lambda t: None)
    url = f"{geo_series_url(gse)}/suppl/{urllib.parse.quote(fname)}"
    log(f"  scarico {fname} ...")
    raw = download_bytes(url)
    log(f"  {len(raw) / 1e6:.1f} MB scaricati, leggo i conteggi...")
    if fname.lower().endswith(".tar"):
        counts, errors = read_raw_tar(raw, samples)
        for e in errors[:5]:
            log(f"  ! {e[:140]}")
        return counts
    return counts_from_table(_read_table(fname, raw), samples)


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


HELP_CONTENT["geo"] = {
    "title": "Aiuto — Dati di ratto da GEO",
    "meaning": (
        "GEO (Gene Expression Omnibus) è l'archivio pubblico dove i ricercatori depositano i loro "
        "esperimenti di espressione genica. Con questa opzione cerchi studi RNA-seq sul ratto "
        "(Rattus norvegicus) per QUALSIASI malattia, non solo i tumori.\n\n"
        "Differenze rispetto al GDC, importanti da conoscere:\n"
        "- Su GEO ogni studio è un progetto separato, con le sue etichette. Quindi prima si sceglie "
        "UNO studio, poi si dividono i suoi campioni in gruppi (es. 'malato' e 'controllo') usando le "
        "caratteristiche che gli autori hanno inserito (genotipo, trattamento, tessuto...).\n"
        "- NCBI ricalcola i conteggi in modo uniforme solo per umano e topo, NON per il ratto. Qui "
        "si usano quindi i file di conteggi caricati dagli autori: il programma riconosce i formati "
        "più comuni, ma può capitare uno studio che non si riesce a leggere.\n"
        "- Ogni studio può usare identificatori di gene diversi (simboli, Ensembl, Entrez). "
        "Il programma lo segnala nel registro, ma non li converte.\n"
        "- Le etichette dei campioni sono scritte dagli autori: controlla sempre che il filtro scelto "
        "corrisponda davvero a ciò che intendi."
    ),
    "usage": (
        "1. Scrivi una malattia o una parola chiave in inglese (o scegli un suggerimento) e premi 'Cerca'.\n"
        "2. Seleziona uno studio dall'elenco (doppio clic o 'Usa lo studio selezionato'): guarda titolo e "
        "numero di campioni.\n"
        "3. Controlla 'File dei conteggi': deve essere un file di conteggi grezzi (raw counts), "
        "non FPKM/TPM.\n"
        "4. Usa i filtri per scegliere i campioni di un gruppo (es. 'treatment = vehicle'), poi "
        "'Aggiungi come gruppo'. Cambia i filtri e aggiungi il secondo gruppo.\n"
        "5. Prosegui con 'Avanti' come per il GDC.\n\n"
        "Si possono unire gruppi di studi diversi solo se usano gli stessi identificatori di gene, "
        "ma è sconsigliato: gli studi diversi hanno differenze tecniche (effetto batch) che si "
        "confondono con la malattia. Per un confronto affidabile usa due gruppi dello stesso studio."
    ),
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

    style.configure("TRadiobutton", background=P["bg"], foreground=P["accent"], font=FONT_BOLD)
    style.map("TRadiobutton", background=[("active", P["bg"])])
    style.configure("Card.TCheckbutton", background=P["bg_card"], foreground=P["text"], font=FONT_BASE)
    style.map("Card.TCheckbutton", background=[("active", P["bg_card"])])
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
        root.title("RNA-Seq Downloader (GDC · GEO)")
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
        self.mode = tk.StringVar(value="gdc")           # "gdc" | "geo"
        self._mode_prev = "gdc"
        self.counter_txt = tk.StringVar(value="File trovati con i filtri attuali\n(1 file = 1 campione)")
        self.geo_query = tk.StringVar()
        self.geo_only_files = tk.BooleanVar(value=True)
        self.geo_msg = tk.StringVar(value="Cerca uno studio per iniziare.")
        self.geo_title_var = tk.StringVar(value="Nessuno studio caricato")
        self.geo_file_var = tk.StringVar()
        self.geo_series = None
        self.geo_sel, self.geo_vars, self.geo_boxes, self.geo_maps = {}, {}, {}, {}

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
        tk.Label(row, text="RNA-Seq Downloader", font=("Trebuchet MS", 21, "bold"),
                 fg="white", bg=PALETTE["accent"]).pack(side=tk.LEFT)
        tk.Label(row, text=" GDC · GEO ", font=("Trebuchet MS", 9, "bold"),
                 fg=PALETTE["accent"], bg=PALETTE["accent_soft"], padx=2
                 ).pack(side=tk.LEFT, padx=(10, 0), pady=(4, 0))
        tk.Label(inner, text="Scarica e unisci dati RNA-Seq: tumori umani dal GDC, oppure ratto "
                             "(tutte le malattie) da GEO, in 3 semplici passi",
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

        src = ttk.Frame(f)
        src.pack(fill=tk.X, padx=12, pady=(4, 0))
        ttk.Label(src, text="Origine dei dati:", style="Muted.TLabel").pack(side=tk.LEFT, padx=(0, 10))
        for text, val in [("Tumori umani · GDC", "gdc"), ("Ratto, tutte le malattie · GEO", "geo")]:
            ttk.Radiobutton(src, text=text, value=val, variable=self.mode,
                            command=self._switch_source).pack(side=tk.LEFT, padx=(0, 18))

        # contatore grande
        counter = tk.Frame(f, bg=PALETTE["accent_pale_2"], highlightthickness=1,
                           highlightbackground=PALETTE["border_strong"])
        counter.pack(fill=tk.X, padx=12, pady=8)
        tk.Label(counter, text="🔎", font=("Segoe UI Emoji", 26), bg=PALETTE["accent_pale_2"],
                 fg=PALETTE["accent"]).pack(side=tk.LEFT, padx=(18, 8), pady=8)
        tk.Label(counter, textvariable=self.total_var, font=FONT_BIG,
                 bg=PALETTE["accent_pale_2"], fg=PALETTE["accent"]).pack(side=tk.LEFT)
        tk.Label(counter, textvariable=self.counter_txt,
                 font=FONT_BASE, justify="left", bg=PALETTE["accent_pale_2"],
                 fg=PALETTE["text"]).pack(side=tk.LEFT, padx=14)
        reset = ttk.Button(counter, text="↺ Azzera tutti i filtri", command=self.reset)
        reset.pack(side=tk.RIGHT, padx=18)
        add_tip(reset, "Riporta tutti i menu su '(qualsiasi)'.")

        self.gdc_frame = ttk.Frame(f)
        basic = ttk.LabelFrame(self.gdc_frame, text="Filtri principali")
        basic.pack(fill=tk.X, padx=12, pady=8)
        self._fill_filters(basic, [x for x in FACETS if not x[3]])

        self.adv_open = False
        self.adv_btn = ttk.Button(self.gdc_frame, text="▸ Mostra filtri avanzati", command=self._toggle_adv)
        self.adv_btn.pack(anchor="w", padx=12, pady=(2, 4))
        self.adv_frame = ttk.LabelFrame(self.gdc_frame, text="Filtri avanzati")
        self._fill_filters(self.adv_frame, [x for x in FACETS if x[3]])
        self._build_geo_panel(f)

        add = ttk.LabelFrame(f, text="Aggiungi questa selezione come gruppo")
        add.pack(fill=tk.X, padx=12, pady=8)
        self.add_frame = add
        self.gdc_frame.pack(fill=tk.X, before=add)
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
            self.adv_frame.pack(fill=tk.X, padx=12, pady=8, after=self.adv_btn)
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

    # ---------------- origine dati: GDC / GEO ----------------
    def _switch_source(self):
        new = self.mode.get()
        if new == self._mode_prev:
            return
        if self.groups and not messagebox.askyesno(
                "Cambia origine dei dati",
                "Cambiando origine i gruppi già aggiunti verranno eliminati.\nVuoi continuare?"):
            self.mode.set(self._mode_prev)
            return
        self.groups.clear()
        self._redraw()
        self._mode_prev = new
        self.gen += 1                                   # scarta le risposte GDC ancora in arrivo
        if new == "geo":
            self.gdc_frame.pack_forget()
            self.geo_frame.pack(fill=tk.X, padx=12, pady=8, before=self.add_frame)
            self.counter_txt.set("Campioni dello studio con i filtri attuali\n(1 campione = 1 animale/prelievo)")
            self.status_var.set("Origine: GEO (ratto). Cerca uno studio e caricalo.")
            self._geo_refresh_filters()
        else:
            self.geo_frame.pack_forget()
            self.gdc_frame.pack(fill=tk.X, before=self.add_frame)
            self.counter_txt.set("File trovati con i filtri attuali\n(1 file = 1 campione)")
            self.refresh()

    # ---------------- GEO: pannello ----------------
    def _build_geo_panel(self, parent):
        self.geo_frame = ttk.Frame(parent)
        a = ttk.LabelFrame(self.geo_frame, text="A · Cerca uno studio sul ratto")
        a.pack(fill=tk.X, pady=(0, 8))
        row = ttk.Frame(a, style="Card.TFrame")
        row.pack(fill=tk.X, padx=12, pady=(10, 4))
        ttk.Label(row, text="Malattia o parola chiave:", style="Card.TLabel").pack(side=tk.LEFT)
        cb = ttk.Combobox(row, textvariable=self.geo_query, values=DISEASE_HINTS, width=34)
        cb.pack(side=tk.LEFT, padx=8)
        cb.bind("<Return>", lambda e: self.geo_search())
        self._no_wheel(cb)
        add_tip(cb, "Scrivi una malattia, un tessuto o un trattamento (in inglese, es. 'diabetes', "
                    "'liver fibrosis') oppure scegli un suggerimento. Vuoto = tutti gli studi sul ratto.")
        go = ttk.Button(row, text="🔍 Cerca", style="Accent.TButton", command=self.geo_search)
        go.pack(side=tk.LEFT)
        self._help_button(row, "geo").pack(side=tk.RIGHT)
        ttk.Checkbutton(a, text="Mostra solo studi con file supplementari (dove stanno i conteggi)",
                        variable=self.geo_only_files, style="Card.TCheckbutton"
                        ).pack(anchor="w", padx=12, pady=(0, 4))

        tf = ttk.Frame(a, style="Card.TFrame")
        tf.pack(fill=tk.X, padx=12, pady=4)
        self.geo_tree = ttk.Treeview(tf, columns=("gse", "year", "n", "title"), show="headings",
                                     height=7, selectmode="browse")
        for c, t, w, anc in [("gse", "Studio", 100, "w"), ("year", "Anno", 60, "center"),
                             ("n", "Campioni", 80, "center"), ("title", "Titolo", 620, "w")]:
            self.geo_tree.heading(c, text=t)
            self.geo_tree.column(c, width=w, anchor=anc)
        sb = ttk.Scrollbar(tf, orient="vertical", command=self.geo_tree.yview)
        self.geo_tree.configure(yscrollcommand=sb.set)
        self.geo_tree.pack(side=tk.LEFT, fill=tk.X, expand=True)
        sb.pack(side=tk.LEFT, fill=tk.Y)
        self.geo_tree.bind("<Double-1>", lambda e: self.geo_load())

        brow = ttk.Frame(a, style="Card.TFrame")
        brow.pack(fill=tk.X, padx=12, pady=(4, 12))
        ttk.Button(brow, text="Usa lo studio selezionato ➜", command=self.geo_load).pack(side=tk.LEFT)
        ttk.Label(brow, textvariable=self.geo_msg, style="CardMuted.TLabel", wraplength=640,
                  justify="left").pack(side=tk.LEFT, padx=12)

        b = ttk.LabelFrame(self.geo_frame, text="B · Scegli i campioni dello studio")
        b.pack(fill=tk.X)
        ttk.Label(b, textvariable=self.geo_title_var, style="CardBold.TLabel", wraplength=900,
                  justify="left").pack(anchor="w", padx=12, pady=(10, 4))
        frow = ttk.Frame(b, style="Card.TFrame")
        frow.pack(fill=tk.X, padx=12, pady=4)
        ttk.Label(frow, text="File dei conteggi:", style="Card.TLabel").pack(side=tk.LEFT)
        self.geo_file_cb = ttk.Combobox(frow, textvariable=self.geo_file_var, state="readonly", width=60)
        self.geo_file_cb.pack(side=tk.LEFT, padx=8)
        self._no_wheel(self.geo_file_cb)
        add_tip(self.geo_file_cb, "File supplementare dello studio che contiene i conteggi grezzi. "
                                  "Il programma sceglie quello più probabile: cambialo solo se serve.")
        self.geo_filters_frame = ttk.Frame(b, style="Card.TFrame")
        self.geo_filters_frame.pack(fill=tk.X, padx=4, pady=(4, 10))
        ttk.Label(self.geo_filters_frame, text="Carica prima uno studio (passo A).",
                  style="CardMuted.TLabel").pack(anchor="w", padx=10, pady=8)

    def geo_search(self):
        q, only = self.geo_query.get().strip(), bool(self.geo_only_files.get())
        self.geo_msg.set("Cerco su GEO...")
        self.status_var.set("Ricerca su GEO in corso...")

        def job():
            try:
                self.q.put(("geo_results", 0, None, search_geo_series(q, only)))
            except Exception as e:
                self.q.put(("geo_err", 0, None, f"Ricerca non riuscita: {e}"))
        threading.Thread(target=job, daemon=True).start()

    def geo_load(self):
        sel = self.geo_tree.selection()
        if not sel:
            messagebox.showinfo("Seleziona uno studio", "Clicca prima su una riga dell'elenco.")
            return
        gse, title = sel[0], self.geo_tree.set(sel[0], "title")
        self.geo_msg.set(f"Carico {gse}...")
        self.status_var.set(f"Carico i metadati di {gse}...")

        def job():
            try:
                df, display = fetch_series_samples(gse)
                try:
                    files = rank_counts_files(geo_list_dir(f"{geo_series_url(gse)}/suppl"))
                except Exception:
                    files = []
                self.q.put(("geo_series", 0, None, {"gse": gse, "title": title, "samples": df,
                                                     "display": display, "files": files}))
            except Exception as e:
                self.q.put(("geo_err", 0, None, f"Impossibile caricare {gse}: {e}"))
        threading.Thread(target=job, daemon=True).start()

    def _geo_show_results(self, payload):
        rows, total = payload
        self.geo_tree.delete(*self.geo_tree.get_children())
        for r in rows:
            self.geo_tree.insert("", "end", iid=r["gse"], values=(r["gse"], r["year"], r["n_samples"], r["title"]))
        if rows:
            self.geo_msg.set(f"{len(rows)} studi mostrati su {total} trovati. Selezionane uno "
                             "(doppio clic o 'Usa lo studio selezionato').")
        else:
            self.geo_msg.set("Nessuno studio trovato: prova un'altra parola chiave (in inglese).")
        self.status_var.set("Ricerca GEO completata.")

    def _geo_apply_series(self, p):
        df, display, files = p["samples"], p["display"], p["files"]
        meta_cols = [c for c in df.columns if c not in ("gsm", "title")]
        fcols = [c for c in meta_cols if 2 <= df[c].nunique() <= GEO_MAX_FILTER_VALUES]
        self.geo_series = {"gse": p["gse"], "title": p["title"], "samples": df, "display": display,
                           "meta_cols": meta_cols, "filter_cols": fcols}
        names = [n for n, _ in files]
        self.geo_file_cb["values"] = names
        self.geo_file_var.set(names[0] if names else "")
        self.geo_title_var.set(f"{p['gse']} — {p['title']}  ({len(df)} campioni)")
        note = f"{p['gse']} caricato."
        if not files:
            note += " ⚠ Nessun file di conteggi riconoscibile tra i supplementari."
        elif files[0][1] <= 0:
            note += " ⚠ I file supplementari sembrano normalizzati (FPKM/TPM) o non di conteggi: controlla la scelta."
        self.geo_msg.set(note)
        self.status_var.set(note)
        self._geo_build_filters()
        self._geo_refresh_filters()

    def _geo_build_filters(self):
        for w in self.geo_filters_frame.winfo_children():
            w.destroy()
        s = self.geo_series
        self.geo_sel = {c: None for c in s["filter_cols"]}
        self.geo_vars, self.geo_boxes, self.geo_maps = {}, {}, {}
        if not s["filter_cols"]:
            ttk.Label(self.geo_filters_frame, style="CardMuted.TLabel",
                      text="Nessuna caratteristica utile per filtrare: puoi aggiungere tutti i campioni "
                           "dello studio come un unico gruppo.").pack(anchor="w", padx=10, pady=8)
            return
        grid = ttk.Frame(self.geo_filters_frame, style="Card.TFrame")
        grid.pack(fill=tk.X, padx=8)
        grid.columnconfigure(0, weight=1)
        grid.columnconfigure(1, weight=1)
        for i, col in enumerate(s["filter_cols"]):
            r, c = divmod(i, 2)
            cell = ttk.Frame(grid, style="Card.TFrame")
            cell.grid(row=r, column=c, sticky="ew", padx=8, pady=6)
            ttk.Label(cell, text=s["display"].get(col, col), style="CardBold.TLabel").pack(anchor="w")
            line = ttk.Frame(cell, style="Card.TFrame")
            line.pack(fill=tk.X, pady=(2, 0))
            var = tk.StringVar(value=ANY)
            cb = ttk.Combobox(line, textvariable=var, values=[ANY], state="readonly")
            cb.pack(side=tk.LEFT, fill=tk.X, expand=True)
            cb.bind("<<ComboboxSelected>>", lambda e, k=col: self._geo_on_select(k))
            self._no_wheel(cb)
            clr = ttk.Button(line, text="✕", style="Clear.TButton", width=3,
                             command=lambda k=col: self._geo_clear_one(k))
            clr.pack(side=tk.LEFT, padx=(4, 0))
            add_tip(clr, "Azzera questo filtro")
            self.geo_vars[col], self.geo_boxes[col] = var, cb

    def _geo_on_select(self, col):
        label = self.geo_vars[col].get()
        self.geo_sel[col] = None if label == ANY else self.geo_maps[col].get(label)
        self._geo_refresh_filters()

    def _geo_clear_one(self, col):
        self.geo_sel[col] = None
        self.geo_vars[col].set(ANY)
        self._geo_refresh_filters()

    def _geo_subset(self, exclude=None):
        df = self.geo_series["samples"]
        mask = pd.Series(True, index=df.index)
        for col, val in self.geo_sel.items():
            if val is not None and col != exclude:
                mask &= df[col] == val
        return df[mask]

    def _geo_refresh_filters(self):
        s = self.geo_series
        if not s:
            self.total = 0
            self.total_var.set("—")
            return
        for _ in range(len(s["filter_cols"]) + 1):          # una scelta può diventare non valida
            changed = False
            for col in s["filter_cols"]:
                counts = self._geo_subset(exclude=col)[col].value_counts()
                mapping = {f"{k}  [{n}]": k for k, n in counts.items()}
                self.geo_maps[col] = mapping
                self.geo_boxes[col]["values"] = [ANY] + list(mapping)
                cur = self.geo_sel[col]
                if cur is None:
                    self.geo_vars[col].set(ANY)
                    continue
                lab = next((l for l, v in mapping.items() if v == cur), None)
                if lab:
                    self.geo_vars[col].set(lab)
                else:
                    self.geo_sel[col] = None
                    self.geo_vars[col].set(ANY)
                    changed = True
            if not changed:
                break
        n = len(self._geo_subset())
        self.total = n
        self.total_var.set(f"{n:,}".replace(",", "."))
        self.status_var.set(f"Campioni di {s['gse']} corrispondenti ai filtri: {n}")

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
        if self.mode.get() == "geo":
            for c in self.geo_sel:
                self.geo_sel[c] = None
                self.geo_vars[c].set(ANY)
            self._geo_refresh_filters()
            return
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
                elif kind == "geo_results":
                    self._geo_show_results(payload)
                elif kind == "geo_series":
                    self._geo_apply_series(payload)
                elif kind == "geo_err":
                    self.geo_msg.set("⚠ " + payload)
                    self.status_var.set(payload[:140])
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
            messagebox.showwarning("Nessun dato", "Nessun campione corrisponde ai filtri correnti.\n"
                                                  "Prova ad azzerare qualche filtro.")
            return
        extra = {}
        if self.mode.get() == "geo":
            s = self.geo_series
            if not s:
                messagebox.showwarning("Nessuno studio", "Carica prima uno studio (passo A).")
                return
            if not self.geo_file_var.get():
                messagebox.showwarning("File dei conteggi", "Per questo studio non c'è nessun file di "
                                       "conteggi riconosciuto: non si possono scaricare i dati.")
                return
            sub = self._geo_subset()
            chosen = {s["display"].get(c, c): v for c, v in self.geo_sel.items() if v}
            desc = f"{s['gse']}: " + ("; ".join(f"{k} = {v}" for k, v in chosen.items()) or "tutti i campioni")
            default_name = " | ".join(chosen.values()) or s["gse"]
            rows = [{"gsm": r["gsm"], "title": r["title"], "source_name": r["source_name"],
                     "meta": {c: r[c] for c in s["meta_cols"]}} for _, r in sub.iterrows()]
            extra = {"source": "geo", "gse": s["gse"], "file": self.geo_file_var.get(),
                     "rows": rows, "desc": desc}
        else:
            chosen = {f: v for f, v in self.sel.items() if v}
            default_name = " | ".join(chosen.values()) or "Tutti"
            extra = {"source": "gdc", "sel": chosen}
        name = self.cond_var.get().strip() or default_name
        try:
            n = int(self.n_var.get())
        except (tk.TclError, ValueError):
            messagebox.showwarning("Numero non valido", "Scrivi un numero intero in 'Quanti campioni'.")
            return
        n = max(1, min(n, self.total))
        self.groups.append({"cond": name, "disp": self.total, "n": n, **extra})
        self._redraw()
        self.cond_var.set("")
        self.status_var.set(f"Gruppo '{name}' aggiunto ({n} campioni). Puoi cambiare i filtri e "
                            "aggiungerne un altro, oppure premere 'Avanti'.")

    def _redraw(self):
        self.tree.delete(*self.tree.get_children())
        for i, g in enumerate(self.groups):
            desc = g.get("desc") or "; ".join(f"{LABEL[f]} = {v}" for f, v in g["sel"].items()) or "nessun filtro"
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
        target = self._download_geo if self.mode.get() == "geo" else self._download
        threading.Thread(target=target, args=(list(self.groups), outdir), daemon=True).start()

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

    def _download_geo(self, groups, outdir):
        log = lambda t: self.q.put(("log", 0, None, t))
        try:
            log("Origine: GEO (ratto). I conteggi sono quelli caricati dagli autori di ogni studio.")
            chosen, seen = [], set()
            for g in groups:
                taken = 0
                for row in g["rows"]:
                    if taken >= g["n"]:
                        break
                    key = (g["gse"], row["gsm"])
                    if key in seen:
                        continue
                    seen.add(key)
                    chosen.append({**row, "cond": g["cond"], "gse": g["gse"], "file": g["file"]})
                    taken += 1
                log(f"Gruppo '{g['cond']}': {taken} campioni di {g['gse']}")
            if not chosen:
                raise RuntimeError("Nessun campione selezionato.")

            batches = {}
            for c in chosen:
                batches.setdefault((c["gse"], c["file"]), []).append(c)
            mats, ok = [], []
            for i, ((gse, fname), items) in enumerate(batches.items(), 1):
                log(f"Studio {gse}: {len(items)} campioni")
                if not fname:
                    raise RuntimeError(f"{gse}: nessun file di conteggi scelto (passo 1, 'File dei conteggi').")
                try:
                    counts = fetch_series_counts(gse, fname, [(c["gsm"], c["title"]) for c in items], log)
                except Exception as e:
                    raise RuntimeError(f"{gse} / {fname}: {e}")
                got = [c for c in items if c["gsm"] in counts.columns]
                lost = [c["gsm"] for c in items if c["gsm"] not in counts.columns]
                if lost:
                    log(f"  ! {len(lost)} campioni non trovati nel file e scartati: {', '.join(lost[:8])}"
                        + (" ..." if len(lost) > 8 else ""))
                if not got:
                    raise RuntimeError(f"{gse}: nessun campione scelto è presente nel file {fname}.")
                mats.append(counts[[c["gsm"] for c in got]])
                ok += got
                self.q.put(("prog", i, len(batches), None))

            warnings = []
            if len(mats) > 1:
                common = mats[0].index
                for m in mats[1:]:
                    common = common.intersection(m.index)
                if len(common) < 0.5 * min(len(m) for m in mats):
                    raise RuntimeError(
                        "Gli studi scelti usano identificatori di gene diversi (es. Ensembl e simboli) "
                        f"e hanno solo {len(common)} geni in comune: non possono essere uniti. "
                        "Usa gruppi dello stesso studio.")
                log(f"Unisco {len(mats)} studi sui {len(common)} geni in comune.")
                warnings.append("i gruppi vengono da studi diversi: attenzione all'effetto batch")
                counts = pd.concat([m.loc[common] for m in mats], axis=1)
            else:
                counts = mats[0]
            counts = counts.fillna(0)
            by_gsm = {c["gsm"]: c for c in ok}
            ok = [by_gsm[g] for g in counts.columns]

            arr = counts.to_numpy(dtype=float)
            if (arr < 0).any():
                warnings.append("ci sono valori negativi (dati trasformati, non conteggi grezzi)")
            frac = float((np.mod(arr, 1) != 0).mean())
            if frac > 0.001:
                warnings.append(f"{frac:.1%} dei valori ha decimali: potrebbero essere FPKM/TPM o stime "
                                "RSEM, non conteggi grezzi")
            for w in warnings:
                log(f"  ⚠ {w}")
            log("Esempio di ID gene: " + ", ".join(map(str, counts.index[:3])))

            log("Unisco i dati e salvo i file...")
            n = len(ok)
            base_cols = ["sample_id", "original_barcode", "condition", "case_id", "sample_type", "file_id"]
            meta = pd.DataFrame({
                "sample_id": [f"S{i}" for i in range(1, n + 1)],
                "original_barcode": [c["gsm"] for c in ok],
                "condition": [c["cond"] for c in ok],
                "case_id": [c["title"] for c in ok],
                "sample_type": [c.get("source_name") or NA_LABEL for c in ok],
                "file_id": [f"{c['gse']}/{c['file']}" for c in ok],
                "series": [c["gse"] for c in ok],
            })
            extra = pd.DataFrame([c["meta"] for c in ok])
            kept = []
            for col in extra.columns:
                if col == "source_name":
                    continue
                s = extra[col].fillna(NA_LABEL)
                if (s == NA_LABEL).all():
                    continue
                name = col + "_geo" if col in base_cols + ["series"] else col
                meta[name] = s.values
                kept.append(name)
            counts.columns = meta["sample_id"].tolist()
            counts.index.name = "gene_id"

            counts.reset_index().to_csv(os.path.join(outdir, "counts_matrix_unito.tsv"),
                                        sep="\t", index=False)
            meta.to_csv(os.path.join(outdir, "metadata_unito.tsv"), sep="\t", index=False)
            txt = (f"Salvati in {outdir}:\n"
                   f"- counts_matrix_unito.tsv ({counts.shape[0]} geni × {n} sample)\n"
                   f"- metadata_unito.tsv (con le colonne: series"
                   + (f", {', '.join(kept)}" if kept else "") + ")")
            if warnings:
                txt += "\n\nAttenzione:\n- " + "\n- ".join(warnings)
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
