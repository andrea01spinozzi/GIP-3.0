# Documentazione Dataset StandHard per Test GIP

Questo set di dati è stato generato dal programma **StandHard (v1.1)** per testare le funzionalità di analisi della piattaforma **GIP**. Il dataset confronta campioni tumorali con campioni di tessuto sano.

---

## 1. Descrizione Generale del Dataset

* **Origine dei Dati**: Combinazione di **GTEx v11** (tessuti sani) e **GDC Data Release 46.0** (tessuti tumorali).
* **Totale Campioni**: 40 campioni divisi in due gruppi bilanciati:
  * **20 campioni di Osteosarcoma** (da GDC: ossa, articolazioni e cartilagini degli arti).
  * **20 campioni di Cellule scheletriche sane** (da GTEx v11: tessuto *Muscle - Skeletal*).
* **Totale Geni**: 55.751 geni analizzati (es. `ENSG00000000003`).
* **Data di Creazione**: 9 Ottobre 2026.

---

## 2. Descrizione dei Singoli File

### 📄 `manifest.json`
È il file di metadati principale che descrive l'intera configurazione e l'esito dei controlli di qualità.
* **Cosa contiene**:
  * Informazioni sul programma e sulle versioni dei database utilizzati.
  * Filtri di selezione e parametri dei gruppi (Osteosarcoma vs Sano).
  * Hash SHA256 e dimensioni dei file generati per la verifica dell'integrità.
  * **Risultato del Controllo Qualità (QC)** e dettaglio dei fattori di variabilità.
  * Informazioni sull'ambiente di esecuzione (Python 3.13.3, pandas 2.3.3, numpy 2.2.6).
  * Mappatura dei 40 campioni (`sample_id` da `S1` a `S40`) con i rispettivi identificativi e barcode originali.

---

### 📊 `counts_matrix_unito.tsv`
È la **matrice di conteggio dei dati di espressione genica (RNA-Seq)**.
* **Cosa contiene**:
  * Un file formattato in formato TSV (*Tab-Separated Values*).
  * **Righe**: Geni identificati tramite Ensembl ID (es. `ENSG00000000003`).
  * **Colonne**: I 40 campioni identificati come `S1`–`S40`.
  * **Valori**: I conteggi grezzi (raw read counts) di espressione genica per ogni gene in ciascun campione.

---

### 📋 `metadata_unito.tsv`
È il file contenente le **informazioni cliniche, demografiche e tecniche** collegate a ciascun campione.
* **Cosa contiene**:
  * Una tabella TSV dove ogni riga corrisponde a uno dei campioni (`S1`–`S40`).
  * Informazioni di annotazione come:
    * Condizione (`Osteosarcoma` vs `Cellule scheletriche sane`).
    * Progetto di origine (`GDC` o `GTEx`).
    * Dati demografici (Età, Sesso, Etnia).
    * Dati di qualità e preparazione del campione (RIN, tempo di ischemia, tipo di decesso Hardy, batch di isolamento/espressione).

---

### ⚠️ `qc_fattori.tsv` (o tabella del Controllo Qualità)
Rappresenta l'analisi dei **fattori di confondimento (effetto batch e confondenti clinici)**.

#### Esito QC: **`CRITICO`**
L'analisi automatica evidenzia che l'analisi di espressione differenziale (DEA) è fortemente influenzata da **effetti batch** e variabili sovrapposte:

| Fattore | Variabilità Spiegata (%) | Rischio | Note / Indicazione |
| :--- | :---: | :---: | :--- |
| **Condizione** | 82.65% | `CONDIZIONE` | Segnale biologico principale. |
| **Origine del dataset** | 82.65% | `CRITICO` | Coincide al 100% con la condizione (GDC = Tumore, GTEx = Sano). Impossibile separare la biologia dal bias di piattaforma. |
| **Progetto (GDC)** | 82.65% | `CRITICO` | Stesso problema dell'origine del dataset. |
| **Organo / sito primario** | 82.65% | `CRITICO` | Sovrapposto alla condizione. |
| **Età** | 72.33% | `CRITICO` | Forte associazione con la condizione (87%). |
| **Tipo di decesso (Hardy)** | 47.13% | `CRITICO` | Spiega molta variabilità; richiede correzione/covariata nel modello. |
| **Sesso, RIN, Ischemia, Etnia** | < 11% | `BASSO` | Nessun impatto critico rilevato. |

> 📌 **Avviso principale**: Poiché i dati tumorali (GDC) e sani (GTEx) sono stati sequenziati ed elaborati con pipeline differenti, le differenze di espressione genica rilevate rischiano di riflettere l'**effetto batch** e non la vera biologia del tumore.
