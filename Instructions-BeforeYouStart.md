## GWAS-QC-AI-Report

*A plain-English tour of what the project does, why it matters, and where it is heading*

---

### 1.  What problem are we solving?

Think of a **Genome-Wide Association Study (GWAS)** as a giant “Where’s Waldo?” hunt—except Waldo is a DNA position linked to disease risk, and the picture is millions of columns (DNA variants) by thousands of rows (people).
Before any real discovery can happen, scientists must

1. **Clean the picture** (remove blurry or missing bits).
2. **Run the search** (statistical tests).
3. **Summarise what they found** (plots, tables, plain-language notes).

Traditionally these steps are glued together with ad-hoc scripts on individual laptops or lab servers. Results are hard to reproduce, QC is manual, and junior researchers spend days writing the same reports over and over.

---

### 2.  What does *GWAS-QC-AI-Report* change?

| Old workflow                                           | With our project                                                      |
| ------------------------------------------------------ | --------------------------------------------------------------------- |
| Dozens of shell & R scripts scattered across machines. | **One reproducible pipeline** driven by Nextflow and Conda.           |
| Manual copy-paste of QC numbers into Word.             | **Automatic PDF/HTML reports** generated straight from JSON metrics.  |
| Human interpretation written from scratch each time.   | **GPT-4/mini** writes the plain-English commentary for you.           |
| Everyone sets up their own environment.                | **Single Conda file** installs the exact same packages for all users. |

(*Container and cloud deployment will arrive in the next milestone—see Road-map.*)

---

### 3.  A gentle step-by-step tour

| Stage                             | What happens (non-technical)                                                                              | Key tool                   | Passes **to next stage**               |
| --------------------------------- | --------------------------------------------------------------------------------------------------------- | -------------------------- | -------------------------------------- |
| **1. Data Ingest**                | Load raw DNA files (VCF/BGEN) and a traits spreadsheet into one big table.                                | Hail                       | A clean “matrix” dataset.              |
| **2. Quality Check (QC)**         | Flag flawed samples or bad variants—think spell-checker for DNA.                                          | Hail & PLINK 2             | `qc_metrics.json`.                     |
| **3. Visualise**                  | Draw Manhattan and QQ plots → figures for papers.                                                         | Matplotlib                 | Publication-ready PNGs.                |
| **4. Association Test**           | Fast statistics find variants linked to the trait.                                                        | SAIGE / BOLT-LMM (wrapped) | List of hits.                          |
| **5. Polygenic Risk Score (PRS)** | Roll many small-effect variants into one “credit-score” per person.                                       | PRS-CS                     | `prs_scores.tsv` & `prs_metrics.json`. |
| **6. AI-Written Report**          | GPT reads the metrics and writes a lay summary: *“10 samples failed QC … PRS explains 12 % of variance.”* | OpenAI GPT-4/mini          | `report_YYYYMMDD.pdf`.                 |
| **7. Continuous Integration**     | Every code change triggers a mini-run on toy data.                                                        | GitHub Actions             | Green badge = pipeline still works.    |

All steps are scheduled by **Nextflow**, the “conductor” that fires each task only when its inputs are ready.

---

### 4.  Why this matters

* **Reproducibility** – everyone runs the same code, same environment, same numbers.
* **Time saving** – a week of QC+report writing shrinks to one overnight run.
* **Accessibility** – labs without dedicated bioinformaticians still get pro-grade outputs.
* **Collaboration** – the clean JSON + PDF outputs drop neatly into clinical or AI pipelines.

---

### 5.  What’s inside the repository?

* `src/` – four focused Python scripts (QC, plots, PRS, report).
* `workflow/` – one Nextflow file describing the order of jobs.
* `environment.yml` – one Conda spec; `conda env create -f …` and you’re set.
* `tests/` – Pytest sanity checks; failures show up in the GitHub badge.
* `report_template/` – Markdown + Jinja2 layout the AI fills in.

---

### 6.  Quick-start demo (local machine)

```bash
# Clone & create environment
git clone https://github.com/francott7/GWAS-QC-AI-report.git
cd GWAS-QC-AI-report
conda env create -f environment.yml      #  ➜ env name: gwas-qc
conda activate gwas-qc

# Add your OpenAI key (optional for AI commentary)
echo '{ "openai_api_key": "sk-..."}' > config.json

# Run the toy dataset (local profile)
nextflow run workflow/main.nf \
  --vcf data/chr22_demo.vcf.gz \
  --gwas data/demo_gwas_stats.tsv \
  --outdir dist/

open dist/report_*.pdf
```

---

### 7.  Current limitations & roadmap

| Limitation in v0.1                 | Planned enhancement                                         |
| ---------------------------------- | ----------------------------------------------------------- |
| No Docker / Singularity image yet. | Provide official container & `-profile docker` flag (v0.2). |
| Local machine or HPC only.         | Terraform modules for AWS EMR / GCP Dataproc (v0.3).        |
| Single-ancestry PRS tested.        | Multi-ancestry PRS-CSx, Tractor pipeline.                   |
| PRS method fixed to PRS-CS.        | Plug-in architecture for LDpred2, Lassosum2.                |
| Minimal UI (PDF only).             | Streamlit dashboard for drag-and-drop runs.                 |

---

### 8.  Who should care?

* **Genetics labs** needing turnkey QC/PRS without writing code.
* **Data-science teams** wanting genetic risk variables inside models.
* **Clinicians** who need a one-page QC certificate before trusting a dataset.
* **Educators** showing students a full GWAS pipeline on laptops.

---

### 9.  Final take-away

*GWAS-QC-AI-Report* bridges modern engineering practice and cutting-edge genomics.
It turns messy DNA files into clear, reproducible insights—automatically and in language anyone can understand.
Containers and cloud scale are on the way, but you can already run the full workflow today with nothing more than Conda and Nextflow.

Fork, star ⭐, or raise an issue—contributions are welcome!
