# Genomic-RAG for RARS1

This project is a Retrieval-Augmented Generation (RAG) system designed to extract the most up-to-date clinical and molecular information for the **RARS1** gene directly from scientific literature. 

## 🧬 Project Overview
This project implements a Retrieval-Augmented Generation (RAG) system that dynamically extracts the most up-to-date clinical and molecular information about the RARS1 gene from scientific literature.
Because genomic knowledge evolves rapidly, static LLM knowledge is insufficient. Therefore, the system queries PubMed in real time, processes medical abstracts, indexes them in a vector database, and uses an LLM to produce structured, citation-grounded answers.

**The system extracts:**

- Genetic variants (e.g., c.5A>G)
- Associated diseases
- Clinical phenotypes/symptoms

Each claim is linked to a valid PubMed ID (PMID).


## ⚙️ System Architecture

 ### 1) Dynamic Data Ingestion
 **Key design decisions:**
  - Multiple query strategies were used to maximize recall:
    - RARS1
    - arginyl-tRNA synthetase
    - ArgRS
    - RARS AND cytoplasmic
    - Up to 200 results per query were retrieved
    - Results were merged and deduplicated
    - Only RARS1-relevant papers were kept using keyword filtering
    - Final dataset capped at 50 most relevant abstracts
           
**Why this approach?**

Initial attempts using a single, narrow PubMed query returned a limited set of only ~12 papers, which did not meet the project's requirement of 20-50 documents. Simply increasing the `retmax` didn't work because the most relevant papers for RARS1 are often indexed under different synonyms or overlapping biological terms.

**To overcome this, I implemented a "Broad Fetch & Strict Filter" strategy in my code:**
- **Expanded Search Depth & Union Query:** Instead of a single search, I implemented four distinct query strategies `(e.g., RARS1, arginyl-tRNA synthetase, ArgRS, and RARS AND cytoplasmic)` in main.py. I increased the retrieval limit to 200 results per query to capture a much larger initial pool of candidates.
- **Heuristic Post-Filtering:** Since fetching 200 papers introduces noise, I developed a custom filtering logic. The system scans the combined results (union) and checks both the Title and Abstract against a specific keyword list (including arginyl-tRNA synthetase 1, argrs, etc.) to ensure high precision.
- **Optimal Dataset Capping:** After filtering, the system selects the top 50 most relevant abstracts for indexing in ChromaDB.
 

 
### 2) Handling PubMed API Rate Limits
  Rate limiting was handled by:
  - Adding delays between API calls
  - Fetching abstracts sequentially
  - Limiting batch size
This prevents exceeding NCBI request thresholds.


### 3) Knowledge Processing & Chunking
  Medical text contains structured mutation names that must not be broken.
  A custom chunking algorithm was implemented to ensure that:
  - HGVS variant names remain intact
    `(e.g., c.2T>C (p.Met1Thr))`
  - Chunk boundaries never split mutations

This satisfies the requirement to preserve variant integrity during chunking.


### 4) Vector Database
  **Vector Store:** ChromaDB is used to store and index high-dimensional embeddings.
  **Embedding Model:** `all-MiniLM-L6-v2` provides a high-performance balance between semantic accuracy in biological text and computational efficiency.
  Reasons for selection:
  - Lightweight and fast
  - Strong performance on scientific text
  - No API cost
  - Suitable for short abstracts

 
 ### 5) LLM Extraction Layer
  The system follows a strict extraction protocol:
  - **Structured Output:** Extracts Variants, Associated Diseases, and Phenotypes into a structured JSON format.
  - **Citation Requirement:** Every molecular or clinical claim MUST include a valid PMID.
    
`Model Choice: Ollama (Local LLM)`
- Initially, OpenAI API was attempted, but access limitations and errors prevented stable use.
- Therefore, the system was switched to `Ollama with Qwen2.5-3B`.
- **Advantages:**
    - Runs locally
    - No API costs
    - No rate limits
    - Suitable for reproducible evaluation


 ### 6) Hallucination Guardrails
  **Safety Guardrail:** If a user query falls outside the retrieved context (e.g., asking about unrelated diseases), the system is programmed to return "I do not know" rather than fabricating an answer.
  To prevent fabricated medical claims:
  - All outputs must cite valid PMIDs
  - Citations must appear in retrieved snippets
  - Invalid mutation formats are rejected



## 🧪 Evaluation & Validation
To verify system reliability, we performed Hallucination Tests with "trick" questions. These tests demonstrate the system's ability to stay within its genomic scope:
| Test Query | System Response | Result |
| :--- | :--- | :---: |
| **COVID-19 Symptoms?** | "I do not know based on the retrieved literature." | **PASS** ✅ |
| **Diabetes Treatment?** | "I do not know based on the retrieved literature." | **PASS** ✅ |
| **RARS2 Variants?** | Identified disease but refused to link RARS1 variants. | **PASS** ✅ |



## 📁 Project Structure

```text
GENOMIC-RAG-RARS1/
├── chroma_db/             # Vector database storage
├── eval_results.json      # Hallucination test results
├── ingest.py              # PubMed API & Data indexing
├── main.py                # The core RAG pipeline, LLM integration, and evaluation logic
├── requirements.txt       # Dependencies
└── README.md              # Documentation
```


## 🚀 Setup and Usage
**1. Installation**
Install the necessary dependencies:
```bash
pip install -r requirements.txt
```
### Prerequisites
- **Python 3.8+**
- **Ollama** (Running `qwen2.5:3b`)
- **Key Libraries:** `biopython`, `chromadb`, `sentence-transformers`, `langchain`

**2. Data Indexing**
Run the indexing command to fetch and process the latest RARS1 data:
```bash
python main.py --index
```

**3. Ask a Question**
Ask a molecular or clinical question:
```bash
python main.py "List the genetic variants of RARS1 and their associated clinical phenotypes and diseases with PMIDs."
```
Example Output (JSON):
The system will return a structured response like this:
```json
{
"variants": [
{
"name": "c.5A>G",
"associated_phenotypes": ["Hypomyelination", "Nystagmus"],
"associated_diseases": ["Pelizaeus-Merzbacher-like disease"],
"citations": ["PMID:33515434"]
}
]
}
```

**4. Running Evaluations**
Generate the hallucination proof file `(eval_results.json)`:
```bash
python main.py --eval
```


**👤 Developed by:** Tuğba Çağla EREN
    
