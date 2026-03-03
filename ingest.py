import time
from typing import List, Dict
from Bio import Entrez #Biopython'ın Entrez modülünü kullanarak PubMed API'sine erişim sağlanır

#PubMed API için mail girilecek
Entrez.email = "youremail@gmail.com"

#Makaleleri arayalım ve PMID'lerini alalım
def search_pubmed( query: str, max_results: int=50) -> List[str]:
    """PubMed'de query için en güncel PMID listesini döndürelim"""
    handle = Entrez.esearch(
        db="pubmed",
        term=query, 
        retmax=max_results, 
        sort="pub date"
    )
    result = Entrez.read(handle) #apiden gelen sonucu okuyarak gelen veriyi sözlük haline getirir
    handle.close()
    print("[INFO] PubMed total matches:", result.get("Count"))
    return result.get("IdList", []) #makale kiimlik numaralarını içeren liste döndürür

#Abstract ve diğer bilgileri alamk için PMID'leri kullanıp makale bilgilerini alalım
def fetch_abstracts(pmids: List[str], sleep_s: float=0.5) -> List[Dict]:
    """PMId listesine göre PubMed'den başlık ve özet (abstract) bilgileri çeker """
    records = [] # makale bilgilerini tutmak için boş liste

    for pmid in pmids: #her makaleyi çekiyoruz ve rettype ile sadece özet bilgisini çektik, xml formatında döndürdük
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="xml" )
            data = Entrez.read(handle) #apiden gelen sonucu okuyarak gelen veriyi sözlük haline getirir
            handle.close()

            article = data["PubmedArticle"][0]["MedlineCitation"]["Article"] #makale bilgilerini içeren sözlük
            title = str(article.get("ArticleTitle", "")).strip() #makale başlığını alır, eğer yoksa boş string döndürür
            if "Abstract" in article and "AbstractText" in article["Abstract"]:
                abstract = article["Abstract"]["AbstractText"] #makale özetini alır
                abstract_text = " ".join([str(x) for x in abstract]).strip() #özet metnini birleştirir ve temizler
            else:
                abstract_text = "" #eğer özet yoksa boş string döndürür
            records.append({
                "pmid":str(pmid).strip(), #makale kimlik numarasını string olarak ekler
                "title":title, #makale başlığını ekler
                "abstract":abstract_text #makale özetini ekler
            })
        except Exception as e:
            print(f"Error fetching PMID {pmid}: {e}") #hata durumunda hata mesajını yazdırır
        time.sleep(sleep_s) #API'ye aşırı yüklenmeyi önlemek için bekleme süresi

    return records #makale bilgilerini içeren liste döndürür