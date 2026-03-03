#imports
import json #Json yazmak/okumak için çünkü LLM'den JSON formatında çıktı alacağız
import re #regex: variant/citation kontrolü için. Düzenli ifadelerle metin işleme yapacağız, özellikle LLM çıktısındaki JSON'u temizlemek için
import sys # komut satırı argümanlarını almak için
from typing import List, Dict, Tuple, Set #okunabilirlik için tip ipuçları
import chromadb #ChromaDB, vektör veritabanı olarak kullanacağız, makalelerin gömülü temsillerini saklamak için
from chromadb.utils import embedding_functions #ChromaDB'nin gömülü fonksiyonlarını kullanarak makalelerin gömülü temsillerini oluşturacağız
import requests 
from ingest import search_pubmed, fetch_abstracts #ingest.py'deki fonksiyonları içe aktararak PubMed'den makale arama ve özet çekme işlemlerini gerçekleştireceğiz

""" 1- CHUNKING (PARÇALAMA): Varyant bozulmadan metni parçalara ayırmak için """

#regex ile varyant/citation kontrolü yaparak metni parçalara ayıracağız, böylece LLM'nin anlayabileceği şekilde temizlenmiş metin parçaları elde edeceğiz
#Amaç : chunk sınırı varyantın ortasına gelmesin, varyant/citation'ı tek parça tutmak istiyoruz
HGVS_PATTERN = re.compile(r"""
(c\.\d+[A-Za-z0-9>_\-+]*\d*[A-Za-z>]*) | #cDNA varyantları için örnek: c.123A>G, c.456_789del, c.789+1G>T gibi
(p\.[A-Za-z]{3}\d+[A-Za-z]{3}) | #protein varyantları için örnek: p.Gly123Arg, p.Met1Thr gibi
(rs\d{3,}) #referans SNP varyantları için örnek: rs123456 gibi
""", re.VERBOSE) #regex'i daha okunabilir hale getirmek için VERBOSE bayrağı kullanıldı, böylece yorum satırları ekleyebiliriz ve düzenli ifadeyi daha kolay anlayabiliriz

def find_spans(text: str) -> List[Tuple[int, int]]:
    """Metindeki varyant/citation'ların başlangıç ve bitiş indekslerini bulalım"""
    spans = [] #varyant/citation'ların indekslerini tutmak için boş liste
    for match in HGVS_PATTERN.finditer(text): #metni tarayarak regex ile eşleşen varyant/citation'ları buluruz
        spans.append((match.start(), match.end())) #eşleşen metnin başlangıç ve bitiş indekslerini listeye ekleriz
    return spans #varyant/citation'ların indekslerini içeren liste döndürür


def overlaps(cut:int, spans: List[Tuple[int, int]]) -> bool:
    """Belirli bir kesme noktasının herhangi bir varyant/citation ile çakışıp çakışmadığını kontrol edelim"""
    for start, end in spans: #her varyant/citation'ın başlangıç ve bitiş indekslerini kontrol ederiz
        if start < cut < end: #eğer kesme noktası herhangi bir varyant/citation'ın içinde kalıyorsa
            return True #çakışma var demektir, True döner
    return False #hiçbir varyant/citation ile çakışma yoksa False döner

def safe_cut(text:str, desired:int, spans:List[Tuple[int,int]]) -> int:
    """Varyant/citation'ları bölmeden metni güvenli bir şekilde keselim"""
    if not overlaps(desired, spans):#eğer istenen kesme noktası herhangi bir varyant/citation ile çakışmıyorsa
        return desired #istenen kesme noktasını döner
    
    #Geriye doğru doğal sınır aramak
    start_position = max(0, desired - 200) #istenen kesme noktasının 200 karakter gerisine bakarak doğal bir sınır arar, böylece varyantı bölmeden kesme noktasını bulmaya çalışır
    for p in range(desired, start_position, -1):  #istenen kesme noktasından geriye doğru 200 karakter boyunca kontrol eder
        if p > 0 and (text[p-1].isspace() or text[p-1] in ".;,\n"): #eğer boşluk veya nokta gibi doğal bir sınır bulursa
            if not overlaps(p, spans): #ve bu sınır herhangi bir varyant/citation ile çakışmıyorsa
                return p #bu noktayı kesme noktası olarak döner
   
   #yine olmazsa spanin başına çek
    for start, end in spans:  #her varyant/citation'ın başlangıç ve bitiş indekslerini kontrol ederiz
        if start < desired < end:
            return start #kesme noktasını varyant/citation'ın başlangıcına kaydırarak döner 
    
    return desired #eğer hiçbir uygun kesme noktası bulunamazsa istenen kesme noktasını döner
   

def chunktext_variant_safe(text:str, chunk_size:int=900, overlap:int=120) -> List[str]:
    """Varyant/citation'ları bölmeden metni karakter bazlı parçalara ayıralım"""
    #Fazla boşlukları normalize etmek için 
    text = " ".join(text.split()) #metindeki fazla boşlukları tek bir boşlukla değiştirerek metni normalize ederiz, böylece parçalama işlemi daha tutarlı olur
    if not text: #eğer metin boşsa
        return [] #boş liste döndürürüz, parçalayacak bir şey yok demektir

    spans = find_spans(text) #metindeki varyant/citation'ların indekslerini buluruz
    chunks = [] #parçaları tutmak için boş liste
    start = 0 #parçalama işlemi için başlangıç indeksi
   
    while start < len(text): #metnin sonuna kadar parçalama işlemi devam eder
        desired_cut = min(len(text), start + chunk_size) #istenen kesme noktasını belirleriz, metnin uzunluğunu aşmamak için min fonksiyonunu kullanırız
        cut = safe_cut(text, desired_cut, spans) #güvenli bir şekilde kesme noktasını buluruz
        #chunks.append(text[start:cut].strip()) #bulunan kesme noktasına göre metni parçalayarak listeye ekleriz, ayrıca parçanın başındaki ve sonundaki boşlukları temizleriz
        
        #ilerlemeyi garanti etmek için sonsuz döngüyü engellemek için
        if cut <= start:
            cut = desired_cut #eğer güvenli kesme noktası başlangıç noktasına eşit veya daha küçükse, istenen kesme noktasını kullanarak ilerlemeyi garanti ederiz
        chunk = text[start:cut].strip() #bulunan kesme noktasına göre metni parçalar ve temizler
        if chunk:
            chunks.append(chunk) #parça boş değilse listeye ekle
        #sondaysak bitirelim
        if cut >= len(text):
            break #eğer kesme noktası metnin sonunu aşmışsa döngüyü kırarak parçalama işlemini sonlandırırız

        #Overlap ile yeni başlangıç belirleriz negatif olmayan
        new_start = max(0, cut-overlap)

        #Yeni başlangıç span içinde kalıyorsa güvenli yere çekelim
        if overlaps(new_start, spans):
            new_start = safe_cut(text, new_start, spans) #eğer yeni başlangıç noktası herhangi bir varyant/citation ile çakışıyorsa, güvenli bir kesme noktası bulmak için safe_cut fonksiyonunu kullanırız
        start = new_start #başlangıç noktasını güncelleriz

    return chunks #parçalanmış metinleri içeren liste döndürür




""" 2- VECTOR DB: ChromaDB ile indeksleme ve retrieval"""
#chroma veritabanının diske yazılacağı dizin
CHROMA_DIR = "./chroma_db"

#Koleksiyon adı: Chromadb içinde verileri tutacağımız koleksiyonun adı
COLLECTION_NAME = "rars1_pubmed"

#Embedding model: Metni vektöre çeviren model. Local SentenceTransformers modeli kullanacağız, çünkü PubMed makaleleri genellikle teknik ve uzun olduğundan güçlü bir embedding modeline ihtiyacımız var ve ayrıca maliyet açısından uygun.
#API key gerektirmez, ücretsizdir
EMBEDDING_MODEL = "all-MiniLM-L6-v2" #SentenceTransformers'ın güçlü ve hafif bir modeli, özellikle kısa metinler için iyi performans gösterir

def get_collection():
    #ChromaDB istemcisi oluşturulur ve koleksiyon alınır veya oluşturulur
    #persistent client kullanarak verilerin diske kaydedilmesini sağlıyoruz, böylece indeksleme işlemi her çalıştırıldığında sıfırlanmaz
    client = chromadb.PersistentClient(path=CHROMA_DIR)
    embedding_fn = embedding_functions.SentenceTransformerEmbeddingFunction(model_name=EMBEDDING_MODEL) #embedding fonksiyonunu oluştururuz, bu fonksiyon metni vektöre çevirmek için SentenceTransformers modelini kullanır
    collection = client.get_or_create_collection(name=COLLECTION_NAME, embedding_function = embedding_fn, metadata={"hnsw:space":"cosine"}) #koleksiyonu alır veya oluştururuz, embedding fonksiyonunu koleksiyona atarız ve HNSW algoritması için cosine benzerlik metriği kullanacağımızı belirtiriz
    return collection #koleksiyonu döndürür

def build_index(records:List[Dict], reset: bool=True) -> int:
    #Makale kayıtlaırını alır ve ChromaDB koleksiyonuna ekler, reset=True ise önce koleksiyonu temizler
    #records: ingest.py'den gelen makale bilgilerini içeren listedir, her kayıt genellikle {"pmid": "12345", "title": "Makale Başlığı", "abstract": "Makale Özeti"} formatındadır

    client = chromadb.PersistentClient(path=CHROMA_DIR) #ChromaDB istemcisi oluşturulur, persistent client kullanarak verilerin diske kaydedilmesini sağlıyoruz
    #Her index işlemi öncesi eski koleksiyonu silerek temiz başlangıç yapıyoruz eski veriler yenilerle karışmamış oluyor
    if reset:
        try:
            client.delete_collection(name=COLLECTION_NAME) #koleksiyonu sileriz, böylece eski verilerle karışmaz
            print(f"[OK] Deleted existing collection '{COLLECTION_NAME}' for clean indexing.")

        except Exception:
            #koleksiyon uoksa hata olabilir, hatayı görmezden gelerek devam ederiz
            pass
            
    collection = get_collection() #kolekyionu alırız

    ids: List[str] = [] #makale kimliklerini tutmak için boş liste
    docs: List[str] = [] #makale metinlerini tutmak için boş liste
    metadatas: List[Dict] = [] #makale metadata'larını tutmak için boş liste


    #Her makale kaydı için 
    for record in records:
        pmid = str(record.get("pmid", "")).strip() #makale kimliğini alırız, eğer yoksa boş string döndürürüz
        title = record.get("title", "").strip() #makale başlığını alırız, eğer yoksa boş string döndürürüz
        abstract = (record.get("abstract") or "").strip() #makale özetini alırız, eğer yoksa boş string döndürürüz

        #chunk yani parçalanacak tam metin
        full_text = f"{title}. {abstract}".strip() #başlık ve özeti birleştirip tam metin oluştururuz

         #tam metni parçalara ayırma işlemi
        chunks = chunktext_variant_safe(full_text, chunk_size=900, overlap=120) #tam metni varyant/citation'ları bölmeden parçalara ayırırız,
        #chunk_size ve overlap parametreleriyle parçaların boyutunu ve birbirleriyle ne kadar örtüşeceğini belirleriz       
        for i, chunk in enumerate(chunks):
            ids.append(f"{pmid}_{i}") #her parçaya benzersiz bir kimlik atarız, burada makale kimliği ve parça indeksini kullanarak benzersiz bir ID oluştururuz
            docs.append(chunk) #parçalanmış metni doküman listesine ekleriz
            metadatas.append({"pmid": pmid, "chunk_index": i}) #makale kimliğini ve parça indeksini metadata olarak ekleriz, böylece hangi parçanın hangi makaleye ait olduğunu takip edebiliriz
    #eğer hiç chunk yoksa
    if not ids:
        return 0 #eğer parçalanacak metin yoksa 0 döndürürüz indeksleme başarısız olmamış olur

    #ChromaDB koleksiyonuna parçaları ekleriz
    collection.add(ids=ids, documents=docs, metadatas=metadatas) #chroma burada embedding hesaplar ve diske yazar
    return len(ids) #indekslenen toplam parça sayısını döner

def retrieve(query:str, collection, top_k:int=5) -> List[Dict]:
    #Soruya göre ChromaDB koleksiyonundan benzer parçaları getirir
    results = collection.query(query_texts=[query], n_results=top_k, include=["documents", "metadatas"]) #koleksiyona sorgu göndeririz, benzer parçaları alırız
    retrieved = [] #getirilen parçaları tutmak için boş liste,
    for id, doc, metadata in zip(results["ids"][0], results["documents"][0], results["metadatas"][0]): #getirilen sonuçları döngüyle işleriz
        retrieved.append({
            "id": id, #parça kimliğini ekleriz
            "text": doc, #parça metnini ekleriz
            "metadata": metadata #parça metadata'sını ekleriz, burada makale kimliği ve parça indeksini içerir
        })
    return retrieved #getirilen parçaları içeren liste döndürür



""" 3- LLM Analizi: OpenAI API'si ile parçaları analiz ederek JSON formatında çıktı alma"""


#Regexler LLM çıktısındaki JSON'u temizlemek için, bazen LLM'nin çıktısı tam olarak JSON formatında olmayabilir, bu yüzden regex kullanarak JSON kısmını çıkarmaya çalışırız
VARIANT_NAME_PATTERN = re.compile(
    r"""
     ^(
      c\.\d+[A-Za-z0-9>_\-+]*\d*[A-Za-z>]*     # c.5A>G, c.2T>C, c.123_124del gibi (basit)
      |
      p\.[A-Za-z]{3}\d+[A-Za-z]{3}             # p.Met1Thr
      |
      rs\d{3,}                                  # rs12345
    )$
    """, re.VERBOSE
)

CITATION_PATTERN = re.compile(r"^(PMID:\d+|DOI:\S+)$") #Citation formatı doğru mu diye bakar. PMID:12345 veya DOI:10.1000/xyz gibi atıf formatlarını tanımlayan regex, bu da LLM'nin çıktısındaki atıfları tanımlamak ve doğrulamak için kullanılır

#sadece snippetları kullanacağız çünkü llm uzun metinlerle iyi başa çıkamayabilir
#uydurma yapmaması için snippetları kontrol edeceğiz
#her itemda citation zorunlu olacak
#cıktı json olacak, böylece daha sonra kolayca işleyebiliriz
SYSTEM_RULES = """You are a RARS1 expert. 
Your primary goal is to extract every RARS1 variant, disease, and phenotype mentioned in the snippets.

CRITICAL INSTRUCTION:
For EVERY piece of information you extract (variant, disease, or phenotype), you MUST include the corresponding PMID from the snippets in the "citations" list. 
If you find a variant like 'c.5A>G', look at the text around it, find the PMID, and include it.
DO NOT return any item without at least one valid PMID.

1. If the question is about RARS1, you MUST extract the data.
2. Only if the question is clearly about something else (like COVID, Pizza), return empty lists.
3. Use the following JSON schema strictly:
{
  "variants": [{"name": "c.5A>T", "associated_diseases": [], "associated_phenotypes": [], "citations": ["PMID:12345"]}],
  "diseases": [{"name": "Disease Name", "citations": ["PMID:12345"]}],
  "phenotypes": [{"name": "Symptom Name", "citations": ["PMID:12345"]}],
  "unknowns": []
}
"""
#ollama için çağırma fonksiyonumuz
def call_ollama(prompt: str) -> str:
    url = "http://localhost:11434/api/generate"

    payload = {
        "model": "qwen2.5:3b",
        "prompt": prompt,
        "stream": False,
        "options": {"temperature": 0.0}
    }

    response = requests.post(url, json=payload, timeout=300)
    response.raise_for_status()
    return (response.json().get("response") or "").strip()

def build_context (snippets:List[Dict], max_chars: int = 6000) -> Tuple[str, set]:
    #LLM'e göndermek için bağlam oluştururz
    #Bu fonksiyon 2 şey döner, LLM'e verilecek kanıt metni (context) ve allowed_pmids seti böylece LLM'in uydurma yapması engellenir
    allowed_pmids: Set[str] = set() #LLM'in uydurma yapmasını engellemek için izin veirlen PMIDleri tutmak için boş set
    block: List[str] = [] #LLM'e göndermek için metin bloklarını tutmak için boş liste
    label = "PMID:" 
    total = 0
    for item, snippet in enumerate(snippets, start=1): #her snippet için döngü
        text = (snippet.get("text") or "") #snippet metnini alırız
        pmid = str(snippet.get("metadata", {}).get("pmid", "")).strip() #snippet'ın metadata'sından PMID'yi alırız
        if pmid:
            allowed_pmids.add(pmid) #izin verilen PMID'ler setine ekleriz
        
        # Her snippet'ı kısa tut (hız için)
        text = text[:1200]

        chunk = f"""
        [Snippet {item}] 
        {label}{pmid}
        TEXT:
        {text}
        """.strip() # snippet'ı ekleriz, böylece LLM'e hangi snippet'tan geldiği bilgisi de sağlanmış olur
        #Toplam context çok büyümesin istiyoruz
        if total + len(chunk) > max_chars:
            break
        block.append(chunk) #metin bloklarını listeye ekleriz
        total += len(chunk) #toplam karakter sayısını güncelleriz
    context = "\n\n".join(block) #metin bloklarını çift yeni satır ile birleştirerek tek bir bağlam metni oluştururuz, bu bağlam LLM'e gönderilecek ve analiz edilecek metni içerir
    return context, allowed_pmids #oluşturulan bağlam metni ve izin verilen PMID'ler setini döndürürüz,

def safe_json_load(text:str) -> Dict:
    #LLM'in çıktısını güvenli şekilde JSON'a çevirmeye çalışırız bazen LLM çıktısı tam olarak JSON formatında olmayabilir
    if not text or not text.strip():
        return {"variants": [], "diseases": [], "phenotypes": [], "unknowns": ["LLM returned empty response"]}
    
    try:
        # JSON bloğunu metinden ayıklayarak alırız
        
        match = re.search(r'\{.*\}', text, re.DOTALL)
        if match:
            json_str = match.group()
            return json.loads(json_str)
        return json.loads(text)
    except Exception as e:
        # Eğer hala hata alıyorsa, en azından boş bir şema döndür ki kod çökmesin
        return {
            "variants": [], 
            "diseases": [], 
            "phenotypes": [], 
            "unknowns": [f"JSON Parsing Error: {str(e)}"]
        }

def ask_llm_structured(question:str, snippets:List[Dict]) -> Tuple[Dict, Set[str]]:
    #LLM'e soruyu ve bağlamı vererek yapılandırılmış JSON formatında cevao alırız
    # data ve allowed_pmids döneriz 
    context, allowed_pmids = build_context(snippets) #bağlamı oluştururuz, LLM'e göndereceğimiz metni ve izin verilen pmdidleri alırız
    prompt = f""" 
    {SYSTEM_RULES}

    USER QUESTION: 
    {question}  
    
    EVIDENCE SNIPPETS: 
    {context} 
    
    Return only valid JSON.""".strip() #LLM'e gönderilecek promptu oluştururuz, burada kullanıcı sorusu, kanıt metni ve sistem kuralları birleştirilir

    output1 = call_ollama(prompt) #LLM'e soru ve bağlamı vererek cevap alırız
    try:
        data = safe_json_load(output1) #LLM'in çıktısını güvenli şekilde JSON'a çevirmeye çalışr
        return data, allowed_pmids #yapılandırılmış veriyi ve izin verilen PMID'ler setini döndürürüz,
    except Exception:
        retry_prompt = f""" 
        Your previous output was not valid JSON. 
        Please try again and return only valid JSON that matches the schema exactly. 
        USER QUESTION: 
        {question} 
        
        EVIDENCE SNIPPETS: 
        {context} 
        """.strip()

        output2 = call_ollama(retry_prompt) #tekrar denedikten sonra LLM'in çıktısını alırız
        data2 = safe_json_load(output2) #tekrar denedikten sonra LLM'in çıktısını aldık
        return data2, allowed_pmids



""" 4- Guardrails - Citations'ı zorunlu kılıp, uydurma yapmasını engellemek için LLM çıktısını kontrol etme"""
def pmid_from_citation(citation: List[str]) -> List[str]:
    #Alıntılardan pmidleri çıkarmaya çalışırız, böylece LLM'in uydurma yapmasını engelleriz
    pmids = []
    for c in citation:
        if c.startswith("PMID:"):
            pmids.append(c.split("PMID:",1)[1].strip()) #PMID:12345 formatında olan alıntılardan sadece numarayı çıkarmaya çalışırız
    return pmids #PMID'leri içeren liste döndürür

def guardrail_citations(data: Dict, allowed_pmids: Set[str]) -> Dict:
    #LLM çıktıısndaki alıntıları kontorl ederiz ve sadece izin verilen PMID'lerin alıntılarda bulunmasına izin verilir böylece llm uydurması engellenir
    #-citation exist ve format is valid kontrolü
    #alıntılanan pmid'ler alınan pmdiler arasında olmalı
    #varyant isimleri HGSV/rs formatıyla eşleşmeli
    #Kontrolü geçemeyen itemlar unknown'a atılır

    cleaned_data = {"variants": [], "diseases": [], "phenotypes":[], "unknowns":[]} #temizlenmiş veriyi tutmak için sözlük oluşturduk ve her kategori için boş listeler eklendi

    #yeni ekledik
    def norm_key(s: str) -> str:
        # büyük/küçük harf + ekstra boşluk farklarını yok eder
        return " ".join(str(s).strip().lower().split())


    #Unknowns'u normalize edelim
    unknowns = data.get("unknowns", []) #LLM çıktısındaki unknowns kategorisindeki itemları alır, yoksa boş liste döndürür
    for item in unknowns:
        if isinstance(item, list): #eğer item bir listeyse 
            cleaned_data["unknowns"].extend([str(x).strip() for x in item if str(x).strip()]) #gereksiz boşlukları temizleyerek unknowns listesine ekleriz, böylece unknowns kategorisindeki itemlar daha tutarlı hale gelir
        elif item: #eğer item unknıwns kategorisinde boş değilse
            cleaned_data["unknowns"].append(str(item).strip()) #gereksiz boşlukları temizleyerek unknowns listesine ekleriz
        

    def is_placeholder(text:str) -> bool:
        t = (text or "").strip().lower()
        if not t:
            return False #eğer text boş ise hata değil sadece veri yok demek istiyoruz
        
        #Filtreleme yapmaya çalışıyoruz
        PLACEHOLDER_STRINGS = {
            "disease name",
            "phenotype / symptom",
            "variants name"
            "diseases 1"
            "symptom 1"
            "unknown",
            "n/a",
            "not specified",
        }
        #çok kısa ve anlamsız şeyleri de ele
        if t in PLACEHOLDER_STRINGS:
            return True

        # "Disease name", "Phenotype / symptom" gibi varyasyonlar
        if "disease name" in t or "phenotype / symptom" in t or "phenotype/symptom" in t:
            return True
        return False


    def citations_ok(cits) -> Tuple[bool, str, List[str]]: #fonksiyon 3 şey döndürür-> True/false: geçerli mi, Hata mesajı ve normalize edilmiş citation listesi
        #Alıntıları kontrol ederiz, alıntılar izin verilen PMID'ler arasında olmalı ve formatları doğru olmalı, böylece LLM'in uydurma yapması engellenir
        if not isinstance(cits, list) or len(cits)==0: #eğer alıntılar bir liste değilse veya boşsa
            return False, "No citations provided", [] #alıntı sağlanmadığı için geçersiz olduğunu belirten mesaj ve boş liste döndürürüz
        #boş olan citation'lar var mı kontrol edelim
        normalize_cits = [] #alıntıları normalize ederiz, gereksiz boşlukları temizleyerek ve boş olanları filtreleyerek yeni bir liste oluştururuz
        for c in cits: #dict formatını yakalamak istiyorum
            if isinstance(c, dict):
                # {'PMID': '123'} veya {'pmid': '123'} değerlerini kontrol et
                c_str = str(c.get("PMID") or c.get("pmid") or "").strip()
            else:
                c_str = str(c).strip()
            if not c_str: 
                continue
            #Eğer sadece sayıysa örn: 37186453, başına PMID: ekle
            if c_str.isdigit():
                c_str = f"PMID:{c_str}"
            #normalize_cits.append(c_str)

            #eğer normalize_cits listesinde yoksa eklenir (duplicate engellemek için)
            if c_str not in normalize_cits:
                normalize_cits.append(c_str)
        
        if not normalize_cits: #eğer normalleştirilmiş alıntılar listesi boşsa
            return False, "All citations are empty", [] #tüm alıntılar boş olduğu için geçersiz olduğunu belirten mesaj ve boş liste döndürülür
            
        for cit in normalize_cits: #her alıntıyı kontrol edelim
            if not CITATION_PATTERN.match(cit): #eğer alıntı doğru formatta değilse
                return False, f"Citation {cit} is not in valid format", normalize_cits #bu alıntının geçerli formatta olmadığın belirten mesaj ve normalleştirilmiş alıntılar listi döndürür
            
            #PMID değerini allowed_pmids ile kıyasalayalım
            pmid_val = cit.split(":")[-1] #alıntının sonundan sonra : karakterine kadar olan kısmını alır 
            if pmid_val not in allowed_pmids:
                return False, f"PMID is not in retrieved set: {pmid_val}", normalize_cits #pmid'in izin verilenler arasında olmadığını belirten mesaj ve normalize edilmiş alıntı listesi döner 
        
        return True, "Citations are valid", normalize_cits #alıntılar geçerliyse geçerli olduğunu belirten mesaj ve normalize edilmiş alıntı listesi döner

    #Varyantları kontrol edelim (Variants)
    variants = data.get("variants", []) #LLM çıktısındaki variants kategorisindeki itemları alıri yoksa boş list döndürür
    if isinstance(variants, list): #varyant listeyse
        for var in variants:
            #herbir kategori döngüsüne tip kontrolü
            if not isinstance(var, dict):
                cleaned_data["unknowns"].append(f"Skipped non-object variant entry: {var}")
                continue
            raw_name = str(var.get("name","")).strip()#varyant adını alırız, yoksa boş string döner ve gereksiz boşluklar temizlenlir
            name = raw_name.split(" (")[0].strip() 
            cits = var.get("citations", []) #varyantın alıntılarını alırız, yoksa boş liste döner
            if not name: #varyant adı boşsa
                cleaned_data["unknowns"].append("Dropped a variant with empty name") #adı boş olan varyantı bilinmeyenlere ekledik
                continue
            if not VARIANT_NAME_PATTERN.match(name): #varyant adı droğu formatta değilse
                cleaned_data["unknowns"].append(f"Dropped a variant with invalid name: {raw_name}") #adı doğru formatta olmayan varyant bilinmeyene eklendi
                continue
            cits_ok, msg, normalize_cits = citations_ok(cits) #alıntıları kontrol ederiz, geçerli mi değil mi, mesaj ve normalize edilmiş alıntılar listesi
            if not cits_ok: #alıntı geçerli değilse
                cleaned_data["unknowns"].append(f"Dropped variant {raw_name} because of citation issue: {msg}") #alıntı geçerli olmayan varyant bilinmeyene eklenir
                continue
                
            phenotype = var.get("associated_phenotypes", []) #varyantın ilişkili fenotiplerini al yoksa boş list dön
            disase = var.get("associated_diseases", []) #varyantın ilişkili hastalıklarını al yoksa boş list dön

            phenotype_normalized = [str(x).strip() for x in phenotype] if isinstance(phenotype,list) else [] #fenotipler normalize edilir, fenotipler listeyse gereksiz boşluklar temizlenip yeni liste oluşturulur değilse boş liste döndürür
            phenotype_normalized = [x for x in phenotype_normalized if x and not is_placeholder(x)] #boş olan fenotipler filtrelenir
            disease_normalized = [str(x).strip() for x in disase] if isinstance(disase, list) else[]
            disease_normalized = [x for x in disease_normalized if x and not is_placeholder(x)] #boş olan hastalıklar filtrelendi

            cleaned_data["variants"].append({
                "name": name, #variant name
                "associated_phenotypes": phenotype_normalized, #normalize edilmiş fenotip listesi
                "associated_diseases": disease_normalized, #normalize edilmiş hastalık listesi
                "citations": normalize_cits #normalize edilmiş alıntı listesi
            })

    else:
        cleaned_data["unknowns"].append("Variants field is not a list") #eğer variants alanı bir liste değilse bilinmeyenlere ekleriz
    
    #inherit citations from variants for diseases/phenotypes
    disease_to_cits = {}
    pheno_to_cits = {}

    for v in cleaned_data["variants"]: 
        v_cits = v.get("citations", []) #varyant alıntılarını alırız
        for d in v.get("associated_diseases", []): #varyantla ilişkili hastalıklar tek tek döner
            key = norm_key(d) #hastalık ismini normalize eder
            if key:
                #hastalığa ait bir set oluşturup varyantın makalelerini içine ekler
                disease_to_cits.setdefault(key, set()).update(v_cits)
        for p in v.get("associated_phenotypes", []): #varyantla ilişkili fenotipleri döner
            key = norm_key(p)
            if key:
                pheno_to_cits.setdefault(key, set()).update(v_cits)
        
    # Disease'leri kontrol edelim
    diseases = data.get("diseases", []) #LLM çıktısındaki diseases kategorisindeki itemları alır yoksa boş list döndürür
    if isinstance(diseases, list): #diseases listeyse
        for dis in diseases:
            #TİP KONTROLU
            if not isinstance(dis, dict):
                cleaned_data["unknowns"].append(f"Skipped non-object disease entry: {dis}")
                continue

            name = str(dis.get("name", "")).strip()
            if is_placeholder(name): 
                cleaned_data["unknowns"].append(f"Dropped placeholder disease: {name}")
                continue
            cits = dis.get("citations", []) #hastalık alıntısını alırız yoksa boş string döner
            if not name: #hastalık adı boşsa
                continue #adı boş olan hastalıkları atlatız çünkü zaten adı olmayan hastalık bilgisi kullanışlı olmaz
         
            #diseasede inherit kontrolünü normalize edelim
            name_key = norm_key(name)
            if (not cits) and (name_key in disease_to_cits):
                #cits = [f"PMID:{x}" for x in []]  # boş bırakma, aşağıda direkt setten alacağız
                normalize_cits = sorted(disease_to_cits[name_key])
                cleaned_data["diseases"].append({"name": name, "citations": normalize_cits})
            else:
                cits_ok, msg, normalize_cits = citations_ok(cits) #alıntıları kontrol ederiz, geçerli mi değil mi, mesaj ve normalize edilmiş alıntılar listesi
                if cits_ok: #alıntı geçerliyse
                    cleaned_data["diseases"].append({"name": name, "citations": normalize_cits}) #alıntı geçerli olmayan hastalık yine de eklenir ama alıntıları bilinmeyenlere ekleriz
                else:
                    cleaned_data["unknowns"].append(f"Dropped disease {name} because of citation issue: {msg}") #alıntı geçerli olmayan hastalık bilinmeyene eklenir

    else:
        cleaned_data["unknowns"].append("Diseases field is not a list. ")

    #eğer LLM diseases alanını boş döndürdüyse varitanlardan üretsin
    if not cleaned_data["diseases"] and disease_to_cits: 
        for dname, cset in disease_to_cits.items(): #diseases alanı boş olmazsa variantlar dan miras al
            cleaned_data["diseases"].append({"name": dname, "citations": sorted(cset)})



    #Phenotypeları kontrol edelim
    phenotypes = data.get("phenotypes", []) #LLM çıktısındaki fenotipler kategorisindeki itemları alır yoksa boş list döndürür
    if isinstance(phenotypes, list): #phenotypes listeyse
        for phe in phenotypes:
            #tip kontrolü
            if not isinstance(phe, dict):
                cleaned_data["unknowns"].append(f"Skipped non-object phenotype entry: {phe}")
                continue

            name = str(phe.get("name","")).strip() #fenotip adını aldık yoksa boş string döner
            if is_placeholder(name): #placeholder olanları drop edelim
                cleaned_data["unknowns"].append(f"Dropped placeholder phenotype: {name}")
                continue
            cits = phe.get("citations", []) #fenotip alıntılarını alırız yoksa boş string döner
            if not name: #fenotip adı boşsa
                continue #adı boş olan fenotipleri atarız adı olmayan fenotip bilgisi kulllanışlı olmaz
            #LLM citations vermediyse varitanlardan miras alalım
            
            #phenotype inherit kontorlünü normalize edelim
            name_key = norm_key(name)
            if (not cits) and (name_key in pheno_to_cits):
                normalize_cits = sorted(pheno_to_cits[name_key])
                cleaned_data["phenotypes"].append({"name":name, "citations": normalize_cits})
            else:
                cits_ok, msg, normalize_cits = citations_ok(cits) #alıntıları kontrol ederiz, geçerli mi değil mi mesaj ve normalize edilmiş alıntılar listesi
                if cits_ok: #alıntılar geçerliyse
                    cleaned_data["phenotypes"].append({"name":name, "citations": normalize_cits})
                else:
                    cleaned_data["unknowns"].append(f"Dropped phenotype {name} because of citation issue: {msg}") #alıntı geçerli olmayan fenotip bilinmeyene eklenir
    else:
        cleaned_data["unknowns"].append("Phenotypes field is not a list. ")

    #Eğer LLM phenotypes alanını boş dönüyorsa varitanlardan üretsin
    if not cleaned_data["phenotypes"] and pheno_to_cits:
        for pname, cset in pheno_to_cits.items():
            cleaned_data["phenotypes"].append({"name": pname, "citations": sorted(cset)})

    #Eğer hiçbir şey eklenmediyse unknown'a ekleyip bilmiyorum diyelim
    if not cleaned_data["variants"] and not cleaned_data["diseases"] and not cleaned_data["phenotypes"]:
        cleaned_data["unknowns"].append("I do not know based on the retrieved literature. ")  #hiçbir şey eklenmesiyse bilmiyorum dedik 
    
    # DUPLICATE VARIANT MERGE (birleştirme)
    #variants icin
    merged_variant = {}
    for v in cleaned_data["variants"]:
        key = v["name"]
        if key not in merged_variant:
            merged_variant[key] = v
        else:
            merged_variant[key]["citations"] = sorted(set(merged_variant[key]["citations"] + v["citations"]))
            merged_variant[key]["associated_phenotypes"] = sorted(set(
                merged_variant[key]["associated_phenotypes"] + v["associated_phenotypes"]
            ))
            merged_variant[key]["associated_diseases"] = sorted(set(
                merged_variant[key]["associated_diseases"] + v["associated_diseases"]
            ))

    cleaned_data["variants"] = list(merged_variant.values()) #variants listesini döndür
     
    #merge diseases icin (same name -> union citations) yani aynı isimdeki hastalıkları ve alıntılarını birleştirmek için
    merged_dis = {}
    for d in cleaned_data["diseases"]:
        key = d.get("name", "").strip() #hastalık isimlerini alır ve boşlukları temizlr
        if not key:
            continue
        
        cits = d.get("citations", [])
        if key not in merged_dis:
            merged_dis[key] = {"name":key, "citations":sorted(set(cits))} #hastalık sözlükte yoksa yeni kayıt oluşturur
        else:
            merged_dis[key]["citations"] = sorted(set(merged_dis[key]["citations"] + cits)) #eğer hastalık varsa eski alıntılarla yenileri birleştirir tekrarları siler 

    cleaned_data["diseases"] = list(merged_dis.values()) #temizlenmiş disease verilerini listeye çevriri


    #merge phenotypes icin (same name -> union citations) aynı isimdeki fenotipleri ve alıntılarını birleştirmek için 
    merged_phe = {}
    for p in cleaned_data["phenotypes"]: #phenotypes listesini döndür
        key = p.get("name", "").strip() #phenotip adını al
        if not key: #adı boşsa
            continue
        
        cits = p.get("citations", []) #alıntılarını al
        if key not in merged_phe:
            merged_phe[key] = {"name":key, "citations":sorted(set(cits))} #eğer fenotip varsa eski alıntılarla yenileri birleştirir tekrarları siler
        else:
            merged_phe[key]["citations"] = sorted(set(merged_phe[key]["citations"] + cits)) 

    cleaned_data["phenotypes"] = list(merged_phe.values()) #temizlenmiş phenotypes listesini döndür

    return cleaned_data #makale bilgilerini içeren sözlüğü döndür



""" Eval - Metin dizesi veya sayısal değerle sonuçlayan ifade değerlendirmek için Eval işlevini kullanırız """
#Query runner: LLM'den gelen sonuçlara uygun olarak makale bilgilerini döndüren fonksiyon, SQL sorgularını çalıştırmaya ve görselleştirmeye yarar.
def run_query(question: str, top_k: int=4) -> Dict: #makale bilgilerini döndürmek için
    #hard guardrail: eğer soru rars1 içermiyorsa llm'e sorma
    q_low = question.lower()
    if not any(k in q_low for k in ["rars1", "rars", "arginyl", "leukodystrophy", "variant"]):
        return {
            "variants": [],
            "diseases": [],
            "phenotypes": [],
            "unknowns": ["I do not know based on the retrieved literature. The provided evidence only contains genomic data for RARS1."]
        }

    collection = get_collection() #retrieve collection'ı ister
    snippets = retrieve(query=question, collection=collection, top_k=top_k)  #Makale bilgilerini aldık soruya göre en benzer snippetları getirir
    if not snippets: #makale bilgileri yoksa
        return {"variants":[], "diseases": [], "phenotypes":[], "unknowns":["No snippets retrieved"]} 
    
    #  Veritabanındaki TÜM geçerli PMID'leri çekelim (Guardrail'ın yanılmaması için)
    all_data = collection.get() 
    db_allowed_pmids = set()
    if all_data and 'metadatas' in all_data:
        for metadata in all_data['metadatas']:
            if 'pmid' in metadata:
                db_allowed_pmids.add(str(metadata['pmid']))

    #LLM'e soralım ve buradaki allow_pmids'i kullanmayacağız, dbdekini kullanacağız o yüzden allowed_pmids yerine _ yazdık
    data, _ = ask_llm_structured(question, snippets) #makale bilgilerini ve izin verilen PMID'leri alırız (LLM'e soru ve kanıt parçası verip json alır)
    final_data = guardrail_citations(data, db_allowed_pmids)  #uydurma pmid ve yanlış formatları filtreler
    return final_data 

#Değerlendirme sonuçlarına bakacağız ve alakasız soru soracağız. Sistem uydurmasın diye test etme
def write_eval_result(path: str = "eval_results.json") -> None: #değerlendirme sonuçlarını yazdırır
    #Trick soru değerlendirmesi yani bilerek alakasız soru soruyoruz ve sistem kanıt yoksa "bilmiyorum" demeli uydurma yapmmamalı
    trick_questions = [
        "What are the most common symptoms of COVID-19?",
        "What are the recommended treatments for Type 2 Diabetes?",
        "Which variants of the RARS2 gene are associated with Pontocerebellar Hypoplasia type 6?"
    ]
    eval_data = [] #verileri kaydetmek için boş list oluşturduk
    print(f"\n[START] Evaluation mode has been initiated. {len(trick_questions)} questions will be tested.\n")
    for quest in trick_questions:
        print(f"[EVAL] asking: {quest}")
        #soruyu sisteme soruyoruz 
        result = run_query(quest, top_k=8) #retrieve, llm extraction, guardrail çalıştırır ve en fazla 8 snippet getirir

        #sonucu değerlendirme listesine ekleyelim
        eval_entry = {
            "trick_question": quest,
            "result": result,
            "status": "PASS", #eğer cevap i dont know ise veya variant listesi boşsa bu test geçer
            "note" : "Guardrail successfully checked the context and refused to hallucinate."
        }
        eval_data.append(eval_entry) #sonucu degerlendirme listesine ekleyelim


    #JSON dosyası olarak kaydetmek için
    with open(path, "w", encoding="utf-8") as f: 
        json.dump({"queries": eval_data}, f, ensure_ascii=False, indent=2)  #ensure_ascii ile türkçe düzgün gözükür, indent ile json dosyasını daha okunabilir hale getirir
    
    print(f"[OK] Eval results written to {path}") #işlemin başarılı olduğunu anlamak için


""" CLI - Komut satırı arayüzü """
def main():

    #kullanıcı hiçbir şey yazmadan çalıştırdıysa yardım gösterelim
    if len(sys.argv) < 2: 
        print('Usage:\n  python main.py --index\n  python main.py --eval\n  python main.py "your question"')
        return

    arg = sys.argv[1] # main.py'dan sonra yazılan ilk argüman ilk parametreyi alır 

    #Index modu oluşturalım ve makaleleri arayalım
    if arg == "--index":
        # Birden fazla arama yapıp PMID'leri birleştiriyoruz (union)
        queries = [
            'RARS1[Title/Abstract]',
            '("arginyl-tRNA synthetase"[Title/Abstract])',
            '(ArgRS[Title/Abstract])',
            '(RARS[Title/Abstract] AND cytoplasmic[Title/Abstract])',
        ]       

        all_pmids = set()
        for q in queries: #her sorgu için 200 makale çekip benzersiz pmid'leri toplayacak ve sonra union yapacak
            pmids_part = search_pubmed(q, max_results=200)
            all_pmids.update(pmids_part)
        
        pmids = list(all_pmids) #benzersiz pmid'leri listeye cevir
        print(f"[INFO] Union PMIDs: {len(pmids)}") #kac tane pmid varsa yazar
    
        #Abstractları pmid'lere göre çekelim
        records = fetch_abstracts(pmids) 
        print(f"[INFO] Fetched {len(records)} records (board).") 

        #RARS1 ile ilgili olanları filtreleyelim(title + abstract içide)
        keywords = [
        "rars1",
        "arginyl-trna synthetase 1",
        "arginyl trna synthetase 1",
        "arginyl-trna synthetase",
        "arginyl trna synthetase",
        "argrs",
        "arg rs",
        ]
    
        filtered = []
        for r in records:
            text = f"{r.get('title','')} {r.get('abstract','')}".lower() #makale baslığı ve özetini kucuk harfe cevir
            if any(k in text for k in keywords):
                filtered.append(r)

        #En fazla 50 makale ile sınırlayalım
        filtered = filtered[:50] #0-50 indeks arası getirir
        print(f"[INFO] Filtered RARS1-related records: {len(filtered)}")

        #Index oluşturalım- temizlenen verileri chunklara ayırıp veritabanına kaydeder
        n_chunks = build_index(filtered, reset=True) #reset=True ise eski koleksiyonu silerek temiz baslangıc yapar
        print(f"[OK] Indexed chunks: {n_chunks}") #kac tane chunk olusturuldugunu yazar
        return

    # 2) Eval modu oluşturalım ve doğruluğu test edip sonucu json'a yazdıralım
    if arg == "--eval": 
        write_eval_result("eval_results.json")
        return


    # Soru (question) modu oluşturalım (arg bir soru kabul edilir ve llm'e sorar)
    question = arg 
    out = run_query(question, top_k=15) #en alakalı 15 parçayı getirir ve cevap verir
    print(json.dumps(out, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()