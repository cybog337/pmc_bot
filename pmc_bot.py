import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from Bio import Entrez
from datetime import datetime

# ================= 사용자 설정 =================
TARGET_EMAIL = "cybog337@gmail.com"
SEARCH_TERM = '"biogems" AND "last 1 day"[pdat]'  # ✅ 60일로 변경
HISTORY_FILE = "sent_list_pubmed.txt"  # ✅ 추가

Entrez.email = TARGET_EMAIL
GMAIL_PASSWORD = os.environ.get("GMAIL_PASSWORD")
# =============================================

def load_sent_history():
    """이전에 발송한 PMC ID 목록 로드"""  # ✅ 추가
    if not os.path.exists(HISTORY_FILE):
        return set()
    
    with open(HISTORY_FILE, 'r', encoding='utf-8') as f:
        return set(line.strip() for line in f if line.strip())

def save_sent_history(pmc_ids):
    """새로 발송한 PMC ID를 이력 파일에 추가"""  # ✅ 추가
    with open(HISTORY_FILE, 'a', encoding='utf-8') as f:
        for pmc_id in pmc_ids:
            f.write(pmc_id + '\n')

def fetch_pmc_new_articles(term):
    try:
        handle = Entrez.esearch(db="pmc", term=term, retmax=50, sort='pub_date')  # ✅ retmax 증가
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        count = int(record["Count"])
        
        print(f"PMC 검색 결과: {count}건")  # ✅ 로그 추가
        
        if count == 0:
            return []

        handle = Entrez.esummary(db="pmc", id=",".join(id_list))
        summary_record = Entrez.read(handle)
        handle.close()
        
        articles = []
        for doc in summary_record:
            pmc_id = doc['Id']
            title = doc.get('Title', 'No Title')
            
            pub_date = doc.get('PubDate', '')
            epub_date = doc.get('EPubDate', '')
            
            if len(epub_date) > len(pub_date):
                display_date = epub_date
            else:
                display_date = pub_date
            
            source = doc.get('Source', '')
            volume = doc.get('Volume', '')
            issue = doc.get('Issue', '')
            pages = doc.get('Pages', '')
            
            citation = f"{source}"
            if volume: citation += f". {volume}"
            if issue:  citation += f"({issue})"
            if pages:  citation += f":{pages}."
            
            link = f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/"
            
            articles.append({
                "pmc_id": pmc_id,  # ✅ ID 추가
                "entry": f"[ {display_date} ]\n{title}\n{citation}\n{link}"
            })
            
        return articles

    except Exception as e:
        print(f"Error fetching data: {e}")
        return []

def filter_new_articles(all_articles, sent_ids):
    """이미 발송한 논문 제외하고 신규만 필터링"""  # ✅ 추가
    new_articles = []
    for article in all_articles:
        if article["pmc_id"] not in sent_ids:
            new_articles.append(article)
    return new_articles

def send_email(articles):
    date_str = datetime.now().strftime("%Y-%m-%d")  # ✅ 추가
    count = len(articles)
    
    msg = MIMEMultipart()
    msg['From'] = TARGET_EMAIL
    msg['To'] = TARGET_EMAIL
    msg['Subject'] = f"[PubMed] {date_str} 신규 논문 알림 ({count}건)"  # ✅ 변경

    if articles:
        body_parts = [article["entry"] for article in articles]
        body = "\n\n".join(body_parts)
    else:
        body = "신규 논문이 없습니다."
    
    msg.attach(MIMEText(body, 'plain', 'utf-8'))

    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.starttls()
        server.login(TARGET_EMAIL, GMAIL_PASSWORD)
        server.send_message(msg)
        server.quit()
        print(f"메일 발송 완료: {count}건")
        return True  # ✅ 추가
    except Exception as e:
        print(f"메일 발송 실패: {e}")
        return False  # ✅ 추가

if __name__ == "__main__":
    if not GMAIL_PASSWORD:
        print("환경변수 설정 필요: GMAIL_PASSWORD")
        exit(1)
    
    # ✅ 증분 보고 로직 추가
    sent_history = load_sent_history()
    print(f"기존 이력: {len(sent_history)}건")
    
    all_articles = fetch_pmc_new_articles(SEARCH_TERM)
    print(f"검색 결과: {len(all_articles)}건")
    
    new_articles = filter_new_articles(all_articles, sent_history)
    print(f"신규 논문: {len(new_articles)}건")
    
    if send_email(new_articles):
        new_ids = [article["pmc_id"] for article in new_articles]
        save_sent_history(new_ids)
        print("이력 업데이트 완료")
