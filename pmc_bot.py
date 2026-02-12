import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from Bio import Entrez
from datetime import datetime

# ================= 사용자 설정 =================
TARGET_EMAIL = "cybog337@gmail.com"
SEARCH_TERM = '"biogems" AND "last 60 days"[pdat]'

Entrez.email = TARGET_EMAIL
GMAIL_PASSWORD = os.environ.get("GMAIL_PASSWORD")
# =============================================

def fetch_pmc_new_articles(term):
    try:
        # 1. 검색 실행
        handle = Entrez.esearch(db="pmc", term=term, retmax=20, sort='pub_date')
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        count = int(record["Count"])
        
        if count == 0:
            return []

        # 2. 상세 정보 가져오기
        handle = Entrez.esummary(db="pmc", id=",".join(id_list))
        summary_record = Entrez.read(handle)
        handle.close()
        
        articles = []
        for doc in summary_record:
            pmc_id = doc['Id']
            title = doc.get('Title', 'No Title')
            
            # [날짜 로직 개선]
            # PubDate(출판일)는 월까지만 있는 경우가 많으므로,
            # EPubDate(전자출판일)가 있으면 그걸 우선적으로 사용합니다.
            pub_date = doc.get('PubDate', '')
            epub_date = doc.get('EPubDate', '')
            
            if epub_date and len(epub_date) > len(pub_date):
                display_date = epub_date
            else:
                display_date = pub_date
            
            # Citation 구성
            source = doc.get('Source', '')
            volume = doc.get('Volume', '')
            issue = doc.get('Issue', '')
            pages = doc.get('Pages', '')
            
            citation = f"{source}"
            if volume: citation += f". {volume}"
            if issue:  citation += f"({issue})"
            if pages:  citation += f":{pages}."
            
            link = f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/"
            
            # [출력 형식 수정]
            # Citation:, Link: 글자 제거
            entry = f"[ {display_date} ]\n{title}\n{citation}\n{link}"
            
            articles.append(entry)
            
        return articles

    except Exception as e:
        print(f"Error fetching data: {e}")
        return []

def send_email(articles):
    msg = MIMEMultipart()
    msg['From'] = TARGET_EMAIL
    msg['To'] = TARGET_EMAIL
    
    # [수정] 결과가 있을 때와 없을 때 제목을 다르게 표시
    if articles:
        msg['Subject'] = f"[PMC 알림] 'biogems' 관련 새 논문 {len(articles)}건"
        body = f"검색어: {SEARCH_TERM}\n\n" + ("-" * 30) + "\n\n"
        body += "\n\n".join(articles)
    else:
        msg['Subject'] = f"[PMC 알림] 'biogems' 관련 새 논문 0건"
        body = f"검색어: {SEARCH_TERM}\n\n최근 1일간 검색된 새로운 논문이 없습니다.\n시스템은 정상 작동 중입니다."

    msg.attach(MIMEText(body, 'plain'))

    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.starttls()
        server.login(TARGET_EMAIL, GMAIL_PASSWORD)
        server.send_message(msg)
        server.quit()
        print(f"✅ 이메일 발송 완료! ({len(articles)}건)")
    except Exception as e:
        print(f"❌ 이메일 발송 실패: {e}")

if __name__ == "__main__":
    if not GMAIL_PASSWORD:
        print("Error: GMAIL_PASSWORD 환경변수가 설정되지 않았습니다.")
    else:
        print(f"PMC Searching for: {SEARCH_TERM}")
        new_articles = fetch_pmc_new_articles(SEARCH_TERM)
        send_email(new_articles)
