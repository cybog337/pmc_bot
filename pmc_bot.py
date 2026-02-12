import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from Bio import Entrez
from datetime import datetime

# ================= 사용자 설정 =================
# 받는 이메일
TARGET_EMAIL = "cybog337@gmail.com"

# 검색어: 보내주신 URL의 조건 그대로 적용 ("biogems" AND "last 14 days"[pdat])
SEARCH_TERM = '"biogems" AND "last 14 days"[pdat]'

# NCBI 접속용 이메일 (필수)
Entrez.email = TARGET_EMAIL

# GitHub Secrets에서 비밀번호 가져오기
GMAIL_PASSWORD = os.environ.get("GMAIL_PASSWORD")
# =============================================

def fetch_pmc_new_articles(term):
    try:
        # 1. PMC 데이터베이스에서 검색 (db="pmc")
        # sort='pub_date': 최신순 정렬 (URL의 &sort=pubdate 반영)
        handle = Entrez.esearch(db="pmc", term=term, retmax=20, sort='pub_date')
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        count = int(record["Count"])
        
        if count == 0:
            return []

        # 2. 요약 정보 가져오기
        handle = Entrez.esummary(db="pmc", id=",".join(id_list))
        summary_record = Entrez.read(handle)
        handle.close()
        
        articles = []
        for doc in summary_record:
            # PMC 데이터에서는 'Id'가 곧 PMCID 번호입니다.
            pmc_id = doc['Id'] 
            title = doc.get('Title', '제목 없음')
            pub_date = doc.get('PubDate', '날짜 미상')
            
            # 보내주신 링크 형식인 pmc.ncbi.nlm.nih.gov 도메인 적용
            link = f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/"
            
            articles.append(f"[{pub_date}] {title}\nLink: {link}")
            
        return articles

    except Exception as e:
        print(f"Error fetching data: {e}")
        return []

def send_email(articles):
    if not articles:
        print("검색 결과가 0건이라 메일을 보내지 않습니다.")
        return

    msg = MIMEMultipart()
    msg['From'] = TARGET_EMAIL
    msg['To'] = TARGET_EMAIL
    msg['Subject'] = f"[PMC 알림] 'biogems' 관련 새 논문 {len(articles)}건"

    body = f"설정하신 검색어: {SEARCH_TERM}\n\n" + "\n\n".join(articles)
    msg.attach(MIMEText(body, 'plain'))

    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.starttls()
        server.login(TARGET_EMAIL, GMAIL_PASSWORD)
        server.send_message(msg)
        server.quit()
        print(f"✅ 이메일 발송 성공! ({len(articles)}건)")
    except Exception as e:
        print(f"❌ 이메일 발송 실패: {e}")

if __name__ == "__main__":
    if not GMAIL_PASSWORD:
        print("Error: GMAIL_PASSWORD 환경변수가 설정되지 않았습니다.")
    else:
        print(f"PMC Searching for: {SEARCH_TERM}")
        new_articles = fetch_pmc_new_articles(SEARCH_TERM)
        send_email(new_articles)
