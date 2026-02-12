import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from Bio import Entrez
from datetime import datetime

# ================= 사용자 설정 =================
# 요청하신 이메일과 검색어 설정
TARGET_EMAIL = "cybog337@gmail.com"
SEARCH_TERM = '"biogems" AND "last 14 days"[pdat]'

# NCBI에 등록될 이메일 (필수)
Entrez.email = TARGET_EMAIL

# GitHub Secrets에서 가져올 앱 비밀번호 (16자리)
GMAIL_PASSWORD = os.environ.get("GMAIL_PASSWORD")
# =============================================

def fetch_pmc_new_articles(term):
    try:
        # 검색 실행
        handle = Entrez.esearch(db="pmc", term=term, retmax=20)
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        if int(record["Count"]) == 0:
            return []

        # 세부 정보 가져오기
        handle = Entrez.esummary(db="pmc", id=",".join(id_list))
        summary_record = Entrez.read(handle)
        handle.close()
        
        articles = []
        for doc in summary_record:
            pmc_id = doc['ArticleIds']['pmc']
            title = doc.get('Title', 'No Title')
            pub_date = doc.get('PubDate', 'Unknown Date')
            link = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/"
            articles.append(f"[{pub_date}] {title}\nLink: {link}")
            
        return articles
    except Exception as e:
        print(f"Error: {e}")
        return []

def send_email(articles):
    if not articles:
        print("No new articles found.")
        return

    msg = MIMEMultipart()
    msg['From'] = TARGET_EMAIL
    msg['To'] = TARGET_EMAIL  # 보내는 사람 = 받는 사람
    msg['Subject'] = f"[BioGems 알림] 새로운 논문 {len(articles)}건 발견"

    body = f"검색어: {SEARCH_TERM}\n\n" + "\n\n".join(articles)
    msg.attach(MIMEText(body, 'plain'))

    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.starttls()
        server.login(TARGET_EMAIL, GMAIL_PASSWORD)
        server.send_message(msg)
        server.quit()
        print(f"Email sent successfully to {TARGET_EMAIL}")
    except Exception as e:
        print(f"Failed to send email: {e}")

if __name__ == "__main__":
    if not GMAIL_PASSWORD:
        print("Error: GMAIL_PASSWORD 환경변수가 설정되지 않았습니다.")
    else:
        print(f"Searching for: {SEARCH_TERM}")
        new_articles = fetch_pmc_new_articles(SEARCH_TERM)
        send_email(new_articles)
