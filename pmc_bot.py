import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from Bio import Entrez
from datetime import datetime

# ================= 사용자 설정 =================
# 받는 이메일
TARGET_EMAIL = "cybog337@gmail.com"

# 검색어: "biogems" AND "last 14 days"[pdat]
SEARCH_TERM = '"biogems" AND "last 14 days"[pdat]'

# NCBI 접속용 이메일
Entrez.email = TARGET_EMAIL

# GitHub Secrets에서 비밀번호 가져오기
GMAIL_PASSWORD = os.environ.get("GMAIL_PASSWORD")
# =============================================

def fetch_pmc_new_articles(term):
    try:
        # 1. 검색 실행 (최신순 정렬)
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
            
            # 정보 추출 (없을 경우 빈칸 처리)
            title = doc.get('Title', 'No Title')
            pub_date = doc.get('PubDate', 'Date Unknown') # 보통 'YYYY Mon DD' 형식
            source = doc.get('Source', '')   # 저널명
            volume = doc.get('Volume', '')   # 권
            issue = doc.get('Issue', '')     # 호
            pages = doc.get('Pages', '')     # 페이지
            
            # Citation 구성 (예: Nature. 2024;10(2):100-110.)
            citation = f"{source}"
            if volume:
                citation += f". {volume}"
            if issue:
                citation += f"({issue})"
            if pages:
                citation += f":{pages}."
            
            link = f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/"
            
            # 요청하신 포맷 적용
            # 1. [ 날짜 ]
            # 2. (줄바꿈) 논문 제목
            # 3. (줄바꿈) Citation
            # 4. (줄바꿈) 링크
            entry = f"[ {pub_date} ]\n{title}\nCitation: {citation}\nLink: {link}"
            
            articles.append(entry)
            
        return articles

    except Exception as e:
        print(f"Error fetching data: {e}")
        return []

def send_email(articles):
    if not articles:
        print("검색 결과가 없습니다.")
        return

    msg = MIMEMultipart()
    msg['From'] = TARGET_EMAIL
    msg['To'] = TARGET_EMAIL
    msg['Subject'] = f"[PMC 알림] 'biogems' 관련 새 논문 {len(articles)}건"

    # 메일 본문 구성 (가독성을 위해 항목 간 줄바꿈 추가)
    body = f"검색어: {SEARCH_TERM}\n\n" + ("-" * 30) + "\n\n"
    body += "\n\n".join(articles)
    
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
