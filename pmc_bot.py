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
            
            # [날짜 로직 수정: 무결성 보장]
            # PubDate가 불완전할 경우를 대비해 History 데이터에서 실제 등록일(YYYY/MM/DD)을 추적합니다.
            display_date = doc.get('EPubDate', '') or doc.get('PubDate', '')
            history = doc.get('History', [])
            for h in history:
                # epublish(전자출판) 또는 pubmed(등록일) 상태의 상세 날짜를 우선 채택
                if h.get('PubStatus') in ['epublish', 'pubmed']:
                    raw_date = h.get('Date', '').split(' ')[0] # '2026/02/12' 형태 추출
                    try:
                        dt = datetime.strptime(raw_date, '%Y/%m/%d')
                        display_date = dt.strftime('%Y %b %d') # '2026 Feb 12'로 변환
                        break
                    except:
                        continue
            
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
            
            # [출력 형식] 캐리지 리턴 반영 및 라벨 제거
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
    
    # [수정] 스콜라와 동일한 날짜 표기 방식 적용
    date_str = datetime.now().strftime("%Y-%m-%d")
    count = len(articles)
    
    if articles:
        # 규격: [PMC] 2026-02-12 신규 논문 알림 (N건)
        msg['Subject'] = f"[PMC] {date_str} 신규 논문 알림 ({count}건)"
        body = f"검색어: {SEARCH_TERM}\n\n" + ("-" * 30) + "\n\n"
        body += "\n\n".join(articles)
    else:
        # 0건일 때도 날짜를 포함하여 시스템 생존 확인
        msg['Subject'] = f"[PMC] {date_str} 신규 논문 알림 (0건)"
        body = f"검색어: {SEARCH_TERM}\n\n최근 검색된 새로운 논문이 없습니다.\n시스템은 정상 작동 중입니다."

    msg.attach(MIMEText(body, 'plain'))

    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.starttls()
        server.login(TARGET_EMAIL, GMAIL_PASSWORD)
        server.send_message(msg)
        server.quit()
        print(f"✅ [{date_str}] PMC 메일 발송 완료! ({count}건)")
    except Exception as e:
        print(f"❌ 이메일 발송 실패: {e}")

if __name__ == "__main__":
    if not GMAIL_PASSWORD:
        print("Error: GMAIL_PASSWORD 환경변수가 설정되지 않았습니다.")
    else:
        print(f"PMC Searching for: {SEARCH_TERM}")
        new_articles = fetch_pmc_new_articles(SEARCH_TERM)
        send_email(new_articles)
