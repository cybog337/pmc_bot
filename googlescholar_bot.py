import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from serpapi import GoogleSearch
from datetime import datetime

# ================= 사용자 설정 =================
TARGET_EMAIL = "cybog337@gmail.com"
SEARCH_QUERY = "biogems -biogem -cjter"

GMAIL_PASSWORD = os.environ.get("GMAIL_PASSWORD") 
SERPAPI_KEY = os.environ.get("SERPAPI_KEY")
# =============================================

def fetch_scholar_data():
    if not SERPAPI_KEY: return []
    all_articles = []
    start_index = 0
    
    while True:
        params = {
            "engine": "google_scholar",
            "q": SEARCH_QUERY,
            "api_key": SERPAPI_KEY,
            "as_ylo": "2026",
            "as_sdt": "0,5",  # 31건 맞추기 위한 필수 옵션 (인용/특허 포함)
            "filter": "0",    # 유사 결과 생략 해제
            "start": start_index,
            "hl": "ko"
        }
        try:
            search = GoogleSearch(params)
            results = search.get_dict()
            organic_results = results.get("organic_results", [])
            if not organic_results: break

            for result in organic_results:
                # [핵심] 실제 날짜 정보 추출 (예: "2026 Jan", "2026 Feb" 등)
                pub_info = result.get("publication_info", {}).get("summary", "2026")
                
                all_articles.append({
                    "title": result.get("title", "No Title"),
                    "link": result.get("link", "No Link"),
                    "info": pub_info
                })

            if "next" in results.get("serpapi_pagination", {}):
                start_index += 10
            else: break
        except: break
    return all_articles

def send_report(articles):
    msg = MIMEMultipart()
    msg['From'] = TARGET_EMAIL
    msg['To'] = TARGET_EMAIL
    
    date_str = datetime.now().strftime("%Y-%m-%d")
    count = len(articles)
    msg['Subject'] = f"[Scholar] {date_str} 신규 논문 알림 ({count}건)"
    
    if articles:
        body_parts = []
        for item in articles:
            # [날짜 정보]를 가공하여 상단에 배치 (하드코딩 제거)
            # pub_info에서 연도/월 정보를 최대한 살려 표시합니다.
            display_date = item['info'].split('-')[0].strip() # 날짜 부분만 추출 시도
            part = f"[ {display_date} ]\n{item['title']}\n{item['info']}\n{item['link']}"
            body_parts.append(part)
        body = "\n\n".join(body_parts)
    else:
        body = "신규 논문이 없습니다."
    
    msg.attach(MIMEText(body, 'plain'))

    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.starttls()
        server.login(TARGET_EMAIL, GMAIL_PASSWORD)
        server.send_message(msg)
        server.quit()
    except: pass

if __name__ == "__main__":
    if GMAIL_PASSWORD and SERPAPI_KEY:
        data = fetch_scholar_data()
        send_report(data)
