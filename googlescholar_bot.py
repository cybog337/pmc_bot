import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from serpapi import GoogleSearch
from datetime import datetime

# ================= 사용자 설정 =================
TARGET_EMAIL = "cybog337@gmail.com"
SEARCH_QUERY = "biogems -biogem -cjter"

GMAIL_PASSWORD = os.environ.get("GMAIL_PASSWORD") # ${{ secrets.MY_PASSWORD }}
SERPAPI_KEY = os.environ.get("SERPAPI_KEY")
# =============================================

def fetch_scholar_full_scan():
    if not SERPAPI_KEY: return []
    all_articles = []
    start_index = 0
    
    while start_index < 50: # 최대 5페이지까지 전수 조사
        params = {
            "engine": "google_scholar",
            "q": SEARCH_QUERY,
            "api_key": SERPAPI_KEY,
            "as_ylo": "2026",
            "start": start_index,
            "hl": "ko"
        }
        try:
            search = GoogleSearch(params)
            results = search.get_dict()
            organic_results = results.get("organic_results", [])
            if not organic_results: break

            for result in organic_results:
                all_articles.append({
                    "title": result.get("title", "No Title"),
                    "link": result.get("link", "#"),
                    "info": result.get("publication_info", {}).get("summary", "No Info")
                })

            if "next" in results.get("serpapi_pagination", {}):
                start_index += 10
            else: break
        except: break
    return all_articles

def send_final_report(articles):
    msg = MIMEMultipart()
    msg['From'] = TARGET_EMAIL
    msg['To'] = TARGET_EMAIL
    
    date_str = datetime.now().strftime("%Y-%m-%d")
    count = len(articles)
    
    # 제목 양식 통일
    msg['Subject'] = f"[Scholar] {date_str} 신규 논문 알림 ({count}건)"
    
    if articles:
        # PMC 예시 양식과 동일하게 구성
        body_parts = []
        for item in articles:
            part = f"[ 2026 Feb ]\n{item['title']}\n{item['info']}\n{item['link']}"
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
        found = fetch_scholar_full_scan()
        send_final_report(found)
