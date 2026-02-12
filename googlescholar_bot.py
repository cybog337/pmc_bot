import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from serpapi import GoogleSearch # 구글 검색용 라이브러리
from datetime import datetime

# ================= 사용자 설정 =================
TARGET_EMAIL = "cybog337@gmail.com"

# 검색어: biogems 포함, biogem 제외, cjter 제외
SEARCH_QUERY = "biogems -biogem -cjter"

# GitHub Secrets에서 키 가져오기
GMAIL_PASSWORD = os.environ.get("GMAIL_PASSWORD")
SERPAPI_KEY = os.environ.get("SERPAPI_KEY") # 구글 스콜라 전용 키
# =============================================

def fetch_scholar_new_articles():
    if not SERPAPI_KEY:
        print("❌ Error: SERPAPI_KEY가 없습니다. GitHub Settings > Secrets에 키를 등록해주세요.")
        return []

    # SerpApi 설정 (구글 스콜라 엔진 사용)
    params = {
        "engine": "google_scholar",
        "q": SEARCH_QUERY,
        "api_key": SERPAPI_KEY,
        "scisbd": "1",  # 1 = 최신순 정렬 (Sort by date) - 가장 중요!
        "num": "10",    # 최근 10개만 가져오기
        "hl": "en"      # 언어 설정 (영어)
    }

    try:
        print(f"Searching Google Scholar for: {SEARCH_QUERY} ...")
        search = GoogleSearch(params)
        results = search.get_dict()
        
        # 검색 결과가 없을 경우 처리
        if "organic_results" not in results:
            return []

        organic_results = results["organic_results"]
        articles = []

        for result in organic_results:
            title = result.get("title", "No Title")
            link = result.get("link", "#")
            snippet = result.get("snippet", "No Snippet")
            
            # 출판 정보 (저자, 저널, 연도 등)
            pub_info = result.get("publication_info", {}).get("summary", "")

            # 메일 본문 포맷
            entry = f"[New] {title}\nInfo: {pub_info}\nSnippet: {snippet}\nLink: {link}"
            articles.append(entry)

        return articles

    except Exception as e:
        print(f"Error fetching Scholar data: {e}")
        return []

def send_email(articles):
    if not articles:
        print("검색 결과가 0건이라 메일을 보내지 않습니다.")
        return

    msg = MIMEMultipart()
    msg['From'] = TARGET_EMAIL
    msg['To'] = TARGET_EMAIL
    msg['Subject'] = f"[Scholar 알림] 'biogems' 관련 최신 논문 {len(articles)}건"

    body = f"검색어: {SEARCH_QUERY}\n(옵션: 최신순 정렬)\n\n" + ("-" * 30) + "\n\n"
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
        print("Error: GMAIL_PASSWORD가 설정되지 않았습니다.")
    elif not SERPAPI_KEY:
        print("Error: SERPAPI_KEY가 설정되지 않았습니다.")
    else:
        new_articles = fetch_scholar_new_articles()
        send_email(new_articles)
