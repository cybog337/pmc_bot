import os
import re
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from serpapi import GoogleSearch
from datetime import datetime

# ================= 사용자 설정 =================
TARGET_EMAIL = "cybog337@gmail.com"
SEARCH_QUERY = "biogems -biogem -cjter"
HISTORY_FILE = "sent_list.txt"

GMAIL_PASSWORD = os.environ.get("GMAIL_PASSWORD") 
SERPAPI_KEY = os.environ.get("SERPAPI_KEY")
# =============================================

def load_sent_history():
    """이전에 발송한 URL 목록 로드"""
    if not os.path.exists(HISTORY_FILE):
        return set()
    
    with open(HISTORY_FILE, 'r', encoding='utf-8') as f:
        return set(line.strip() for line in f if line.strip())

def save_sent_history(urls):
    """새로 발송한 URL을 이력 파일에 추가"""
    with open(HISTORY_FILE, 'a', encoding='utf-8') as f:
        for url in urls:
            f.write(url + '\n')

def extract_date_info(pub_info, snippet=""):
    """
    publication_info와 snippet에서 동적으로 날짜 추출
    예: "2026 Jan", "2026 Feb" 등
    """
    combined_text = pub_info + " " + snippet
    
    match = re.search(r'(202[0-9])\s+(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)', combined_text)
    if match:
        return f"{match.group(1)} {match.group(2)}"
    
    match_year = re.search(r'202[0-9]', combined_text)
    if match_year:
        return match_year.group(0)
    
    return "2026"

def fetch_scholar_data():
    """Google Scholar에서 전체 데이터 수집"""
    if not SERPAPI_KEY: 
        return []
    
    all_articles = []
    start_index = 0
    
    while True:
        params = {
            "engine": "google_scholar",
            "q": SEARCH_QUERY,
            "api_key": SERPAPI_KEY,
            "as_ylo": "2026",
            "as_sdt": "0,5",
            "filter": "0",
            "start": start_index,
            "hl": "ko"
        }
        
        try:
            search = GoogleSearch(params)
            results = search.get_dict()
            organic_results = results.get("organic_results", [])
            
            if not organic_results: 
                break

            for result in organic_results:
                pub_info = result.get("publication_info", {}).get("summary", "")
                snippet = result.get("snippet", "")
                
                date_str = extract_date_info(pub_info, snippet)
                
                all_articles.append({
                    "title": result.get("title", "No Title"),
                    "link": result.get("link", "No Link"),
                    "info": pub_info if pub_info else "정보 없음",
                    "date": date_str
                })

            if "next" in results.get("serpapi_pagination", {}):
                start_index += 10
            else: 
                break
                
        except Exception as e:
            print(f"Error fetching data: {e}")
            break
    
    return all_articles

def filter_new_articles(all_articles, sent_urls):
    """이미 발송한 논문 제외하고 신규만 필터링"""
    new_articles = []
    for article in all_articles:
        if article["link"] not in sent_urls:
            new_articles.append(article)
    return new_articles

def send_report(articles):
    """메일 발송 (PMC 양식)"""
    msg = MIMEMultipart()
    msg['From'] = TARGET_EMAIL
    msg['To'] = TARGET_EMAIL
    
    date_str = datetime.now().strftime("%Y-%m-%d")
    count = len(articles)
    msg['Subject'] = f"[Scholar] {date_str} 신규 논문 알림 ({count}건)"
    
    if articles:
        body_parts = []
        for item in articles:
            part = f"[ {item['date']} ]\n{item['title']}\n{item['info']}\n{item['link']}"
            body_parts.append(part)
        
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
        return True
    except Exception as e:
        print(f"메일 발송 실패: {e}")
        return False

if __name__ == "__main__":
    if not GMAIL_PASSWORD or not SERPAPI_KEY:
        print("환경변수 설정 필요: GMAIL_PASSWORD, SERPAPI_KEY")
        exit(1)
    
    sent_history = load_sent_history()
    print(f"기존 이력: {len(sent_history)}건")
    
    all_articles = fetch_scholar_data()
    print(f"검색 결과: {len(all_articles)}건")
    
    new_articles = filter_new_articles(all_articles, sent_history)
    print(f"신규 논문: {len(new_articles)}건")
    
    if send_report(new_articles):
        new_urls = [article["link"] for article in new_articles]
        save_sent_history(new_urls)
        print("이력 업데이트 완료")
