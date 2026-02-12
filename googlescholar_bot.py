import os
import re
import smtplib
import requests
from bs4 import BeautifulSoup
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime
import time

# ================= 사용자 설정 =================
TARGET_EMAIL = "cybog337@gmail.com"
SEARCH_QUERY = "biogems -biogem -cjter"
HISTORY_FILE = "sent_list_scholar.txt"

GMAIL_PASSWORD = os.environ.get("GMAIL_PASSWORD")
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

def extract_date_info(text):
    """날짜 정보 추출"""
    match = re.search(r'(202[0-9])\s*-?\s*(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)?', text)
    if match:
        year = match.group(1)
        month = match.group(2) if match.group(2) else ""
        return f"{year} {month}".strip()
    return "2026"

def fetch_scholar_data_html():
    """Google Scholar HTML 직접 파싱"""
    all_articles = []
    start_index = 0
    page_num = 1
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }
    
    while True:
        url = f"https://scholar.google.com/scholar?q={SEARCH_QUERY}&hl=ko&as_sdt=0,5&as_ylo=2026&filter=0&start={start_index}"
        
        try:
            response = requests.get(url, headers=headers, timeout=10)
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'html.parser')
            
            results = soup.find_all('div', class_='gs_ri')
            
            print(f"Page {page_num}: {len(results)}건 수집 (start={start_index})")
            
            if not results:
                break
            
            for result in results:
                # 제목 추출
                title_tag = result.find('h3', class_='gs_rt')
                if title_tag:
                    # 링크 추출
                    link_tag = title_tag.find('a')
                    title = title_tag.get_text(strip=True)
                    link = link_tag['href'] if link_tag and link_tag.has_attr('href') else "No Link"
                else:
                    title = "No Title"
                    link = "No Link"
                
                # 저자/저널 정보 추출
                info_tag = result.find('div', class_='gs_a')
                info = info_tag.get_text(strip=True) if info_tag else "정보 없음"
                
                # 날짜 추출
                date_str = extract_date_info(info)
                
                all_articles.append({
                    "title": title,
                    "link": link,
                    "info": info,
                    "date": date_str
                })
            
            # 다음 페이지 확인
            next_button = soup.find('button', {'aria-label': '다음 페이지'}) or soup.find('a', string='다음')
            if next_button and not next_button.get('disabled'):
                start_index += 10
                page_num += 1
                time.sleep(2)  # Google 차단 방지
            else:
                break
                
        except Exception as e:
            print(f"Error fetching page {page_num}: {e}")
            break
    
    print(f"총 수집 건수: {len(all_articles)}건")
    return all_articles

def filter_new_articles(all_articles, sent_urls):
    """이미 발송한 논문 제외하고 신규만 필터링"""
    new_articles = []
    for article in all_articles:
        if article["link"] not in sent_urls:
            new_articles.append(article)
    return new_articles

def send_report(articles):
    """메일 발송"""
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
    if not GMAIL_PASSWORD:
        print("환경변수 설정 필요: GMAIL_PASSWORD")
        exit(1)
    
    sent_history = load_sent_history()
    print(f"기존 이력: {len(sent_history)}건")
    
    all_articles = fetch_scholar_data_html()
    print(f"검색 결과: {len(all_articles)}건")
    
    new_articles = filter_new_articles(all_articles, sent_history)
    print(f"신규 논문: {len(new_articles)}건")
    
    if send_report(new_articles):
        new_urls = [article["link"] for article in new_articles]
        save_sent_history(new_urls)
        print("이력 업데이트 완료")
