import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from Bio import Entrez
from datetime import datetime, timezone, timedelta

# ================= 사용자 설정 =================
SENDER_EMAIL = os.environ.get("MY_EMAIL")
APP_PASSWORD = os.environ.get("MY_PASSWORD")
RECEIVER_EMAIL = os.environ.get("TO_EMAIL")

# 검색어 설정
SEARCH_QUERY = "biogems"

Entrez.email = SENDER_EMAIL 

def fetch_pmc_articles():
    """최근 24시간 논문 검색"""
    print(f"🔍 검색어 [{SEARCH_QUERY}] 검색 시작...")
    
    try:
        # 최근 1일(reldate=1), 등록일 기준(datetype="edat")
        handle = Entrez.esearch(db="pmc", term=SEARCH_QUERY, reldate=1, datetype="edat", retmax=20)
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        total_count = int(record["Count"])
        
        if not id_list:
            return [], 0
            
        handle = Entrez.esummary(db="pmc", id=",".join(id_list))
        summaries = Entrez.read(handle)
        handle.close()
        
        extracted_data = []
        for item in summaries:
            title = item.get("Title", "제목 없음")
            journal = item.get("Source", "저널명 미상")
            pub_date = item.get("PubDate", "")
            authors = item.get("AuthorList", [])
            author_str = f"{authors[0]} et al." if len(authors) > 1 else (authors[0] if authors else "저자 미상")
            
            pmc_id = item.get("ArticleIds", {}).get("pmcid", "")
            link = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/" if pmc_id else "#"
            
            extracted_data.append({
                "title": title,
                "citation": f"{author_str}, <i>{journal}</i> ({pub_date})",
                "link": link
            })
            
        return extracted_data, total_count

    except Exception as e:
        print(f"❌ 에러 발생: {e}")
        return None, 0

def send_email(articles, count):
    """메일 발송 (UTF-8 인코딩 수정됨)"""
    msg = MIMEMultipart()
    msg['From'] = SENDER_EMAIL
    msg['To'] = RECEIVER_EMAIL
    
    kst_now = datetime.now(timezone.utc) + timedelta(hours=9)
    date_str = kst_now.strftime("%Y-%m-%d")

    if articles:
        msg['Subject'] = f"[PMC 알림] '{SEARCH_QUERY}' 신규 논문 {count}건 ({date_str})"
        html_body = f"""
        <h2 style="color: #2c3e50;">📅 {date_str} 검색 결과</h2>
        <p>검색어 <b>'{SEARCH_QUERY}'</b> 관련 {count}건의 논문이 발견되었습니다.</p>
        <hr>
        """
        for art in articles:
            html_body += f"""
            <div style="padding: 10px; border-bottom: 1px solid #ddd;">
                <b>{art['title']}</b><br>
                <span style="color: #666;">{art['citation']}</span><br>
                <a href="{art['link']}" style="color: #007bff;">🔗 원문 보기</a>
            </div>
            """
    else:
        msg['Subject'] = f"[PMC 알림] '{SEARCH_QUERY}' 신규 논문 없음 ({date_str})"
        html_body = f"""
        <h3>📅 {date_str} 검색 결과</h3>
        <p>'{SEARCH_QUERY}' 관련 새로운 논문이 없습니다.</p>
        <p style="color: gray; font-size: 12px;">내일 아침 7시에 다시 확인합니다.</p>
        """

    # [중요 수정] 'utf-8'을 명시하여 한글 깨짐 방지
    msg.attach(MIMEText(html_body, 'html', 'utf-8'))

    try:
        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.starttls()
        server.login(SENDER_EMAIL, APP_PASSWORD)
        
        # [중요 수정] sendmail 대신 send_message 사용 (자동으로 인코딩 처리)
        server.send_message(msg)
        
        server.quit()
        print(f"✅ '{RECEIVER_EMAIL}'로 메일 발송 완료!")
    except Exception as e:
        print(f"❌ 메일 발송 실패: {e}")

if __name__ == "__main__":
    data, count = fetch_pmc_articles()
    if data is not None:
        send_email(data, count)
