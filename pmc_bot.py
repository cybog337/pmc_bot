from Bio import Entrez
import smtplib
from email.mime.text import MIMEText
from datetime import datetime, timedelta

def check_pmc_updates(keyword, email_address):
    Entrez.email = "your_email@example.com"  # NCBI 정책상 이메일 필수
    
    # 지난 1일간의 데이터 검색
    handle = Entrez.esearch(db="pmc", term=f"{keyword} AND ("last 1 days"[dp])", retmax=10)
    record = Entrez.read(handle)
    handle.close()
    
    if int(record["Count"]) > 0:
        ids = record["IdList"]
        # 여기에 메일 발송 로직 구현 (smtplib 활용)
        print(f"새로운 논문 {len(ids)}건 발견! 메일을 발송합니다.")
    else:
        print("새로운 업데이트가 없습니다.")
