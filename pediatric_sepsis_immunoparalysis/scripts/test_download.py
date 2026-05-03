#!/usr/bin/env python3
"""测试GEO数据下载"""
import requests
from pathlib import Path

# 测试GSE26378下载
gse_id = "GSE26378"
url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file"
print(f"Testing download URL: {url}")

dest_path = Path(f"./长期计划/儿童脓毒症免疫瘫痪研究/data/raw/{gse_id}_series_matrix.txt.gz")

try:
    print("Starting download...")
    response = requests.get(url, stream=True, timeout=60)
    print(f"Status: {response.status_code}")
    print(f"Content-Length: {response.headers.get('content-length', 'N/A')}")
    
    if response.status_code == 200:
        with open(dest_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        print(f"Downloaded to: {dest_path}")
        print(f"File size: {dest_path.stat().st_size} bytes")
    else:
        print(f"Download failed: {response.status_code}")
except Exception as e:
    print(f"Error: {e}")
