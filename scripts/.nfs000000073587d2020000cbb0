import pandas as pd
import os, sys
from datetime import datetime
import base64

title = 'CompleteWGS analysis report'

def read_xls(file_path):
    try:
        # 直接读取文本文件
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # 解析数据
        data = []
        headers = []
        for line in lines:
            line = line.strip()
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    data.append(parts)
                else:
                    # 如果只有一列，可能是标题行
                    headers.append(line)
        
        if data:
            # 如果有数据行，使用第一行作为表头
            df = pd.DataFrame(data[1:], columns=data[0])
        else:
            # 如果只有标题行，创建一个单列DataFrame
            df = pd.DataFrame(headers, columns=['Content'])
        
        # 转换为HTML表格
        return df.to_html(index=False, classes='data-table', border=1)
    except Exception as e:
        return f'<div style="color: red;">无法读取文件: {str(e)}</div>'

def image_to_base64(image_path):
    try:
        with open(image_path, 'rb') as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except:
        return ''

def generate_html():
    # 获取当前工作目录的绝对路径
    data_dir = os.path.join(current_dir, 'data')
    
    # 读取所有数据文件
    fqstats_html = read_xls(os.path.join(data_dir, 'fqstats.xls'))
    bamstat_html = read_xls(os.path.join(data_dir, 'demo1.bamstat.xls'))
    vcfstat_html = read_xls(os.path.join(data_dir, 'demo1.vcfstat.xls'))

    # 转换图片为Base64
    logo_base64 = image_to_base64(os.path.join(data_dir, 'logo.png'))
    base_base64 = image_to_base64(os.path.join(data_dir, 'demo1.base.png'))
    qual_base64 = image_to_base64(os.path.join(data_dir, 'demo1.qual.png'))
    insertsize_base64 = image_to_base64(os.path.join(data_dir, 'demo1.insertsize.png'))
    hist_base64 = image_to_base64(os.path.join(data_dir, 'demo1.histPlot.png'))
    cumu_base64 = image_to_base64(os.path.join(data_dir, 'demo1.cumuPlot.png'))

    # HTML模板
    html_content = f'''
    <!DOCTYPE html>
    <html lang="zh">
    <head>
        <meta charset="UTF-8">
        <title>{title}</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 0;
                padding: 20px;
                background-color: #f5f5f5;
            }}
            .header {{
                background-color: #1e88e5;
                color: white;
                padding: 20px;
                display: flex;
                align-items: center;
                margin-bottom: 30px;
            }}
            .header img {{
                height: 60px;
                margin-right: 20px;
            }}
            .header h1 {{
                margin: 0;
                font-size: 28px;
            }}
            .section {{
                background-color: white;
                border-radius: 8px;
                padding: 20px;
                margin-bottom: 20px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            .section h2 {{
                color: #1e88e5;
                border-bottom: 2px solid #1e88e5;
                padding-bottom: 10px;
                margin-top: 0;
            }}
            .subsection {{
                margin-bottom: 20px;
            }}
            .subsection h3 {{
                color: #333;
                margin-bottom: 10px;
            }}
            .data-table {{
                width: 100%;
                border-collapse: collapse;
                margin: 10px 0;
                font-size: 14px;
            }}
            .data-table th, .data-table td {{
                padding: 12px;
                border: 1px solid #ddd;
                text-align: left;
            }}
            .data-table th {{
                background-color: #1e88e5;
                color: white;
            }}
            .data-table tr:nth-child(even) {{
                background-color: #f9f9f9;
            }}
            .data-table tr:hover {{
                background-color: #f5f5f5;
            }}
            img {{
                max-width: 100%;
                height: auto;
                border-radius: 4px;
                margin: 10px 0;
            }}
            .image-container {{
                display: flex;
                flex-wrap: wrap;
                gap: 20px;
                margin: 10px 0;
            }}
            .image-container img {{
                flex: 1;
                min-width: 300px;
            }}
        </style>
    </head>
    <body>
        <div class="header">
            <img src="data:image/png;base64,{logo_base64}" alt="Logo">
            <h1>{title}</h1>
        </div>

        <div class="section">
            <h2>Metrics</h2>
                {fqstats_html}
            <div class="subsection">
                <h3>2) 碱基分布</h3>
                <img src="data:image/png;base64,{base_base64}" alt="碱基分布">
            </div>
            <div class="subsection">
                <h3>3) 碱基质量分布</h3>
                <img src="data:image/png;base64,{qual_base64}" alt="碱基质量分布">
            </div>
        </div>

        <div class="section">
            <h2>Phase block cumulative coverage plot </h2>
            <div class="subsection">
                <h3>1) 比对结果</h3>
                {bamstat_html}
            </div>
            <div class="subsection">
                <h3>2) 插入片段</h3>
                <img src="data:image/png;base64,{insertsize_base64}" alt="插入片段">
            </div>
            <div class="subsection">
                <h3>3) 深度分布</h3>
                <div class="image-container">
                    <img src="data:image/png;base64,{hist_base64}" alt="深度分布直方图">
                    <img src="data:image/png;base64,{cumu_base64}" alt="深度分布累积图">
                </div>
            </div>
        </div>

        <div class="section">
            <h2>Phase ideogram</h2>
            <div class="subsection">
                <h3>1) 变异统计</h3>
                {vcfstat_html}
            </div>
        </div>

        <div style="text-align: right; margin-top: 20px; color: #666;">
            {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        </div>
    </body>
    </html>
    '''

    # 写入HTML文件
    with open(sample + '_report.html', 'w', encoding='utf-8') as f:
        f.write(html_content)

if __name__ == '__main__':
    path = sys.argv[1]
    samples = [ entry.name for entry in os.scandir(path) if entry.is_dir() ]
    for sample in samples:
        generate_html(path, sample)
