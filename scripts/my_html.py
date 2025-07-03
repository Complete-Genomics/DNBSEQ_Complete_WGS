import pandas as pd
import os, sys, glob
from datetime import datetime
import base64

title = 'CompleteWGS report'

def read_xls(csv):
    try:
        data = [line.strip.split(',') for line in open(csv)]
        df = pd.DataFrame(data[1:], columns=data[0])
        return df.to_html(index=False, classes='data-table', border=1)
    except Exception as e:
        return f'<div style="color: red;">无法读取文件: {str(e)}</div>'

def image_to_base64(image_path):
    try:
        with open(image_path, 'rb') as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except:
        return ''

def generate_html(outdir, sample):
    csv = glob.glob(os.path.join(outdir, '*.report'))[0]

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
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>{title}</title>
    </head>
    <body>
        <div class="header">
            <h1>{title}</h1>
        </div>

        <div class="section">
            <h2>Metrics</h2>
                {csv}
        </div>

        <div class="section">
            <h2>Phase block cumulative coverage plot </h2>
            <div class="subsection">
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
    outdir = sys.argv[1] # params.outdir/report/
    samples = [ entry.name for entry in os.scandir(path) if entry.is_dir() ]
    for sample in samples:
        outdir1 = os.path.join(outdir, sample)
        generate_html(outdir1, sample)
