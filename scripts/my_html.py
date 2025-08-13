import pandas as pd
import os, sys, glob
from datetime import datetime
import base64
import pdfkit

title = 'CompleteWGS report'

def parse_report(csv_file):
    data = [line.strip().split('\t') for line in open(csv_file)]
    df = pd.DataFrame(data[1:], columns=data[0])
    return df.to_html(index=False, classes='data-table', border=1)

def parse_hlala(csv_file):
    data = [line.strip().split('\t')[:3] for line in open(csv_file)]
    df = pd.DataFrame(data[1:], columns=data[0])
    return df.to_html(index=False, classes='data-table', border=1)

def image_to_base64(image_path):
    try:
        with open(image_path, 'rb') as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except:
        return ''

def generate_html(outdir, sample):
    csv                 = glob.glob(os.path.join(outdir, '*.report'))[0]
    hlala_path          = os.path.join(outdir, 'hlala_out', sample, 'hla', 'R1_bestguess_G.txt')

    # table
    metrics     = parse_report(csv)
    hla         = parse_hlala(hlala_path)

    var_class_table           = parse_report(os.path.join(outdir, 'var_class.csv'))
    cons_type_severe_table    = parse_report(os.path.join(outdir, 'cons_type_severe.csv'))
    cons_type_all_table       = parse_report(os.path.join(outdir, 'cons_type_all.csv'))
    coding_cons_type_table    = parse_report(os.path.join(outdir, 'coding_cons_type.csv'))


    # png
    cumuplot = image_to_base64(os.path.join(outdir, 'cumulative_coverage_plot.png'))
    ideogram = image_to_base64(os.path.join(outdir, 'chromosome.png'))

    pangenie_png            = image_to_base64(os.path.join(outdir, 'pangenie_var_plot.png'))

    var_class_png           = image_to_base64(os.path.join(outdir, 'var_class.png'))
    cons_type_severe_png    = image_to_base64(os.path.join(outdir, 'cons_type_severe.png'))
    cons_type_all_png       = image_to_base64(os.path.join(outdir, 'cons_type_all.png'))
    coding_cons_type_png    = image_to_base64(os.path.join(outdir, 'coding_cons_type.png'))
    var_chrom_png           = image_to_base64(os.path.join(outdir, 'var_chrom.png'))

    # HTML模板
    html_content = f'''
    <!DOCTYPE html>
    <html lang="en">
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
                width: 50%;
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
                max-width: 50%;
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
            <h1>{title}</h1>
        </div>

        <div class="section">
            <h2>Metrics</h2>
            {metrics}
        </div>

        <div class="section">
            <h2>Phase ideogram</h2>
            <img src="data:image/png;base64,{ideogram}" alt="ideogram">
        </div>

        <div class="section">
            <h2>Phase block cumulative coverage plot </h2>
            <img src="data:image/png;base64,{cumuplot}" alt="plot">
        </div>

        <div class="section">
            <h2>Pangenie result</h2>
            <img src="data:image/png;base64,{pangenie_png}" alt="plot">
        </div>

        <div class="section">
            <h2>VEP result </h2>
            <h3>Variant classes</h3>
            <img src="data:image/png;base64,{var_class_png}" alt="plot">
            {var_class_table}

            <h3>Consequences (most severe)</h3>
            <img src="data:image/png;base64,{cons_type_severe_png}" alt="plot">
            {cons_type_severe_table}

            <h3>Consequences (all)</h3>
            <img src="data:image/png;base64,{cons_type_all_png}" alt="plot">
            {cons_type_all_table}

            <h3>Coding consequences</h3>
            <img src="data:image/png;base64,{coding_cons_type_png}" alt="plot">
            {coding_cons_type_table}

            <h3>Variants by chromosome</h3>
            <img src="data:image/png;base64,{var_chrom_png}" alt="plot">
        </div>

        <div class="section">
            <h2>HLA result</h2>
            {hla}
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
    
    # pdfkit.from_file(sample + '_report.html', sample + '_report.pdf')

if __name__ == '__main__':
    outdir = sys.argv[1] # params.outdir/report/
    samples = [ entry.name for entry in os.scandir(outdir) if entry.is_dir() ]
    for sample in samples:
        outdir1 = os.path.join(outdir, sample)
        generate_html(outdir1, sample)
