import pandas as pd
import os, sys, glob
from datetime import datetime
import base64
# import pdfkit

title = 'cWGS Report'

SKIPPED = set(filter(None, sys.argv[2].split(','))) if len(sys.argv) > 2 else set()

def fail_box(label, reason='process failed or produced no output; see .nextflow.log / task work dir'):
    return f'<div class="failTitle">FAILED: {label} ({reason})</div>'

def skipped_box(label):
    return f'<div class="noDataTitle">SKIPPED: {label} (disabled by parameters)</div>'

def parse_report(csv_file, label='result'):
    try:
        data = [line.strip().split('\t') for line in open(csv_file)]
        df = pd.DataFrame(data[1:], columns=data[0])
        return df.to_html(index=False, classes='data-table', border=1)
    except Exception:
        return fail_box(label)

def parse_specimmune(result_dir, category):
    if category in SKIPPED:
        return skipped_box(category)
    # marker written by the specimmune process when the tool itself failed
    marker = os.path.join(result_dir, category, 'FAILED')
    if os.path.exists(marker):
        reason = open(marker).read().strip() or 'tool exited with error'
        return fail_box(category, reason)
    try:
        patterns = [
            os.path.join(result_dir, category, '**', '*.hap.alleles.txt'),
            os.path.join(result_dir, category, '**', '*.final.type.result.formatted.txt'),
        ]
        matches = []
        for pattern in patterns:
            matches.extend(glob.glob(pattern, recursive=True))
        result_file = next((path for path in matches if os.path.getsize(path) > 0), None)
        if not result_file:
            raise FileNotFoundError(category)
        data = [line.rstrip('\n').split('\t') for line in open(result_file)]
        df = pd.DataFrame(data[1:], columns=data[0])
        return df.to_html(index=False, classes='data-table', border=1)
    except Exception:
        return fail_box(category)

def image_to_base64(image_path):
    try:
        with open(image_path, 'rb') as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except Exception:
        return ''

def img_tag(image_path, label):
    b64 = image_to_base64(image_path)
    if b64:
        return f'<img src="data:image/png;base64,{b64}" alt="{label}">'
    # a *.FAILED marker next to the expected png carries the reason (written by the plot processes)
    marker = os.path.splitext(image_path)[0] + '.FAILED'
    if os.path.exists(marker):
        try:
            reason = open(marker).read().strip().replace('<', '&lt;') or 'process reported failure'
        except Exception:
            reason = 'process reported failure'
        return fail_box(label, reason)
    return fail_box(label)

def report_image_path(sample_dir, filename):
    report_dir = os.path.dirname(sample_dir.rstrip(os.sep))
    candidates = [
        os.path.join(sample_dir, filename),
        os.path.join(report_dir, filename),
    ]
    for image_path in candidates:
        if os.path.exists(image_path):
            return image_path
    return candidates[0]

def generate_html(outdir, sample):
    csv_matches         = glob.glob(os.path.join(outdir, sample + '.lariat.dv.report'))
    csv                 = csv_matches[0] if csv_matches else None
    specimmune_dir      = os.path.join(outdir, 'specimmune_out')

    # table
    metrics     = parse_report(csv) if csv else fail_box('Metrics (stLFRQC)', 'report file not found')
    hla         = parse_specimmune(specimmune_dir, 'HLA')

    var_class_table           = parse_report(os.path.join(outdir, 'var_class.csv'), 'VEP var_class table')
    cons_type_severe_table    = parse_report(os.path.join(outdir, 'cons_type_severe.csv'), 'VEP cons_type_severe table')
    cons_type_all_table       = parse_report(os.path.join(outdir, 'cons_type_all.csv'), 'VEP cons_type_all table')
    coding_cons_type_table    = parse_report(os.path.join(outdir, 'coding_cons_type.csv'), 'VEP coding_cons_type table')


    # png
    cumuplot = img_tag(report_image_path(outdir, 'cumulative_coverage_plot.png'), 'cumulative_coverage_plot')
    ideogram = img_tag(report_image_path(outdir, 'chromosome_sv.png'), 'chromosome_sv')

    pangenie_png = img_tag(report_image_path(outdir, 'pangenie_var_plot.png'), 'pangenie_var_plot')

    var_class_png = img_tag(report_image_path(outdir, 'var_class.png'), 'var_class')
    cons_type_severe_png = img_tag(report_image_path(outdir, 'cons_type_severe.png'), 'cons_type_severe')
    cons_type_all_png = img_tag(report_image_path(outdir, 'cons_type_all.png'), 'cons_type_all')
    coding_cons_type_png = img_tag(report_image_path(outdir, 'coding_cons_type.png'), 'coding_cons_type')
    var_chrom_png = img_tag(report_image_path(outdir, 'var_chrom.png'), 'var_chrom')

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
            .failTitle {{
                color: #b71c1c;
                background-color: #ffebee;
                border: 1px solid #b71c1c;
                border-radius: 4px;
                padding: 10px;
                margin: 10px 0;
                font-weight: bold;
            }}
            .noDataTitle {{
                color: #666;
                padding: 10px 0;
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
            <h2>Phase ideogram and >10k SV</h2>
            {ideogram}
        </div>

        <div class="section">
            <h2>Phase block cumulative coverage plot </h2>
            {cumuplot}
        </div>

        <div class="section">
            <h2>Pangenie result</h2>
            {pangenie_png}
        </div>

        <div class="section">
            <h2>VEP result </h2>
            <h3>Variant classes</h3>
            {var_class_png}
            {var_class_table}

            <h3>Consequences (most severe)</h3>
            {cons_type_severe_png}
            {cons_type_severe_table}

            <h3>Consequences (all)</h3>
            {cons_type_all_png}
            {cons_type_all_table}

            <h3>Coding consequences</h3>
            {coding_cons_type_png}
            {coding_cons_type_table}

            <h3>Variants by chromosome</h3>
            {var_chrom_png}
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
        try:
            generate_html(outdir1, sample)
        except Exception as e:
            print(sample, e)
            continue
