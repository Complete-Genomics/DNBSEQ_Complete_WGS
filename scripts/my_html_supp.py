import argparse
import base64
import glob
import os
from datetime import datetime

import pandas as pd

try:
    from html import escape as html_escape
except ImportError:
    try:
        from cgi import escape as html_escape
    except ImportError:
        html_escape = None


TITLE = 'cWGS Report'


def escape_html(value):
    value = str(value)
    if html_escape:
        return html_escape(value)
    return (
        value.replace('&', '&amp;')
        .replace('<', '&lt;')
        .replace('>', '&gt;')
        .replace('"', '&quot;')
    )


def missing_html(label, expected_path):
    path = escape_html(expected_path)
    return (
        '<div class="noDataTitle">{label} was not generated for this sample.'
        '<br><span class="missingPath">Expected file: {path}</span></div>'
    ).format(
        label=label,
        path=path
    )


def normalize_report_dir(outdir):
    outdir = os.path.abspath(outdir)
    if os.path.basename(outdir.rstrip(os.sep)) == 'report':
        return outdir
    report_dir = os.path.join(outdir, 'report')
    if os.path.isdir(report_dir):
        return report_dir
    return outdir


def parse_report(csv_file, label='Result', expected_path=None):
    expected_path = expected_path or csv_file or 'unknown'
    if not csv_file or not os.path.exists(csv_file):
        return missing_html(label, expected_path)
    try:
        data = [line.strip().split('\t') for line in open(csv_file)]
        df = pd.DataFrame(data[1:], columns=data[0])
        return df.to_html(index=False, classes='data-table', border=1)
    except Exception:
        return missing_html(label, expected_path)


def parse_hlala(csv_file):
    if not os.path.exists(csv_file):
        return missing_html('HLA result', csv_file)
    try:
        data = [line.strip().split('\t')[:3] for line in open(csv_file)]
        df = pd.DataFrame(data[1:], columns=data[0])
        return df.to_html(index=False, classes='data-table', border=1)
    except Exception:
        return missing_html('HLA result', csv_file)


def image_to_base64(image_path):
    if not os.path.exists(image_path):
        return ''
    try:
        with open(image_path, 'rb') as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except Exception:
        return ''


def image_html(sample_dir, filename, alt, label):
    image_path = os.path.join(sample_dir, filename)
    image_data = image_to_base64(image_path)
    if not image_data:
        return missing_html(label, image_path)
    return '<img src="data:image/png;base64,{image_data}" alt="{alt}">'.format(
        image_data=image_data,
        alt=alt
    )


def find_metrics_report(sample_dir, sample):
    matches = glob.glob(os.path.join(sample_dir, sample + '.lariat.dv.report'))
    if matches:
        return matches[0]
    matches = glob.glob(os.path.join(sample_dir, '*.lariat.dv.report'))
    return matches[0] if matches else None


def sample_dirs(report_dir, selected_sample=None):
    if selected_sample:
        sample_dir = os.path.join(report_dir, selected_sample)
        if os.path.isdir(sample_dir):
            return [(selected_sample, sample_dir)]
        raise IOError('Sample report directory does not exist: {0}'.format(sample_dir))

    dirs = []
    for name in os.listdir(report_dir):
        sample_dir = os.path.join(report_dir, name)
        if os.path.isdir(sample_dir):
            dirs.append((name, sample_dir))
    return sorted(dirs)


def generate_html(report_dir, sample, sample_dir):
    metrics_report = find_metrics_report(sample_dir, sample)
    metrics_expected = os.path.join(sample_dir, sample + '.lariat.dv.report')
    hlala_path = os.path.join(sample_dir, 'hlala_out', sample, 'hla', 'R1_bestguess_G.txt')

    metrics = parse_report(metrics_report, 'stLFRQC result', metrics_expected)
    hla = parse_hlala(hlala_path)

    var_class_table = parse_report(os.path.join(sample_dir, 'var_class.csv'), 'Variant class table')
    cons_type_severe_table = parse_report(os.path.join(sample_dir, 'cons_type_severe.csv'), 'Most severe consequence table')
    cons_type_all_table = parse_report(os.path.join(sample_dir, 'cons_type_all.csv'), 'All consequence table')
    coding_cons_type_table = parse_report(os.path.join(sample_dir, 'coding_cons_type.csv'), 'Coding consequence table')

    ideogram = image_html(sample_dir, 'chromosome_sv.png', 'ideogram', 'Phase ideogram')
    cumuplot = image_html(sample_dir, 'cumulative_coverage_plot.png', 'plot', 'Phase block cumulative coverage plot')
    pangenie_png = image_html(sample_dir, 'pangenie_var_plot.png', 'plot', 'Pangenie result')

    var_class_png = image_html(sample_dir, 'var_class.png', 'plot', 'Variant class plot')
    cons_type_severe_png = image_html(sample_dir, 'cons_type_severe.png', 'plot', 'Most severe consequence plot')
    cons_type_all_png = image_html(sample_dir, 'cons_type_all.png', 'plot', 'All consequence plot')
    coding_cons_type_png = image_html(sample_dir, 'coding_cons_type.png', 'plot', 'Coding consequence plot')
    var_chrom_png = image_html(sample_dir, 'var_chrom.png', 'plot', 'Variants by chromosome plot')

    html_content = '''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>{TITLE}</title>
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
            .noDataTitle {{
                color: #666;
                padding: 12px 0;
            }}
            .missingPath {{
                display: inline-block;
                margin-top: 6px;
                font-family: monospace;
                font-size: 12px;
                color: #333;
                word-break: break-all;
            }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>{TITLE}</h1>
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
            <h2>Phase block cumulative coverage plot</h2>
            {cumuplot}
        </div>

        <div class="section">
            <h2>Pangenie result</h2>
            {pangenie_png}
        </div>

        <div class="section">
            <h2>VEP result</h2>
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
            {timestamp}
        </div>
    </body>
    </html>
    '''.format(
        TITLE=TITLE,
        metrics=metrics,
        ideogram=ideogram,
        cumuplot=cumuplot,
        pangenie_png=pangenie_png,
        var_class_png=var_class_png,
        var_class_table=var_class_table,
        cons_type_severe_png=cons_type_severe_png,
        cons_type_severe_table=cons_type_severe_table,
        cons_type_all_png=cons_type_all_png,
        cons_type_all_table=cons_type_all_table,
        coding_cons_type_png=coding_cons_type_png,
        coding_cons_type_table=coding_cons_type_table,
        var_chrom_png=var_chrom_png,
        hla=hla,
        timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    )

    output_html = os.path.join(report_dir, sample + '_report.html')
    with open(output_html, 'w', encoding='utf-8') as f:
        f.write(html_content)
    return output_html


def main():
    parser = argparse.ArgumentParser(
        description='Regenerate cWGS HTML reports from per-sample report output directories.'
    )
    parser.add_argument('outdir', help='Pipeline outdir or the outdir/report directory')
    parser.add_argument('--sample', help='Regenerate one sample only')
    args = parser.parse_args()

    report_dir = normalize_report_dir(args.outdir)
    for sample, sample_dir in sample_dirs(report_dir, args.sample):
        output_html = generate_html(report_dir, sample, sample_dir)
        print(output_html)


if __name__ == '__main__':
    main()
