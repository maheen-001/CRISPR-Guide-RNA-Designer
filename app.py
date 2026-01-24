# This project was made by Maheen Abbasi on Jan 2026.

# TODO:
# - UI obv it is super ugly rn

# app.py
#--------------
# Create a web interface for the CRISPR gRNA designer using Flask.
# Connect user input from an HTML form to the bioinformatics logic defined in crisp_logic.py

# Main
from flask import Flask, render_template, request
import crispr_logic as cl

# CSV
import csv
from flask import Response

# PDF
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from io import BytesIO

# ------------------------------------------------------------------------------------ #
# Main Logic
# ------------------------------------------------------------------------------------ #

# Flask application instance
app = Flask(__name__)
MAX_RESULTS = 100

# Main route: supports GET (display empty form) and POST (process DNA input and display results)
@app.route("/", methods=["GET", "POST"])
def index():
    # Vars passed to the HTML template
    results = None
    error = None

    if request.method == "GET":
        results = None
        error = None
        request.form = {}

    # User submits the form
    if request.method == "POST":
        try:
            dna = ""
            enzyme = request.form.get("enzyme", "SpCas9")

            # FASTA upload
            if "fasta_file" in request.files:
                file = request.files["fasta_file"]
                if file and file.filename:
                    content = file.read().decode("utf-8")
                    dna = cl.normalize(content)

            # TEXT input
            if not dna:
                raw_input = request.form.get("dna", "")
                dna = cl.normalize(raw_input)

            # No sequence provided
            if not dna:
                raise ValueError("No DNA sequence provided.")
            
            cl.validate_dna(dna)

            guides = cl.find_all_grnas(dna, enzyme = enzyme)
            ranked = cl.rank_grnas(guides, dna)
            results = ranked[:MAX_RESULTS]
            print("Selected enzyme:", enzyme)
        
        except Exception as e:
            error = str(e)

    # Render the HTMLK page and pass results/errors to it
    return render_template("index.html", results = results, error = error)

# ------------------------------------------------------------------------------------ #
# Exporting
# ------------------------------------------------------------------------------------ #

# CSV
@app.route("/export/csv", methods = ["POST"])
def export_csv():
    dna = request.form["dna"]
    enzyme = request.form["enzyme"]

    guides = cl.find_all_grnas(dna, enzyme)
    ranked = cl.rank_grnas(guides, dna)

    # Helper
    def generate():
        yield "Guide,PAM,PAM_Class,GC,OffTargets,Score,Position,Strand,Enzyme\n"
        for g in ranked:
            yield f"{g['guide']},{g['pam']},{g['pam_class']},{g['gc']},{g['off_targets']},{g['score']},{g['position']},{g['strand']},{g['enzyme']}\n"

    return Response(
        generate(),
        mimetype = "text/csv",
        headers = {"Content-Disposition": "attachment; filename = grna_results.csv"}
    )

# PDF
@app.route("/export/pdf", methods = ["POST"])
def export_pdf():
    dna = request.form["dna"]
    enzyme = request.form["enzyme"]

    guides = cl.find_all_grnas(dna, enzyme)
    ranked = cl.rank_grnas(guides, dna)[:10]

    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize = letter)
    styles = getSampleStyleSheet()
    elements = []

    elements.append(Paragraph(
        f"<b>CRISPR gRNA Report</b><br/>Enzyme: {enzyme}<br/>Sequence length: {len(dna)}",
        styles["Normal"]
    ))

    table_data = [["Guide", "PAM", "GC", "Off-targets", "Score", "Strand"]]

    for g in ranked:
        table_data.append([
            g["guide"],
            g["pam"],
            g["gc"],
            g["off_targets"],
            g["score"],
            g["strand"]
        ])

    table = Table(table_data)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.lightgrey),
        ("GRID", (0,0), (-1,-1), 1, colors.black),
        ("FONT", (0,0), (-1,0), "Helvetica-Bold"),
        ("ALIGN", (2,1), (-1,-1), "CENTER")
    ]))

    elements.append(table)
    doc.build(elements)

    buffer.seek(0)
    return Response(
        buffer,
        mimetype = "application/pdf",
        headers={"Content-Disposition": "attachment; filename = grna_report.pdf"}
    )

if __name__ == "__main__":
    app.run(debug = True)