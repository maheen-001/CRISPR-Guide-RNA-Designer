# This project was made by Maheen Abbasi on Jan 2026.

# app.py
#--------------
# Create a web interface for the CRISPR gRNA designer using Flask.
# Connect user input from an HTML form to the bioinformatics logic defined in crisp_logic.py

from flask import Flask, render_template, request
import crispr_logic as cl

# Flask application instance
app = Flask(__name__)

# Main route: supports GET (display empty form) and POST (process DNA input and display results)
@app.route("/", methods=["GET", "POST"])
def index():
    # Vars passed to the HTML template
    results = None
    error = None

    # User submits the form
    if request.method == "POST":
        try:
            dna = ""

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

            guides = cl.find_grnas(dna)
            ranked = cl.rank_grnas(guides, dna)
            results = ranked[:5]
        
        except Exception as e:
            error = str(e)
    
    # Render the HTMLK page and pass results/errors to it
    return render_template("index.html", results = results, error = error)

if __name__ == "__main__":
    app.run(debug = True)