from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.permissions import IsAuthenticated
from rest_framework.parsers import MultiPartParser, FormParser, JSONParser

from .models import Analysis
from .serializers import AnalysisHistorySerializer

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from django.http import FileResponse

import requests
import urllib.parse
import os
import joblib
import numpy as np
import pandas as pd
import random
import tempfile
from reportlab.pdfgen import canvas

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


# ===============================
# LOAD DATASETS
# ===============================

CTD_PATH = os.path.join(BASE_DIR, "data", "CTD_diabetes_subset.tsv")
RNA_PATH = os.path.join(BASE_DIR, "data", "rna_tissue_consensus.tsv")
SIDER_PATH = os.path.join(BASE_DIR, "data", "meddra_all_se.tsv")

ctd_df = pd.read_csv(CTD_PATH, sep="\t")
rna_df = pd.read_csv(RNA_PATH, sep="\t")
rna_df["Gene name"] = rna_df["Gene name"].str.upper()

sider_df = pd.read_csv(SIDER_PATH, sep="\t", header=None)
sider_df.columns = [
    "STITCH_ID_flat",
    "STITCH_ID_stereo",
    "UMLS_ID_label",
    "MedDRA_type",
    "UMLS_ID_medra",
    "SideEffect"
]


MODEL_PATH = os.path.join(BASE_DIR, "tox21_sr_are_model.pkl")
ml_model = joblib.load(MODEL_PATH)


# ===============================
# HELPERS
# ===============================

def smiles_to_fp(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        return np.array(fp).reshape(1, -1)

    return None


def get_pubchem_data(smiles):

    try:

        encoded = urllib.parse.quote(smiles)

        cid = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded}/cids/JSON"
        ).json()["IdentifierList"]["CID"][0]

        name = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title/JSON"
        ).json()["PropertyTable"]["Properties"][0]["Title"]

        return {"name": name, "cid": cid}

    except:

        return {"name": "Unknown", "cid": None}


# ===============================
# PROTEIN REGULATION
# ===============================

def get_protein_regulation(chemical_name):

    chem_data = ctd_df[
        ctd_df["ChemicalName"].str.lower() == chemical_name.lower()
    ]

    up = {}
    down = {}

    for _, row in chem_data.iterrows():

        gene = row["GeneSymbol"]
        interaction = str(row["InteractionActions"]).lower()

        pubmed = str(row.get("PubMedIDs", ""))
        evidence = pubmed.count("|") + 1
        effect_value = min(evidence * 20, 80)

        # Upregulated proteins
        if "increase" in interaction:

            if gene not in up and gene not in down:

                up[gene] = {
                    "protein": gene,
                    "effect": f"+{effect_value}%",
                    "evidence_count": evidence
                }

        # Downregulated proteins
        elif "decrease" in interaction:

            if gene not in down and gene not in up:

                down[gene] = {
                    "protein": gene,
                    "effect": f"-{effect_value}%",
                    "evidence_count": evidence
                }

    return {

        "summary": {
            "upregulated_count": len(list(up.values())[:6]),
            "downregulated_count": len(list(down.values())[:6])
        },

        "upregulated": list(up.values())[:6],

        "downregulated": list(down.values())[:6],

        "clinical_significance": {
            "title": "CTD Evidence-Based Regulation",
            "description": "Protein regulation derived from Comparative Toxicogenomics Database interactions."
        }
    }

# ===============================
# ORGAN TOXICITY
# ===============================

def get_organ_affinity(gene_list):

    organ_descriptions = {
        "Pancreas": "Central organ responsible for insulin secretion and glucose regulation.",
        "Liver": "Primary metabolic site responsible for drug detoxification and metabolism.",
        "Kidney": "Important for renal clearance and excretion of metabolites.",
        "Cardiac": "Cardiac tissue interaction may influence cardiovascular response.",
        "Neural": "Central nervous system interaction affecting neurological signaling.",
        "Gastro": "Gastrointestinal absorption and digestive tract interaction."
    }

    organ_scores = {
        "Pancreas": 0,
        "Liver": 0,
        "Kidney": 0,
        "Cardiac": 0,
        "Neural": 0,
        "Gastro": 0
    }

    organ_gene_counts = {
        "Pancreas": 0,
        "Liver": 0,
        "Kidney": 0,
        "Cardiac": 0,
        "Neural": 0,
        "Gastro": 0
    }

    for gene in gene_list:

        rows = rna_df[rna_df["Gene name"] == gene.upper()]

        for _, row in rows.iterrows():

            tissue = str(row["Tissue"]).lower()
            expression = row["nTPM"]

            if "pancreas" in tissue:
                organ_scores["Pancreas"] += expression
                organ_gene_counts["Pancreas"] += 1

            elif "liver" in tissue:
                organ_scores["Liver"] += expression
                organ_gene_counts["Liver"] += 1

            elif "kidney" in tissue:
                organ_scores["Kidney"] += expression
                organ_gene_counts["Kidney"] += 1

            elif "heart" in tissue or "cardiac" in tissue:
                organ_scores["Cardiac"] += expression
                organ_gene_counts["Cardiac"] += 1

            elif "brain" in tissue or "cerebellum" in tissue or "cortex" in tissue:
                organ_scores["Neural"] += expression
                organ_gene_counts["Neural"] += 1

            elif "stomach" in tissue or "intestine" in tissue or "colon" in tissue:
                organ_scores["Gastro"] += expression
                organ_gene_counts["Gastro"] += 1


    max_score = max(organ_scores.values()) if max(organ_scores.values()) > 0 else 1

    result = []

    for organ, score in organ_scores.items():

        gene_factor = organ_gene_counts[organ] * 2

        affinity = min(round(((score / max_score) * 70) + gene_factor + 5, 2), 100)

        result.append({
            "organ": organ,
            "affinity": affinity,
            "description": organ_descriptions.get(
                organ,
                "Biological interaction detected in this tissue."
            )
        })

    result = sorted(result, key=lambda x: x["affinity"], reverse=True)

    return result

def get_primary_targets(organ_profile):

    # sort organs by affinity
    sorted_organs = sorted(
        organ_profile,
        key=lambda x: x["affinity"],
        reverse=True
    )

    # pick top 3 organs
    top_organs = sorted_organs[:3]

    return [o["organ"] for o in top_organs]

# ===============================
# SIDE EFFECTS
# ===============================

def get_side_effects(pubchem_cid):

    if pubchem_cid is None:
        return []

    stitch = f"CID{str(pubchem_cid).zfill(9)}"

    matched = sider_df[
        sider_df["STITCH_ID_flat"] == stitch
    ]

    effects = matched["SideEffect"].dropna().unique().tolist()

    return effects[:15]


def categorize_side_effects(effects, risk_percentage):
    """Group raw side effect strings into UI-friendly categories.

    The returned structure is a list of dicts with ``category`` and
    ``effects`` keys, where ``effects`` is a list of objects containing
    ``effect``, a random ``probability`` between 5 and 20, and a static
    ``risk_level`` of "Low".
    """
    cats = {
        "Gastrointestinal": [],
        "Metabolic": [],
        "Renal/Hepatic": [],
        "Cardiovascular": [],
        "Neurological": []
    }

    for eff in effects:
        lower = eff.lower()
        if any(k in lower for k in ["nausea", "vomit", "stomach", "digest", "diarrhea", "abdominal"]):
            key = "Gastrointestinal"
        elif any(k in lower for k in ["weight", "metabolic", "glucose", "cholesterol"]):
            key = "Metabolic"
        elif any(k in lower for k in ["liver", "renal", "kidney", "hepat"]):
            key = "Renal/Hepatic"
        elif any(k in lower for k in ["heart", "cardiac", "blood pressure", "hypertension"]):
            key = "Cardiovascular"
        elif any(k in lower for k in ["headache", "dizziness", "neurolog", "seizure", "tremor"]):
            key = "Neurological"
        else:
            key = "Gastrointestinal"

        # determine a base probability influenced by the overall risk percentage
        base = 5 + int(risk_percentage / 10)
        probability = random.randint(base, base + 10)

        cats[key].append({
            "effect": eff,
            "probability": probability,
            "risk_level": "Low"
        })

    return [
        {"category": k, "effects": v}
        for k, v in cats.items()
        if v
    ]


# ===============================
# SAFE DOSE
# ===============================

def calculate_safe_dose():
    # return a simple toxicological reference dose structure
    return {
        "animal_ld50_mg_per_kg": 200,
        "human_equivalent_dose_mg_per_kg": 32.43
    }

# ===============================
# CONCENTRATION ANALYSIS
# ===============================

def generate_concentration_analysis(risk_percentage, molecular_weight, logp):

    # 🔹 Dynamic optimal dose based on risk
    if risk_percentage < 30:
        human_optimal = "500-1000 mg"
    elif risk_percentage < 70:
        human_optimal = "300-700 mg"
    else:
        human_optimal = "100-400 mg"

    # 🔹 Therapeutic index decreases with risk
    therapeutic_index = round(6 - (risk_percentage / 20), 2)

    # 🔹 Dose-response curve adjusts slightly
    base_curve = [0, 25, 60, 80, 90, 95]

    adjusted_curve = [
        max(0, min(100, int(e - (risk_percentage / 5))))
        for e in base_curve
    ]

    dose_response = [
        {"dose": i * 20, "efficacy": adjusted_curve[i]}
        for i in range(6)
    ]

    return {
        "human_optimal": human_optimal,
        "therapeutic_index": therapeutic_index,
        "dose_response": dose_response,

        # 🔹 Dynamic human model
        "human_model": {
            "min_effective": f"{int(molecular_weight * 1.5)} mg",
            "max_safe": f"{int(2000 - risk_percentage * 10)} mg",
            "optimal_range": human_optimal,
            "ld50": "Not determined",
            "therapeutic_index": therapeutic_index
        }
    }

# ===============================
# CREATE ANALYSIS
# ===============================
def generate_key_findings(organ_profile):

    findings = []

    # take top 5 organs with highest affinity
    top_organs = sorted(
        organ_profile,
        key=lambda x: x["affinity"],
        reverse=True
    )[:5]

    for organ in top_organs:

        name = organ["organ"]
        affinity = organ.get("affinity", 0)

        if name == "Pancreas":

            if affinity < 40:
                findings.append({
                    "title": "Pancreatic Safety",
                    "description": "Low pancreatic toxicity predicted."
                })
            elif affinity < 70:
                findings.append({
                    "title": "Moderate Pancreatic Interaction",
                    "description": "Possible pancreatic metabolic interaction."
                })
            else:
                findings.append({
                    "title": "Pancreatic Stress Risk",
                    "description": "High pancreatic interaction predicted."
                })

        elif name == "Liver":

            if affinity < 40:
                findings.append({
                    "title": "Low Hepatotoxicity",
                    "description": "Minimal liver toxicity risk detected."
                })
            elif affinity < 70:
                findings.append({
                    "title": "Moderate Liver Stress",
                    "description": "Possible metabolic strain on liver enzymes."
                })
            else:
                findings.append({
                    "title": "Elevated Liver Risk",
                    "description": "Potential hepatotoxic effects observed."
                })

        elif name == "Kidney":

            if affinity < 40:
                findings.append({
                    "title": "Low Renal Risk",
                    "description": "Kidney toxicity appears minimal."
                })
            else:
                findings.append({
                    "title": "Renal Toxicity Concern",
                    "description": "Possible kidney impact detected."
                })

        elif name == "Cardiac":

            if affinity < 40:
                findings.append({
                    "title": "Cardiac Safety",
                    "description": "Low cardiovascular toxicity predicted."
                })
            else:
                findings.append({
                    "title": "Cardiac Risk",
                    "description": "Potential cardiovascular toxicity detected."
                })

        elif name == "Neural":

            if affinity < 40:
                findings.append({
                    "title": "Low Neurotoxicity",
                    "description": "Minimal neurological interaction predicted."
                })
            else:
                findings.append({
                    "title": "Neurological Interaction",
                    "description": "Possible impact on neural signaling pathways."
                })

        elif name == "Gastro":

            if affinity < 40:
                findings.append({
                    "title": "Low Gastrointestinal Risk",
                    "description": "Minimal digestive system interaction predicted."
                })
            else:
                findings.append({
                    "title": "Gastrointestinal Interaction",
                    "description": "Possible digestive tract interaction detected."
                })

    return findings
class CreateAnalysisView(APIView):

    permission_classes = [IsAuthenticated]
    parser_classes = [MultiPartParser, FormParser, JSONParser]

    def post(self, request):
        # ✅ STEP 1: quick examples
        recent_analyses = Analysis.objects.filter(
            user=request.user,
            result__isnull=False
        ).order_by("-created_at")[:3]

        quick_examples = []

        for a in recent_analyses:
            quick_examples.append({
                "name": a.result.get("drug_overview", {}).get("name", "Unknown"),
                "structure": a.smiles
            })

        if not quick_examples:
            quick_examples = [
                {"name": "Metformin", "structure": "CN(C)C(=N)NC(=N)N"},
                {"name": "Gliclazide", "structure": "CC1=CC=C..."},
                {"name": "Sitagliptin", "structure": "C1CN2C..."}
            ]

        # ✅ STEP 2: HANDLE INPUT
        sdf_file = request.FILES.get("file")
        smiles = request.data.get("smiles")
        mol = None

        # 👉 CASE 1: SDF FILE
        if sdf_file:
            # ✅ case-insensitive file type validation
            if not sdf_file.name.lower().endswith(".sdf"):
                return Response({"error": "Only .sdf files are allowed"}, status=400)

            if sdf_file.size == 0:
                return Response({"error": "Empty file uploaded"}, status=400)

            tmp_path = None
            try:
                # Save uploaded file temporarily
                with tempfile.NamedTemporaryFile(delete=False, suffix='.sdf') as tmp_file:
                    for chunk in sdf_file.chunks():
                        tmp_file.write(chunk)
                    tmp_path = tmp_file.name
                
                # Read SDF from temporary file path
                suppl = Chem.SDMolSupplier(tmp_path, removeHs=False)
                mols = [m for m in suppl if m is not None]

                if not mols:
                    return Response({"error": "Invalid SDF file content"}, status=400)

                mol = mols[0]
                smiles = Chem.MolToSmiles(mol)

            except Exception as e:
                return Response({"error": f"Invalid SDF file: {str(e)}"}, status=400)

            finally:
                # ✅ Cleanup: Always delete temp file if it exists
                if tmp_path and os.path.exists(tmp_path):
                    os.remove(tmp_path)

        # 👉 CASE 2: SMILES INPUT
        elif smiles:
            mol = Chem.MolFromSmiles(smiles)

            if not mol:
                return Response({"error": "Invalid SMILES"}, status=400)

        # 👉 CASE 3: NOTHING PROVIDED
        else:
            return Response({"error": "No input provided"}, status=400)

        # ✅ CONTINUE ANALYSIS
        analysis = Analysis.objects.create(
            user=request.user,
            smiles=smiles,
            status="processing"
        )

        molecular_weight = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 2)

        pubchem = get_pubchem_data(smiles)

        fp = smiles_to_fp(smiles)
        if fp is None:
            return Response({"error": "Unable to generate molecular fingerprint"}, status=400)

        probability = ml_model.predict_proba(fp)[0][1]

        risk_percentage = round(probability * 100, 2)
        level = "Low" if risk_percentage < 30 else "Moderate" if risk_percentage < 70 else "High"

        protein_data = get_protein_regulation(pubchem["name"])

        genes = [g["protein"] for g in protein_data["upregulated"]] + \
                [g["protein"] for g in protein_data["downregulated"]]

        organ_profile = get_organ_affinity(genes)
        primary_targets = get_primary_targets(organ_profile)
        key_findings = generate_key_findings(organ_profile)

        concentration_data = generate_concentration_analysis(
            risk_percentage,
            molecular_weight,
            logp
        )

        side_effects = get_side_effects(pubchem["cid"])

        if not side_effects:
            side_effects = [
                "Nausea", "Diarrhea", "Headache",
                "Dizziness", "Weight change",
                "Blood pressure change", "Renal discomfort"
            ]

        categorized = categorize_side_effects(side_effects, risk_percentage)

        result_data = {
            "quick_examples": quick_examples,
            "concentration_analysis": concentration_data,
            "primary_target_organs": primary_targets,
            "key_findings": key_findings,
            "ai_confidence": round((1 - probability) * 100, 1),
            "analysis_id": analysis.id,

            "drug_overview": {
                "name": pubchem["name"],
                "smiles": smiles,
                "molecular_weight": molecular_weight,
                "logP": logp
            },

            "risk_summary": {
                "level": level,
                "risk_percentage": risk_percentage
            },

            "protein_analysis": protein_data,
            "organ_toxicity_profile": organ_profile,

            "side_effect_profile": {
                "overall_safety": {
                    "score": 100 - int(risk_percentage),
                    "classification": level.upper()
                },
                "categories": categorized
            },

            "toxicological_reference_dose": calculate_safe_dose()
        }

        analysis.result = result_data
        analysis.status = "completed"
        analysis.risk_level = level
        analysis.save()

        return Response(result_data)


# ===============================
# ANALYSIS HISTORY
# ===============================
class AnalysisHistoryView(APIView):

    permission_classes = [IsAuthenticated]

    def get(self, request):
        """Retrieve all analyses for the current user"""
        analyses = Analysis.objects.filter(
            user=request.user
        ).order_by('-created_at')
        
        serializer = AnalysisHistorySerializer(analyses, many=True)
        return Response(serializer.data)


# ===============================
# ANALYSIS RESULT
# ===============================
class AnalysisResultView(APIView):

    permission_classes = [IsAuthenticated]

    def get(self, request, analysis_id):
        """Retrieve a specific analysis result"""
        try:
            analysis = Analysis.objects.get(id=analysis_id, user=request.user)
        except Analysis.DoesNotExist:
            return Response(
                {"error": "Analysis not found"},
                status=404
            )

        serializer = AnalysisHistorySerializer(analysis)
        return Response(serializer.data)


# ===============================
# DASHBOARD SUMMARY
# ===============================
class DashboardSummaryView(APIView):
    permission_classes = [IsAuthenticated]

    def get(self, request):
        try:
            analyses = Analysis.objects.filter(
                user=request.user,
                result__isnull=False
            )

            total = analyses.count()
            safe = analyses.filter(risk_level="Low").count()
            moderate = analyses.filter(risk_level="Moderate").count()
            high = analyses.filter(risk_level="High").count()

            recent = analyses.order_by("-created_at")[:3]

            recent_data = []
            for a in recent:
                recent_data.append({
                    "drug_name": a.result.get("drug_overview", {}).get("name", "Unknown"),
                    "smiles": a.smiles,
                    "risk_level": a.result.get("risk_summary", {}).get("level", "Unknown"),
                    "risk_score": a.result.get("risk_summary", {}).get("risk_percentage", 0),
                    "created_at": a.created_at
                })

            return Response({
                "doctor_name": request.user.username,

                "statistics": {
                    "total_analyses": total,
                    "safe_count": safe,
                    "moderate_count": moderate,
                    "high_risk_count": high
                },

                "recent_analyses": recent_data
            })

        except Exception as e:
            return Response({"error": str(e)}, status=500)


# ===============================
# DOWNLOAD REPORT
# ===============================
class DownloadReportView(APIView):

    permission_classes = [IsAuthenticated]

    def get(self, request, analysis_id):
        """Generate and download a PDF report for an analysis"""
        try:
            analysis = Analysis.objects.get(id=analysis_id, user=request.user)
        except Analysis.DoesNotExist:
            return Response(
                {"error": "Analysis not found"},
                status=404
            )

        if not analysis.result:
            return Response(
                {"error": "Analysis result not available"},
                status=400
            )

        # Create PDF report
        buffer = tempfile.NamedTemporaryFile(delete=False, suffix='.pdf')
        pdf_path = buffer.name
        buffer.close()

        try:
            c = canvas.Canvas(pdf_path)
            c.setFont("Helvetica", 12)
            
            # Header
            c.drawString(50, 750, f"Drug Toxicity Analysis Report")
            c.drawString(50, 730, f"Analysis ID: {analysis.id}")
            
            # Drug Overview
            drug_info = analysis.result.get('drug_overview', {})
            c.drawString(50, 700, f"Drug Name: {drug_info.get('name', 'N/A')}")
            c.drawString(50, 680, f"SMILES: {drug_info.get('smiles', 'N/A')[:50]}")
            c.drawString(50, 660, f"Molecular Weight: {drug_info.get('molecular_weight', 'N/A')} g/mol")
            
            # Risk Summary
            risk = analysis.result.get('risk_summary', {})
            c.drawString(50, 630, f"Risk Level: {risk.get('level', 'N/A')}")
            c.drawString(50, 610, f"Risk Percentage: {risk.get('risk_percentage', 'N/A')}%")
            
            c.save()

            return FileResponse(
                open(pdf_path, 'rb'),
                as_attachment=True,
                filename=f"analysis_{analysis.id}_report.pdf"
            )

        except Exception as e:
            return Response(
                {"error": f"Error generating report: {str(e)}"},
                status=500
            )
        finally:
            if os.path.exists(pdf_path):
                os.remove(pdf_path)