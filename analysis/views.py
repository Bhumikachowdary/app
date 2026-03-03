from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.permissions import IsAuthenticated
from rest_framework.parsers import MultiPartParser, FormParser, JSONParser

from .models import Analysis
from .serializers import AnalysisHistorySerializer

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

import requests
import urllib.parse
import os
import joblib
import numpy as np
import pandas as pd


# =====================================================
# BASE DIRECTORY
# =====================================================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


# =====================================================
# LOAD CTD (Protein Regulation)
# =====================================================

CTD_PATH = os.path.join(BASE_DIR, "data", "CTD_diabetes_subset.tsv")

try:
    ctd_df = pd.read_csv(CTD_PATH, sep="\t")
    print("CTD loaded")
except:
    ctd_df = None


DIABETES_GENES = [
    "INSR", "IRS1", "IRS2", "AKT1", "PRKAA1",
    "TNF", "IL6", "IL1B",
    "SOD1", "CAT", "GPX1", "NFE2L2"
]


def get_protein_regulation(chemical_name):

    if ctd_df is None or chemical_name == "Unknown":
        return {
            "upregulated_genes": [],
            "downregulated_genes": [],
            "upregulated_count": 0,
            "downregulated_count": 0
        }

    chem_data = ctd_df[
        (ctd_df["ChemicalName"].str.lower() == chemical_name.lower()) &
        (ctd_df["GeneSymbol"].isin(DIABETES_GENES))
    ]

    up, down = [], []

    for _, row in chem_data.iterrows():
        interaction = str(row["InteractionActions"]).lower()
        gene = row["GeneSymbol"]

        if "increase" in interaction:
            up.append(gene)
        elif "decrease" in interaction:
            down.append(gene)

    return {
        "upregulated_genes": list(set(up)),
        "downregulated_genes": list(set(down)),
        "upregulated_count": len(set(up)),
        "downregulated_count": len(set(down))
    }


# =====================================================
# LOAD RNA TISSUE CONSENSUS (REAL ORGAN MAPPING)
# =====================================================

RNA_PATH = os.path.join(BASE_DIR, "data", "rna_tissue_consensus.tsv")

try:
    rna_df = pd.read_csv(RNA_PATH, sep="\t")
    print("RNA loaded")
except:
    rna_df = None


def get_organ_affinity(gene_list):

    if rna_df is None or not gene_list:
        return {"primary_targets": []}

    organ_scores = {}

    # Clean RNA gene column once
    rna_df["Gene name"] = rna_df["Gene name"].astype(str).str.strip().str.upper()

    for gene in gene_list:

        gene_clean = gene.strip().upper()

        gene_data = rna_df[rna_df["Gene name"] == gene_clean]

        for _, row in gene_data.iterrows():
            tissue = row["Tissue"]
            expression = row["nTPM"]

            organ_scores[tissue] = organ_scores.get(tissue, 0) + expression

    sorted_organs = sorted(
        organ_scores.items(),
        key=lambda x: x[1],
        reverse=True
    )

    top_organs = [org[0] for org in sorted_organs[:4]]

    return {"primary_targets": top_organs}
# =====================================================
# LOAD SIDER (Side Effects)
# =====================================================

SIDER_PATH = os.path.join(BASE_DIR, "data", "meddra_all_se.tsv")

try:
    sider_df = pd.read_csv(SIDER_PATH, sep="\t", header=None)
    sider_df.columns = [
        "STITCH_ID_flat",
        "STITCH_ID_stereo",
        "UMLS_ID_label",
        "MedDRA_type",
        "UMLS_ID_medra",
        "SideEffect"
    ]
    print("SIDER loaded")
except:
    sider_df = None


def get_side_effects(pubchem_cid):

    if sider_df is None or pubchem_cid is None:
        return []

    stitch1 = f"CID{str(pubchem_cid).zfill(9)}"
    stitch2 = f"CID1{pubchem_cid}"

    matched = sider_df[
        (sider_df["STITCH_ID_flat"] == stitch1) |
        (sider_df["STITCH_ID_stereo"] == stitch1) |
        (sider_df["STITCH_ID_flat"] == stitch2) |
        (sider_df["STITCH_ID_stereo"] == stitch2)
    ]

    return matched["SideEffect"].dropna().unique().tolist()[:15]


# =====================================================
# LD50 BASED SAFE DOSE (Allometric Scaling - FDA Km Method)
# =====================================================

def calculate_safe_dose_from_ld50(animal_ld50_mg_per_kg, species="rat"):

    km_values = {
        "mouse": 3,
        "rat": 6,
        "human": 37
    }

    if species not in km_values:
        return None

    animal_km = km_values[species]
    human_km = km_values["human"]

    hed = animal_ld50_mg_per_kg * (animal_km / human_km)

    return {
        "reference_animal": species,
        "animal_ld50_mg_per_kg": animal_ld50_mg_per_kg,
        "human_equivalent_dose_mg_per_kg": round(hed, 2),
        "conversion_method": "Allometric scaling using FDA Km factors",
        "note": "Toxicological reference value for research purposes only"
    }
# =====================================================
# LOAD ML MODEL
# =====================================================

MODEL_PATH = os.path.join(BASE_DIR, "tox21_sr_are_model.pkl")
ml_model = joblib.load(MODEL_PATH)


def smiles_to_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        return np.array(fp).reshape(1, -1)
    return None


# =====================================================
# PUBCHEM LOOKUP
# =====================================================

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


def risk_color(level):
    return {
        "Low": "#22c55e",
        "Moderate": "#facc15",
        "High": "#ef4444"
    }.get(level, "#22c55e")


# =====================================================
# MAIN ANALYSIS VIEW
# =====================================================

class CreateAnalysisView(APIView):
    permission_classes = [IsAuthenticated]
    parser_classes = [MultiPartParser, FormParser, JSONParser]

    def post(self, request):

        smiles = request.data.get("smiles")

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return Response({"error": "Invalid SMILES"}, status=400)

        analysis = Analysis.objects.create(
            user=request.user,
            smiles=smiles,
            status="processing"
        )

        molecular_weight = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 2)

        pubchem = get_pubchem_data(smiles)

        fp = smiles_to_fp(smiles)
        probability = ml_model.predict_proba(fp)[0][1] if fp is not None else 0

        risk_percentage = round(probability * 100, 2)
        level = "Low" if risk_percentage < 30 else "Moderate" if risk_percentage < 70 else "High"

        protein_regulation = get_protein_regulation(pubchem["name"])

        all_genes = (
            protein_regulation["upregulated_genes"] +
            protein_regulation["downregulated_genes"]
        )

        organ_affinity = get_organ_affinity(all_genes)

        side_effects = get_side_effects(pubchem["cid"])

        # Example LD50 reference value (can later fetch from PubChem)
        reference_ld50 = 200  # mg/kg (example literature value for demonstration)

        safe_dose = calculate_safe_dose_from_ld50(
        animal_ld50_mg_per_kg=reference_ld50,
        species="rat"
)

        result_data = {
            "model": "Type 2 Diabetes Specific Human Toxicity Model",
            "analysis_id": analysis.id,
            "status": "completed",

            "drug_overview": {
                "name": pubchem["name"],
                "smiles": smiles,
                "molecular_weight": molecular_weight,
                "logP": logp,
                "pubchem_cid": pubchem["cid"]
            },

            "risk_summary": {
                "level": level,
                "score": round(risk_percentage / 10, 1),
                "risk_percentage": risk_percentage,
                "color": risk_color(level)
            },

            "protein_regulation": protein_regulation,
            "organ_affinity": organ_affinity,

            "side_effect_profile": {
                "reported_effects": side_effects,
                "effect_count": len(side_effects)
            },
              "toxicological_reference_dose": safe_dose
            }
        

        analysis.result = result_data
        analysis.status = "completed"
        analysis.save()

        return Response(result_data, status=201)
# =====================================================
# HISTORY VIEW
# =====================================================

class AnalysisHistoryView(APIView):
    permission_classes = [IsAuthenticated]

    def get(self, request):
        analyses = Analysis.objects.filter(
            user=request.user
        ).order_by('-created_at')

        serializer = AnalysisHistorySerializer(analyses, many=True)
        return Response(serializer.data)


# =====================================================
# RESULT VIEW
# =====================================================

class AnalysisResultView(APIView):
    permission_classes = [IsAuthenticated]

    def get(self, request, analysis_id):
        try:
            analysis = Analysis.objects.get(
                id=analysis_id,
                user=request.user
            )
        except Analysis.DoesNotExist:
            return Response({"detail": "Not found"}, status=404)

        return Response(analysis.result)


# =====================================================
# DASHBOARD SUMMARY
# =====================================================

class DashboardSummaryView(APIView):
    permission_classes = [IsAuthenticated]

    def get(self, request):
        analyses = Analysis.objects.filter(
            user=request.user,
            result__isnull=False
        )

        total = analyses.count()
        safe = analyses.filter(result__risk_summary__level="Low").count()
        moderate = analyses.filter(result__risk_summary__level="Moderate").count()
        high = analyses.filter(result__risk_summary__level="High").count()

        return Response({
            "total_analyses": total,
            "safe_count": safe,
            "moderate_count": moderate,
            "high_risk_count": high
        })