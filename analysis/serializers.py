from rest_framework import serializers
from .models import Analysis


class AnalysisCreateSerializer(serializers.ModelSerializer):

    class Meta:
        model = Analysis
        fields = ['input_type', 'smiles', 'sdf_file']

    def validate(self, data):
        input_type = data.get('input_type')

        if input_type == 'smiles' and not data.get('smiles'):
            raise serializers.ValidationError("SMILES is required.")

        if input_type == 'sdf' and not data.get('sdf_file'):
            raise serializers.ValidationError("SDF file is required.")

        return data


from rest_framework import serializers
from .models import Analysis


class AnalysisHistorySerializer(serializers.ModelSerializer):

    name = serializers.SerializerMethodField()
    smiles = serializers.SerializerMethodField()
    risk_level = serializers.SerializerMethodField()
    risk_score = serializers.SerializerMethodField()
    risk_color = serializers.SerializerMethodField()

    class Meta:
        model = Analysis
        fields = [
            "id",
            "name",
            "smiles",
            "risk_level",
            "risk_score",
            "risk_color",
            "status",
            "created_at",
        ]

    def get_name(self, obj):
        if obj.result:
            return obj.result.get("drug_overview", {}).get("name")
        return None

    def get_smiles(self, obj):
        return obj.smiles

    def get_risk_level(self, obj):
        if obj.result:
            return obj.result.get("risk_summary", {}).get("level")
        return None

    def get_risk_score(self, obj):
        if obj.result:
            return obj.result.get("risk_summary", {}).get("score")
        return None

    def get_risk_color(self, obj):
        if obj.result:
            return obj.result.get("risk_summary", {}).get("color")
        return None