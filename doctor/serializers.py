from django.contrib.auth.models import User
from django.contrib.auth import authenticate
from rest_framework import serializers
from rest_framework_simplejwt.tokens import RefreshToken
from .models import PasswordResetOTP
import re
from rest_framework import serializers
from analysis.models import Analysis

class AnalysisHistorySerializer(serializers.ModelSerializer):

    name = serializers.SerializerMethodField()
    risk_level = serializers.SerializerMethodField()
    risk_score = serializers.SerializerMethodField()

    class Meta:
        model = Analysis
        fields = [
            "id",
            "name",
            "smiles",
            "risk_level",
            "risk_score",
            "status",
            "created_at"
        ]

    def get_name(self, obj):
        if obj.result:
            return obj.result.get("drug_overview", {}).get("name")
        return None

    def get_risk_level(self, obj):
        if obj.result:
            return obj.result.get("risk_summary", {}).get("level")
        return None

    def get_risk_score(self, obj):
        if obj.result:
            return obj.result.get("risk_summary", {}).get("risk_percentage")
        return None


# ===============================
# REGISTER SERIALIZER
# ===============================

class RegisterSerializer(serializers.ModelSerializer):
    password = serializers.CharField(write_only=True)

    class Meta:
        model = User
        fields = ['email', 'username', 'password']

    # ✅ Email validation
    def validate_email(self, value):
        if User.objects.filter(email=value).exists():
            raise serializers.ValidationError(
                "Email already registered."
            )
        return value

    # ✅ Username validation
    def validate_username(self, value):
        if len(value) < 3:
            raise serializers.ValidationError(
                "Username must be at least 3 characters long."
            )
        if not value.isalnum():
            raise serializers.ValidationError(
                "Username should contain only letters and numbers."
            )
        return value

    # ✅ Password validation (your rules)
    def validate_password(self, value):

        if len(value) < 8:
            raise serializers.ValidationError(
                "Password must be at least 8 characters long."
            )

        if not value[0].isupper():
            raise serializers.ValidationError(
                "Password must start with a capital letter."
            )

        if not re.search(r"[a-z]", value):
            raise serializers.ValidationError(
                "Password must include a lowercase letter."
            )

        if not re.search(r"[0-9]", value):
            raise serializers.ValidationError(
                "Password must include a number."
            )

        if not re.search(r"[!@#$%^&*(),.?\":{}|<>]", value):
            raise serializers.ValidationError(
                "Password must include a special character."
            )

        return value

    # ✅ Create user
    def create(self, validated_data):
        user = User.objects.create_user(
            email=validated_data['email'],
            username=validated_data['username'],
            password=validated_data['password']
        )
        return user


# ===============================
# LOGIN SERIALIZER
# ===============================

class LoginSerializer(serializers.Serializer):
    username = serializers.CharField()
    password = serializers.CharField(write_only=True)

    def validate(self, data):
        user = authenticate(
            username=data['username'],
            password=data['password']
        )

        if user is None:
            raise serializers.ValidationError("Invalid credentials")

        refresh = RefreshToken.for_user(user)

        return {
            'refresh': str(refresh),
            'access': str(refresh.access_token),
        }


# ===============================
# FORGOT PASSWORD
# ===============================

class ForgotPasswordSerializer(serializers.Serializer):
    email = serializers.EmailField()


class VerifyOTPSerializer(serializers.Serializer):
    email = serializers.EmailField()
    otp = serializers.CharField(max_length=6)


# ===============================
# RESET PASSWORD WITH VALIDATION
# ===============================

class ResetPasswordSerializer(serializers.Serializer):
    email = serializers.EmailField()
    new_password = serializers.CharField(write_only=True)

    def validate_new_password(self, value):

        if len(value) < 8:
            raise serializers.ValidationError(
                "Password must be at least 8 characters long."
            )

        if not value[0].isupper():
            raise serializers.ValidationError(
                "Password must start with a capital letter."
            )

        if not re.search(r"[a-z]", value):
            raise serializers.ValidationError(
                "Password must include a lowercase letter."
            )

        if not re.search(r"[0-9]", value):
            raise serializers.ValidationError(
                "Password must include a number."
            )

        if not re.search(r"[!@#$%^&*(),.?\":{}|<>]", value):
            raise serializers.ValidationError(
                "Password must include a special character."
            )

        return value