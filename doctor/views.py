from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from rest_framework.permissions import IsAuthenticated
from django.core.mail import send_mail
from django.contrib.auth.models import User
from .serializers import (
    RegisterSerializer,
    LoginSerializer,
    ForgotPasswordSerializer,
    VerifyOTPSerializer,
    ResetPasswordSerializer
)
from .models import PasswordResetOTP
import random


# =======================
# REGISTER
# =======================
class RegisterView(APIView):

    def post(self, request):
        serializer = RegisterSerializer(data=request.data)

        if serializer.is_valid():
            serializer.save()
            return Response(
                {"message": "User registered successfully"},
                status=status.HTTP_201_CREATED
            )

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


# =======================
# LOGIN
# =======================
class LoginView(APIView):

    def post(self, request):
        serializer = LoginSerializer(data=request.data)

        if serializer.is_valid():
            return Response(serializer.validated_data, status=status.HTTP_200_OK)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


# =======================
# FORGOT PASSWORD (SEND OTP)
# =======================
class ForgotPasswordView(APIView):

    def post(self, request):
        serializer = ForgotPasswordSerializer(data=request.data)

        if serializer.is_valid():
            email = serializer.validated_data['email']

            # Safe lookup (avoids MultipleObjectsReturned)
            user = User.objects.filter(email=email).first()

            if not user:
                return Response({"error": "User not found"}, status=400)

            # Generate OTP
            otp = str(random.randint(100000, 999999))

            # Delete old OTPs
            PasswordResetOTP.objects.filter(user=user).delete()

            # Save new OTP
            PasswordResetOTP.objects.create(user=user, otp=otp)

            try:
                send_mail(
                    'Password Reset OTP',
                    f'Your OTP is {otp}',
                    'kankanampatibhumika1111@gmail.com',
                    [email],
                    fail_silently=False,
                )
            except Exception as e:
                print("Email error:", e)
                return Response(
                    {
                        "message": "OTP created but email sending failed.",
                        "otp_for_testing": otp
                    },
                    status=200
                )

            return Response({"message": "OTP sent to email"}, status=200)

        return Response(serializer.errors, status=400)


# =======================
# VERIFY OTP
# =======================
class VerifyOTPView(APIView):

    def post(self, request):
        serializer = VerifyOTPSerializer(data=request.data)

        if serializer.is_valid():
            email = serializer.validated_data['email']
            otp = serializer.validated_data['otp']

            user = User.objects.filter(email=email).first()

            if not user:
                return Response({"error": "User not found"}, status=400)

            otp_obj = PasswordResetOTP.objects.filter(user=user, otp=otp).last()

            if otp_obj:
                return Response({"message": "OTP verified"}, status=200)
            else:
                return Response({"error": "Invalid or expired OTP"}, status=400)

        return Response(serializer.errors, status=400)


# =======================
# RESET PASSWORD
# =======================
class ResetPasswordView(APIView):

    def post(self, request):
        serializer = ResetPasswordSerializer(data=request.data)

        if serializer.is_valid():
            email = serializer.validated_data['email']
            new_password = serializer.validated_data['new_password']

            user = User.objects.filter(email=email).first()

            if not user:
                return Response({"error": "User not found"}, status=400)

            user.set_password(new_password)
            user.save()

            # Delete OTP after reset
            PasswordResetOTP.objects.filter(user=user).delete()

            return Response({"message": "Password reset successful"}, status=200)

        return Response(serializer.errors, status=400)


# =======================
# PROFILE (JWT Protected)
# =======================
class ProfileView(APIView):

    permission_classes = [IsAuthenticated]

    def get(self, request):
        user = request.user

        return Response({
            "id": user.id,
            "username": user.username,
            "email": user.email,
            "organization": getattr(user, "organization", ""),
            "profile_image": request.build_absolute_uri(user.profile_image.url) if hasattr(user, "profile_image") and user.profile_image else None,
            "date_joined": user.date_joined,
        })

    def put(self, request):
        user = request.user

        # ✅ update email
        user.email = request.data.get("email", user.email)

        # ✅ update organization (if exists)
        if hasattr(user, "organization"):
            user.organization = request.data.get("organization", user.organization)

        # ✅ update profile image
        if "profile_image" in request.FILES:
            user.profile_image = request.FILES["profile_image"]

        user.save()

        return Response({"message": "Profile updated successfully"})