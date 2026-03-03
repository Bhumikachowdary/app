from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from .serializers import RegisterSerializer, LoginSerializer
from django.core.mail import send_mail
from django.contrib.auth.models import User
from .models import PasswordResetOTP
from .serializers import ForgotPasswordSerializer, VerifyOTPSerializer, ResetPasswordSerializer
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
import random
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.permissions import IsAuthenticated


class RegisterView(APIView):

    def post(self, request):
        serializer = RegisterSerializer(data=request.data)

        if serializer.is_valid():
            serializer.save()
            return Response({
                "message": "User registered successfully"
            }, status=status.HTTP_201_CREATED)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class LoginView(APIView):

    def post(self, request):
        serializer = LoginSerializer(data=request.data)

        if serializer.is_valid():
            return Response(serializer.validated_data, status=status.HTTP_200_OK)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
class ForgotPasswordView(APIView):

    def post(self, request):
        serializer = ForgotPasswordSerializer(data=request.data)
        if serializer.is_valid():

            email = serializer.validated_data['email']

            try:
                user = User.objects.get(email=email)
            except User.DoesNotExist:
                return Response({"error": "User not found"}, status=404)

            otp = str(random.randint(100000, 999999))

            PasswordResetOTP.objects.create(user=user, otp=otp)

            send_mail(
                'Password Reset OTP',
                f'Your OTP is {otp}',
                'yourgmail@gmail.com',
                [email],
                fail_silently=False,
            )

            return Response({"message": "OTP sent to email"}, status=200)

        return Response(serializer.errors, status=400)
class VerifyOTPView(APIView):

    def post(self, request):
        serializer = VerifyOTPSerializer(data=request.data)

        if serializer.is_valid():
            email = serializer.validated_data['email']
            otp = serializer.validated_data['otp']

            try:
                user = User.objects.get(email=email)
                otp_obj = PasswordResetOTP.objects.filter(user=user, otp=otp).last()
            except:
                return Response({"error": "Invalid OTP"}, status=400)

            if otp_obj:
                return Response({"message": "OTP verified"}, status=200)
            else:
                return Response({"error": "Invalid OTP"}, status=400)

        return Response(serializer.errors, status=400)
class ResetPasswordView(APIView):

    def post(self, request):
        serializer = ResetPasswordSerializer(data=request.data)

        if serializer.is_valid():
            email = serializer.validated_data['email']
            new_password = serializer.validated_data['new_password']

            try:
                user = User.objects.get(email=email)
                user.set_password(new_password)
                user.save()
                return Response({"message": "Password reset successful"}, status=200)
            except:
                return Response({"error": "User not found"}, status=404)

        return Response(serializer.errors, status=400)



class ProfileView(APIView):
    permission_classes = [IsAuthenticated]

    def get(self, request):
        user = request.user

        return Response({
            "id": user.id,
            "username": user.username,
            "email": user.email,
            "date_joined": user.date_joined,
        })