from django.urls import path
from .views import ProfileView, RegisterView, LoginView
from rest_framework_simplejwt.views import TokenRefreshView
from .views import ForgotPasswordView, VerifyOTPView, ResetPasswordView

urlpatterns = [
    path('register/', RegisterView.as_view(), name='register'),
    path('login/', LoginView.as_view(), name='login'),
    path('token/refresh/', TokenRefreshView.as_view(), name='token_refresh'),
    path('forgot-password/', ForgotPasswordView.as_view()),
    path('verify-otp/', VerifyOTPView.as_view()),
    path('reset-password/', ResetPasswordView.as_view()),
    path('profile/', ProfileView.as_view(), name='profile'),
]