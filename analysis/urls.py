from django.urls import path
from .views import CreateAnalysisView, AnalysisHistoryView, AnalysisResultView, DashboardSummaryView

urlpatterns = [
    path('create/', CreateAnalysisView.as_view(), name='create-analysis'),
    path('history/', AnalysisHistoryView.as_view(), name='analysis-history'),
    path('result/<int:analysis_id>/', AnalysisResultView.as_view(), name='analysis-result'),
    path('dashboard/', DashboardSummaryView.as_view(), name='dashboard-summary'),
]