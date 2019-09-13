from tool.api import viewsets
from rest_framework.routers import DefaultRouter

router = DefaultRouter()
router.register(r'api/job', viewsets.JobViewSet, basename='job')
urlpatterns = router.urls

"""
from django.urls import path
from . import views
from tool.api import viewsets

urlpatterns = [
    path('api/backgroundformat/', views.BackgroundFormatListCreate.as_view()),
    path('api/job/', viewsets.Job),
]
"""