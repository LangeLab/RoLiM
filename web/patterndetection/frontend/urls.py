from django.urls import path
from . import views

urlpatterns = [
	path('patterndetection1/', views.index),
]
