from django.urls import path
from . import views

urlpatterns = [
	path('patterndetection/', views.index),
	path('', views.splash),
]
