from django.urls import path
from . import views

urlpatterns = [
	path('patterndetection/', views.index),
	path('prealigned_text_file', views.prealigned_text_file),
	path('text_file_peptide_list', views.text_file_peptide_list),
]
