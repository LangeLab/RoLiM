from django.urls import path
from . import views

urlpatterns = [
	path('patterndetection/', views.index),
	path('patterndetection/textfile/', views.prealigned_text_file),
	path('patterndetection/peptidelist/', views.text_file_peptide_list),
]
