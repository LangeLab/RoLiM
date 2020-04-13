from django.urls import path
from . import views

urlpatterns = [
	path('patterndetection/', views.index),
	path('textfile/', views.prealigned_text_file),
	path('peptidelist/', views.text_file_peptide_list),
]
